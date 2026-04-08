// --- Linalg Includes ---
#include "packages/solvers/inc/MaskedIdentityOperator.hpp"
#include "packages/solvers/impl/MaskedIdentityOperator_impl.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"


namespace cie::linalg {


#define CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(T,TMask)                                                                   \
    template <>                                                                                                             \
    void MaskedIdentityOperator<DefaultSpace<T>,DefaultSpace<TMask>>::product(                                              \
        typename DefaultSpace<T>::Value inScale,                                                                            \
        typename DefaultSpace<T>::ConstVectorView in,                                                                       \
        typename DefaultSpace<T>::Value outScale,                                                                           \
        typename DefaultSpace<T>::VectorView out) {                                                                         \
            const auto kernel = [in, out, this] (std::size_t iRow, const auto& op) -> void {                                \
                op(                                                                                                         \
                    out[iRow],                                                                                              \
                    _mask[iRow] < _threshold ? in[iRow] : static_cast<T>(0));                                               \
            }; /*kernel*/                                                                                                   \
                                                                                                                            \
            const std::size_t rowCount = _pSpace->size(in);                                                                 \
            const auto job = [this, &kernel, rowCount] (const auto& op) -> void {                                           \
                OptionalRef<mp::ThreadPoolBase> maybeThreads = _pSpace->getThreads();                                       \
                if (maybeThreads.has_value()) {                                                                             \
                    mp::ParallelFor<std::size_t>(maybeThreads.value()).execute(                                             \
                        rowCount,                                                                                           \
                        [&op, &kernel] (std::size_t iRow) -> void {kernel(iRow, op);});                                     \
                } else {                                                                                                    \
                    for (std::size_t iRow=0ul; iRow<rowCount; ++iRow)                                                       \
                        kernel(iRow, op);                                                                                   \
                }                                                                                                           \
            }; /*job*/                                                                                                      \
                                                                                                                            \
            if (inScale == static_cast<T>(1)) {                                                                             \
                if (outScale == static_cast<T>(1)) {                                                                        \
                    job([] (Ref<T> rLeft, T right) -> void {rLeft += right;});                                              \
                } /*if outScale == 1*/ else if (outScale == static_cast<T>(-1)) {                                           \
                    job([] (Ref<T> rLeft, T right) -> void {rLeft -= right;});                                              \
                } /*if outScale == -1*/ else {                                                                              \
                    job([outScale] (Ref<T> rLeft, T right) -> void {rLeft += outScale * right;});                           \
                }                                                                                                           \
            } /*if inScale == 1*/ else if (inScale == static_cast<T>(0)) {                                                  \
                if (outScale == static_cast<T>(1)) {                                                                        \
                    job([] (Ref<T> rLeft, T right) -> void {rLeft = right;});                                               \
                } /*if outScale == 1*/ else if (outScale == static_cast<T>(-1)) {                                           \
                    job([] (Ref<T> rLeft, T right) -> void {rLeft = -right;});                                              \
                } /*if outScale == -1*/ else {                                                                              \
                    job([outScale] (Ref<T> rLeft, T right) -> void {rLeft = outScale * right;});                            \
                }                                                                                                           \
            } /*if inScale == 0*/ else if (inScale == static_cast<T>(-1)) {                                                 \
                if (outScale == static_cast<T>(1)) {                                                                        \
                    job([] (Ref<T> rLeft, T right) -> void {rLeft = right - rLeft;});                                       \
                } /*if outScale == 1*/ else if (outScale == static_cast<T>(-1)) {                                           \
                    job([] (Ref<T> rLeft, T right) -> void {rLeft = right - rLeft;});                                       \
                } /*if outScale == -1*/ else {                                                                              \
                    job([outScale] (Ref<T> rLeft, T right) -> void {rLeft = outScale * right - rLeft;});                    \
                }                                                                                                           \
            } /*if inScale == -1*/ else {                                                                                   \
                if (outScale == static_cast<T>(1)) {                                                                        \
                    job([inScale] (Ref<T> rLeft, T right) -> void {rLeft = inScale * rLeft + right;});                      \
                } /*if outScale == 1*/ else if (outScale == static_cast<T>(-1)) {                                           \
                    job([inScale] (Ref<T> rLeft, T right) -> void {rLeft = inScale * rLeft - right;});                      \
                } /*if outScale == -1*/ else {                                                                              \
                    job([inScale, outScale] (Ref<T> rLeft, T right) -> void {rLeft = inScale * rLeft + outScale * right;}); \
                }                                                                                                           \
            } /*else*/                                                                                                      \
    }                                                                                                                       \
    template class MaskedIdentityOperator<DefaultSpace<T>,DefaultSpace<TMask>>;


CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(float, std::uint16_t)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(float, int)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(float, std::size_t)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(double, std::uint16_t)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(double, int)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(double, std::size_t)


#undef CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR


#ifdef CIE_ENABLE_SYCL

#define CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(T,TMask)                                                           \
    template <>                                                                                                     \
    void MaskedIdentityOperator<SYCLSpace<T>,SYCLSpace<TMask>>::product(                                            \
        typename SYCLSpace<T>::Value inScale,                                                                       \
        typename SYCLSpace<T>::ConstVectorView in,                                                                  \
        typename SYCLSpace<T>::Value outScale,                                                                      \
        typename SYCLSpace<T>::VectorView out) {                                                                    \
            CIE_CHECK(in.size() == out.size(), "")                                                                  \
            CIE_BEGIN_EXCEPTION_TRACING                                                                             \
                const auto job = [&] (const auto& op) -> void {                                                     \
                    Ptr<const T> pInBegin = in.get();                                                               \
                    Ptr<T> pOutBegin = out.get();                                                                   \
                    Ptr<const TMask> pMaskBegin = _mask.get();                                                      \
                    const TMask threshold = _threshold;                                                             \
                    _pSpace->getQueue()->parallel_for(                                                              \
                        sycl::range<1>(in.size()),                                                                  \
                        [=] (sycl::item<1> it) -> void {                                                            \
                            const std::size_t iRow = it.get_linear_id();                                            \
                            if (pMaskBegin[iRow] < threshold)                                                       \
                                pOutBegin[iRow] = op(pOutBegin[iRow], pInBegin[iRow]);                              \
                            else                                                                                    \
                                pOutBegin[iRow] = 0;                                                                \
                        }).wait_and_throw();                                                                        \
                }; /*job*/                                                                                          \
                                                                                                                    \
                if (inScale == static_cast<T>(1)) {                                                                 \
                if (outScale == static_cast<T>(1)) {                                                                \
                    job([] (T, T right) -> T {return right;});                                                      \
                } /*if outScale == 1*/ else if (outScale == static_cast<T>(-1)) {                                   \
                    job([] (T, T right) -> T {return right;});                                                      \
                } /*if outScale == -1*/ else {                                                                      \
                    job([outScale] (T, T right) -> T {return outScale * right;});                                   \
                }                                                                                                   \
            } /*if inScale == 1*/ else if (inScale == static_cast<T>(0)) {                                          \
                if (outScale == static_cast<T>(1)) {                                                                \
                    job([] (T, T right) -> T {return right;});                                                      \
                } /*if outScale == 1*/ else if (outScale == static_cast<T>(-1)) {                                   \
                    job([] (T, T right) -> T {return -right;});                                                     \
                } /*if outScale == -1*/ else {                                                                      \
                    job([outScale] (T, T right) -> T {return outScale * right;});                                   \
                }                                                                                                   \
            } /*if inScale == 0*/ else if (inScale == static_cast<T>(-1)) {                                         \
                if (outScale == static_cast<T>(1)) {                                                                \
                    job([] (T left, T right) -> T {return right - left;});                                          \
                } /*if outScale == 1*/ else if (outScale == static_cast<T>(-1)) {                                   \
                    job([] (T left, T right) -> T {return right - left;});                                          \
                } /*if outScale == -1*/ else {                                                                      \
                    job([outScale] (T left, T right) -> T {return outScale * right - left;});                       \
                }                                                                                                   \
            } /*if inScale == -1*/ else {                                                                           \
                if (outScale == static_cast<T>(1)) {                                                                \
                    job([inScale] (T left, T right) -> T {return inScale * left + right;});                         \
                } /*if outScale == 1*/ else if (outScale == static_cast<T>(-1)) {                                   \
                    job([inScale] (T left, T right) -> T {return inScale * left - right;});                         \
                } /*if outScale == -1*/ else {                                                                      \
                    job([inScale, outScale] (T left, T right) -> T {return inScale * left + outScale * right;});    \
                }                                                                                                   \
            } /*else*/                                                                                              \
            CIE_END_EXCEPTION_TRACING                                                                               \
        }                                                                                                           \
        template class MaskedIdentityOperator<SYCLSpace<T>,SYCLSpace<TMask>>;

CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(float, std::uint16_t)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(float, int)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(float, std::size_t)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(double, std::uint16_t)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(double, int)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(double, std::size_t)

#undef CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR

#endif


} // namespace cie::linalg
