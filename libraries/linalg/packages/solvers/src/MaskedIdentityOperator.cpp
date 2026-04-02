// --- Linalg Includes ---
#include "packages/solvers/inc/MaskedIdentityOperator.hpp"
#include "packages/solvers/impl/MaskedIdentityOperator_impl.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"


namespace cie::linalg {


#define CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(T,TTag)                                                                    \
    template <>                                                                                                             \
    void MaskedIdentityOperator<DefaultSpace<T,TTag>>::product(                                                             \
        typename DefaultSpace<T,TTag>::Value inScale,                                                                       \
        typename DefaultSpace<T,TTag>::ConstVectorView in,                                                                  \
        typename DefaultSpace<T,TTag>::Value outScale,                                                                      \
        typename DefaultSpace<T,TTag>::VectorView out) {                                                                    \
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
    template class MaskedIdentityOperator<DefaultSpace<T,TTag>>;


CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(float, tags::Serial)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(float, tags::SMP)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(double, tags::Serial)
CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR(double, tags::SMP)


#undef CIE_INSTANTIATE_MASKED_IDENTITY_OPERATOR


} // namespace cie::linalg
