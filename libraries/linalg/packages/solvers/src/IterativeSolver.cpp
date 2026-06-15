// --- Linalg Includes ---
#include "packages/solvers/inc/IterativeSolver.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TS>
IterativeSolver<TS>::IterativeSolver(Ref<const Status> rConfiguration)
    : LoggedOperator<TS>() {
        this->configure(rConfiguration);
}


template <LinalgSpaceLike TS>
typename IterativeSolver<TS>::Status IterativeSolver<TS>::configuration() const {
    return this->streamLogger().linearSystemConfiguration();
}


template <LinalgSpaceLike TS>
typename IterativeSolver<TS>::Status IterativeSolver<TS>::status() const {
    return this->streamLogger().linearSystemStatus();
}


template <LinalgSpaceLike TS>
void IterativeSolver<TS>::configure(Ref<const Status> rConfiguration) {
    this->streamLogger().configureLinearSystem(rConfiguration);
}


template <LinalgSpaceLike TS>
void IterativeSolver<TS>::updateStatus(Ref<const Status> rStatus) {
    this->streamLogger().submitLinearSystemStatus(rStatus);
}


#define CIE_INSTANTIATE_ITERATIVE_SOLVER(T)         \
    template class IterativeSolver<DefaultSpace<T>>;


CIE_INSTANTIATE_ITERATIVE_SOLVER(float)
CIE_INSTANTIATE_ITERATIVE_SOLVER(double)


#undef CIE_INSTANTIATE_ITERATIVE_SOLVER


#ifdef CIE_ENABLE_SYCL

template class IterativeSolver<SYCLSpace<float>>;
template class IterativeSolver<SYCLSpace<double>>;

#endif


} // namespace cie::linalg
