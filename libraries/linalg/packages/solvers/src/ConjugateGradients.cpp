// --- Linalg Includes ---
#include "packages/solvers/impl/ConjugateGradients_impl.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"


namespace cie::linalg {


template class ConjugateGradients<DefaultSpace<float,tags::Serial>>;
template class ConjugateGradients<DefaultSpace<float,tags::SMP>>;
template class ConjugateGradients<DefaultSpace<double,tags::Serial>>;
template class ConjugateGradients<DefaultSpace<double,tags::SMP>>;


#ifdef CIE_ENABLE_SYCL

template class ConjugateGradients<SYCLSpace<float>>;
template class ConjugateGradients<SYCLSpace<double>>;

#endif


} // namespace cie::linalg
