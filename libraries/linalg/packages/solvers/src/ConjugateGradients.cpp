// --- Linalg Includes ---
#include "packages/solvers/impl/ConjugateGradients_impl.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"


namespace cie::linalg {


template class ConjugateGradients<DefaultSpace<float,tags::Serial>>;
template class ConjugateGradients<DefaultSpace<float,tags::SMP>>;
template class ConjugateGradients<DefaultSpace<double,tags::Serial>>;
template class ConjugateGradients<DefaultSpace<double,tags::SMP>>;


} // namespace cie::linalg
