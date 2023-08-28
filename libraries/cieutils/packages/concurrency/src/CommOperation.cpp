// --- Utility Includes ---
#include "packages/concurrency/inc/CommOperation.hpp"


namespace cie::mpi {


CommOperation::CommOperation(RightRef<std::shared_ptr<MPI>> r_mpi)
    : _p_mpi(std::move(r_mpi))
{
}


Ref<const MPI> CommOperation::getMPI() const noexcept
{
    return *_p_mpi;
}


Ref<MPI> CommOperation::getMPI() noexcept
{
    return *_p_mpi;
}


std::shared_ptr<const MPI> CommOperation::shareMPI() const
{
    return _p_mpi;
}


std::shared_ptr<MPI> CommOperation::shareMPI()
{
    return _p_mpi;
}


} // namespace cie::mpi
