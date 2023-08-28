// --- Utility Includes ---
#include "packages/concurrency/inc/MPI.hpp"
#include "packages/concurrency/inc/MPISingleton.hpp"
#include "packages/concurrency/inc/MPISerialImpl.hpp"
#include "packages/concurrency/inc/MPIDistributedImpl.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <sstream>
#include <mutex>


namespace cie::mpi {


MPI::MPI(Ptr<MPI_Comm> p_comm, MPIKey)
{
    #ifdef CIEUTILS_COMPILE_MPI
    if (p_comm) {
        _p_impl.reset(new MPIDistributedImpl(p_comm));
    } else {
        _p_impl.reset(new MPISerialImpl());
    }
    #else
    CIE_CHECK(!p_comm, "CiE is compiled without MPI support")
    _p_impl.reset(new MPISerialImpl());
    #endif
}


MPI::MPI(RightRef<std::unique_ptr<MPIImpl>> rp_impl)
    : _p_impl(std::move(rp_impl))
{
}


MPI::~MPI()
{
}


Size MPI::size() const
{
    return _p_impl->size();
}


RankID MPI::getRankID() const noexcept
{
    return _p_impl->getRankID();
}


RankID MPI::getMasterRankID() const
{
    return _p_impl->getMasterRankID();
}


Ptr<MPI_Comm> MPI::getComm()
{
    return _p_impl->getComm();
}


Ref<MPI> MPI::getRoot()
{
    return MPISingleton::get();
}


bool MPI::isMaster() const
{
    return this->getRankID() == this->getMasterRankID();
}


} // namespace cie::mpi
