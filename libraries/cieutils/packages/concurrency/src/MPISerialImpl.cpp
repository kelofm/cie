// --- Utility Includes ---
#include "packages/concurrency/inc/MPISerialImpl.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <cstring>


namespace cie::mpi {


MPISerialImpl::MPISerialImpl() noexcept
    : MPIImpl()
{
}


void MPISerialImpl::send(Out, RankID, MessageTag)
{
    CIE_THROW(Exception, "Attempted send operation on serial MPI")
}


void MPISerialImpl::receive(In, RankID, MessageTag)
{
    CIE_THROW(Exception, "Attempted receive operation on serial MPI")
}


void MPISerialImpl::sendAndReceive(Out send,
                                   RankID sendTo,
                                   In receive,
                                   RankID receiveFrom,
                                   MessageTag tag)
{
    CIE_CHECK(
        !sendTo && !receiveFrom,
        "send-receive operation from " << receiveFrom
        << " to " << sendTo << " with tag " << tag
        << " on serial MPI"
    )
    CIE_OUT_OF_RANGE_CHECK(send.second == receive.second)
    std::memcpy(receive.first, send.first, receive.second);
}


void MPISerialImpl::broadcast(In,
                              RankID source,
                              MessageTag tag)
{
    CIE_CHECK(
        !source,
        "broadcast operation from " << source
        << " with tag " << tag << " on serial MPI"
    )
}


void MPISerialImpl::barrier()
{
}


Size MPISerialImpl::size() const
{
    return 1u;
}


RankID MPISerialImpl::getRankID() const noexcept
{
    return RankID(0);
}


RankID MPISerialImpl::getMasterRankID() const
{
    return RankID(0);
}


Ptr<MPI_Comm> MPISerialImpl::getComm()
{
    return nullptr;
}


bool MPISerialImpl::isDistributed() const noexcept
{
    return false;
}


} // namespace cie::mpi
