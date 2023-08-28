// --- Utility Includes ---
#include "packages/concurrency/inc/MPIDistributedImpl.hpp"
#include "packages/concurrency/inc/MPISingleton.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/logging.hpp"


namespace cie::mpi {


#ifdef CIEUTILS_COMPILE_MPI


MPIDistributedImpl::MPIDistributedImpl(Ptr<MPI_Comm> p_comm)
    : MPIImpl(),
      _comm()
{
    int isInitialized;
    MPI_Initialized(&isInitialized);

    if (!isInitialized) [[likely]] {
        const auto threadingRequirement = MPI_THREAD_SINGLE;

        // Initialize the root MPI interface
        {
            int argc = 0;
            char** argv = nullptr;
            int threadingSupport;
            const auto error = MPI_Init_thread(&argc, &argv, threadingRequirement, &threadingSupport);
            CIE_CHECK(!error, "MPI initialization failed on rank " << this->getRankID() << " with error code " << error)
        }

        // Duplicate and set the incoming communicator
        {
            CIE_CHECK(p_comm, "MPIDistributedImpl got nullptr for MPI_Comm*")
            const auto error = MPI_Comm_dup(*p_comm, &_comm);
            CIE_CHECK(!error, "Failed to duplicate MPI communicator and returned code " << error)
        }

        // Check the threading support MPI provided
        {
            int threadingSupport;
            MPI_Query_thread(&threadingSupport);
            if (threadingSupport < threadingRequirement)
                CIE_THROW(Exception, "MPI was initialized with threading support level " << threadingSupport << " but CiE requires " << threadingRequirement)
        }
    } else {
        CIE_THROW(Exception, "Rank " << this->getRankID() << " attempted to reinitialize MPI")
    } // isInitialized
}


MPIDistributedImpl::~MPIDistributedImpl()
{
    if (MPI_Comm_free(&_comm)) std::cout << "Failed to free MPI communicator on rank " << this->getRankID() << std::endl;

    // MPISingleton can only be uninitialized when an MPI
    // destructor is called, if the MPI to be destroyed is
    // the static interface => this is the last MPI interface
    // and MPI_Finalize must be called.
    if (!MPISingleton::isInitialized()) {
        int isFinalized;
        MPI_Finalized(&isFinalized);

        if (!isFinalized)
        {
            const auto error = MPI_Finalize();
            if (error)
                std::cout << "MPI finalization failed on rank " << this->getRankID() << std::endl;
        }
    }
}


void MPIDistributedImpl::send(Out message, RankID destination, MessageTag tag)
{
    CIE_PROFILE_SCOPE
    this->checkBlockingMessage(destination);
    const auto error = MPI_Send(message.first, message.second, MPI_CHAR, destination.getInt(), tag, _comm);
    if (error)
        CIE_THROW(Exception, "MPI send operation from " << this->getRankID() << " to " << destination << " with tag " << tag << " returned code " << error)
}


void MPIDistributedImpl::receive(In message, RankID source, MessageTag tag)
{
    CIE_PROFILE_SCOPE
    this->checkBlockingMessage(source);
    MPI_Status status;
    const auto error = MPI_Recv(message.first, message.second, MPI_CHAR, source.getInt(), tag, _comm, &status);
    if (error)
        CIE_THROW(Exception, "MPI receive operation from " << this->getRankID() << " to " << source << " with tag " << tag << " returned code " << error)
}


void MPIDistributedImpl::sendAndReceive(Out send, RankID sendTo, In receive, RankID receiveFrom, MessageTag tag)
{
    CIE_PROFILE_SCOPE
    MPI_Status status;
    const auto error = MPI_Sendrecv(send.first,         // send buffer
                                    send.second,        // send size
                                    MPI_CHAR,           // send type
                                    sendTo,             // destination
                                    tag,                // send tag
                                    receive.first,      // receive buffer
                                    receive.second,     // receive size
                                    MPI_CHAR,           // receive type
                                    receiveFrom,        // source
                                    tag,                // receive tag
                                    _comm,              // communicator
                                    &status);           // status

    CIE_CHECK(
        !error,
        "MPI send-receive operation from " << receiveFrom
        << " to " << sendTo << " with tag " << tag
        << " returned code " << error
    )
}


void MPIDistributedImpl::broadcast(In message, RankID source, MessageTag tag)
{
    CIE_PROFILE_SCOPE
    const auto error = MPI_Bcast(message.first, message.second, MPI_CHAR, source, _comm);
    CIE_CHECK(
        !error,
        "MPI broadcast operation from " << source
        << " on " << this->getRankID()
        << " with tag " << tag
        << " returned code " << error
    )
}


void MPIDistributedImpl::barrier()
{
    const auto error = MPI_Barrier(_comm);
    if (error)
        CIE_THROW(Exception, "MPI barrier operation on " << this->getRankID() << " returned code " << error)
}


Size MPIDistributedImpl::size() const
{
    int size;
    MPI_Comm_size(_comm, &size);
    return size;
}


RankID MPIDistributedImpl::getRankID() const noexcept
{
    int rankID;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankID);
    return RankID(rankID);
}


RankID MPIDistributedImpl::getMasterRankID() const
{
    return RankID(0);
}


Ptr<MPI_Comm> MPIDistributedImpl::getComm()
{
    return &_comm;
}


bool MPIDistributedImpl::isDistributed() const noexcept
{
    return true;
}


inline void MPIDistributedImpl::checkBlockingMessage(RankID targetRankID)
{
    // A blocking message's recipient/sender can't be itself
    CIE_CHECK(targetRankID != this->getRankID(), "Rank " << targetRankID << " attempted blocking communication with itself")
}


#endif // CIEUTILS_COMPILE_MPI


} // namespace cie::mpi
