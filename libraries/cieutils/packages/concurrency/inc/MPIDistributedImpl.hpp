#ifndef CIE_UTILS_MPI_DISTRIBUTED_IMPL_HPP
#define CIE_UTILS_MPI_DISTRIBUTED_IMPL_HPP

// --- External Includes ---
#ifdef CIEUTILS_COMPILE_MPI
    #include "mpi.h"
#endif

// --- Utility Includes ---
#include "packages/concurrency/inc/mpi_fwd.hpp"
#include "packages/concurrency/inc/MPIImpl.hpp"


namespace cie::mpi {


#ifdef CIEUTILS_COMPILE_MPI


class MPIDistributedImpl final : public MPIImpl
{
public:
    ~MPIDistributedImpl() override;

    /// @name generic Interface
    /// @{

    void send(Out message, RankID destination, MessageTag tag) override;

    void receive(In message, RankID source, MessageTag tag) override;

    void sendAndReceive(Out send,
                        RankID sendTo,
                        In receive,
                        RankID receiveFrom,
                        MessageTag tag) override;

    void broadcast(In message,
                   RankID source,
                   MessageTag tag) override;

    /// @}
    /// @name Synchronization
    /// @{

    void barrier() override;

    /// @}
    /// @name Queries
    /// @{

    Size size() const override;

    RankID getRankID() const noexcept override;

    RankID getMasterRankID() const override;

    Ptr<MPI_Comm> getComm() override;

    bool isDistributed() const noexcept override;

    /// @}

private:
    MPIDistributedImpl(Ptr<MPI_Comm>);

    MPIDistributedImpl(MPIDistributedImpl&& r_rhs) = delete;

    MPIDistributedImpl(const MPIDistributedImpl& r_rhs) = delete;

    void checkBlockingMessage(RankID targetRankID);

private:
    MPI_Comm _comm;

private:
    friend class MPI;
}; // class MPIDistributedImpl


#else
struct MPIDistributedImpl {};
#endif


} // namespace cie::mpi


#endif
