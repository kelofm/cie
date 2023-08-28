#ifndef CIE_UTILS_MPI_SERIAL_IMPL_HPP
#define CIE_UTILS_MPI_SERIAL_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/MPIImpl.hpp"


namespace cie::mpi {


class MPISerialImpl final : public MPIImpl
{
private:
    using Buffer = std::stringstream;

public:
    /// @name Generic Interface
    /// @{

    void send(Out message, RankID destination, MessageTag tag) override;

    void receive(In message, RankID source, MessageTag tag) override;

    void sendAndReceive(Out send,
                        RankID sendTo,
                        In receive,
                        RankID receiveFrom,
                        MessageTag tag) override;

    void broadcast(In message, RankID source, MessageTag tag) override;

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
    /// @brief Root constructor.
    MPISerialImpl() noexcept;

private:
    friend class MPI;
}; // class MPISerialImpl


} // namespace cie::mpi

#endif
