#ifndef CIE_MPI_IMPL_BASE_HPP
#define CIE_MPI_IMPL_BASE_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/RankID.hpp"
#include "packages/concurrency/inc/MessageTag.hpp"
#include "packages/concurrency/inc/mpi_fwd.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::mpi {


class MPI;


class MPIImpl
{
public:
    using Container = std::string;

    // {begin, size}
    using Out = std::pair<Ptr<const char>, const int>;

    // {begin, size}
    using In = std::pair<Ptr<char>, const int>;

public:
    virtual ~MPIImpl() {}

    /// @name Generic Interface
    /// @{

    virtual void send(Out message,
                      RankID destination,
                      MessageTag tag) = 0;

    virtual void receive(In message,
                         RankID source,
                         MessageTag tag) = 0;

    virtual void sendAndReceive(Out send,
                                RankID sendTo,
                                In receive,
                                RankID receiveFrom,
                                MessageTag tag) = 0;

    virtual void broadcast(In message,
                           RankID source,
                           MessageTag tag) = 0;

    /// @}
    /// @name Synchronization
    /// @{

    virtual void barrier() = 0;

    /// @}
    /// @name Queries
    /// @{

    virtual Size size() const = 0;

    virtual RankID getRankID() const noexcept = 0;

    virtual RankID getMasterRankID() const = 0;

    virtual Ptr<MPI_Comm> getComm() = 0;

    virtual bool isDistributed() const noexcept = 0;

    /// @}

protected:
    MPIImpl() noexcept = default;

    MPIImpl(MPIImpl&&) noexcept = default;

    MPIImpl& operator=(MPIImpl&&) noexcept = default;

private:
    MPIImpl(const MPIImpl&) = delete;

    MPIImpl& operator=(const MPIImpl&) = delete;
}; // class MPI::Impl


} // namespace cie::mpi


#endif
