#pragma once

// --- Utility Includes ---
#include "packages/concurrency/inc/mpi_fwd.hpp"
#include "packages/io/inc/Serializer.hpp"
#include "packages/concurrency/inc/RankID.hpp"
#include "packages/concurrency/inc/MessageTag.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::mpi {


class MPIImpl;


class MPI {
private:
    struct MPIKey {};

public:
    /** @brief Construct the root MPI interface.
     *  @param pComm: Ptr to the MPI_Comm to use, or nullptr.
     *  @param key: private-access key ensuring that this constructor cannot be called by outsiders.
     *  @details If @a pComm is not null, this constructor initializes MPI and duplicates
     *           the passed communicator.
     */
    MPI(Ptr<MPI_Comm> pComm, MPIKey key);

    ~MPI();

    /// @name Value Return Interface
    /// @{

    template <class T>
    T receive(RankID source, MessageTag tag = 0);

    /// @}
    /// @name Ref Interface
    /// @{

    template <TriviallySerializable T>
    void send(
        Ref<const T> rMessage,
        RankID destination,
        MessageTag tag = 0);

    template <TriviallySerializable T>
    void receive(
        Ref<T> rMessage,
        RankID source,
        MessageTag tag = 0);

    template <TriviallySerializable T>
    void sendAndReceive(
        Ref<const T> r_send,
        RankID sendTo,
        Ref<T> rReceive,
        RankID receiveFrom,
        MessageTag tag = 0);

    template <TriviallySerializable T>
    void broadcast(
        Ref<T> rMessage,
        RankID source,
        MessageTag tag = 0);

    /// @}
    /// @name Synchronization
    /// @{

    void barrier();

    /// @}
    /// @name Queries
    /// @{

    Size size() const;

    RankID getRankID() const noexcept;

    RankID getMasterRankID() const;

    /// @brief Get the static MPI interface.
    static Ref<MPI> getRoot();

    Ptr<MPI_Comm> getComm();

    bool isMaster() const;

    bool isDistributed() const noexcept;

    bool isRoot() const noexcept;

    /// @}

private:
    MPI(RightRef<std::unique_ptr<MPIImpl>> rp_impl);

    MPI(MPI&&) = default;

    MPI(const MPI&&) = delete;

    MPI& operator=(MPI&&) = delete;

    MPI& operator=(const MPI&) = delete;

private:
    friend class MPIImpl;

    friend class MPISerialImpl;

    friend class MPIDistributedImpl;

    friend class MPISingleton;

    std::unique_ptr<MPIImpl> _p_impl;
}; // class MPI


} // namespace cie::mpi

#include "packages/concurrency/impl/MPI_impl.hpp"
