#pragma once

// --- Utility Includes ---
#include "packages/concurrency/inc/MPI.hpp"
#include "packages/concurrency/inc/MPIImpl.hpp"


namespace cie::mpi {


template <class T>
inline T MPI::receive(RankID source, MessageTag tag) {
    T object;
    this->receive(
        object,
        source,
        tag);
}


template <TriviallySerializable T>
inline void MPI::send(
    Ref<const T> rMessage,
    RankID destination,
    MessageTag tag) {
        _p_impl->send(
            MPIImpl::Out {
                reinterpret_cast<const char*>(&rMessage),
                sizeof(T)},
            destination,
            tag);
}


template <TriviallySerializable T>
inline void MPI::receive(
    Ref<T> rMessage,
    RankID source,
    MessageTag tag) {
        _p_impl->receive(
            MPIImpl::In {
                reinterpret_cast<char*>(&rMessage),
                sizeof(T)},
            source,
            tag);
}


template <TriviallySerializable T>
inline void MPI::sendAndReceive(
    Ref<const T> rSend,
    RankID sendTo,
    Ref<T> rReceive,
    RankID receiveFrom,
    MessageTag tag) {
    _p_impl->sendAndReceive(
        MPIImpl::Out {
            reinterpret_cast<const char*>(&rSend),
            sizeof(T)},
        sendTo,
        MPIImpl::In {
            reinterpret_cast<char*>(&rReceive),
            sizeof(T)},
        receiveFrom,
        tag);
}


template <TriviallySerializable T>
inline void MPI::broadcast(
    Ref<T> rMessage,
    RankID source,
    MessageTag tag) {
        _p_impl->broadcast(
            MPIImpl::In {
                reinterpret_cast<char*>(&rMessage),
                sizeof(T)},
            source,
            tag);
}


inline bool MPI::isDistributed() const noexcept {
    return _p_impl->isDistributed();
}


} // namespace cie::mpi
