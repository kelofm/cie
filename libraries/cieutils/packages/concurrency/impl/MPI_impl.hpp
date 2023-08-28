#ifndef CIE_MPI_IMPL_HPP
#define CIE_MPI_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/MPI.hpp"
#include "packages/concurrency/inc/MPIImpl.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::mpi {


template <class T>
inline T MPI::receive(RankID source, MessageTag tag)
{
    T object;
    this->receive(object, source, tag);
}


template <concepts::TriviallySerializable T>
inline void MPI::send(Ref<const T> r_message, RankID destination, MessageTag tag)
{
    _p_impl->send(MPIImpl::Out {reinterpret_cast<const char*>(&r_message), sizeof(T)},
                  destination,
                  tag);
}


template <concepts::TriviallySerializable T>
inline void MPI::receive(Ref<T> r_message, RankID source, MessageTag tag)
{
    _p_impl->receive(MPIImpl::In {reinterpret_cast<char*>(&r_message), sizeof(T)},
                     source,
                     tag);
}


template <concepts::TriviallySerializable T>
inline void MPI::sendAndReceive(Ref<const T> r_send,
                                RankID sendTo,
                                Ref<T> r_receive,
                                RankID receiveFrom,
                                MessageTag tag)
{
    _p_impl->sendAndReceive(MPIImpl::Out {reinterpret_cast<const char*>(&r_send), sizeof(T)},
                            sendTo,
                            MPIImpl::In {reinterpret_cast<char*>(&r_receive), sizeof(T)},
                            receiveFrom,
                            tag);
}


template <concepts::TriviallySerializable T>
inline void MPI::broadcast(Ref<T> r_message,
                           RankID source,
                           MessageTag tag)
{
    _p_impl->broadcast(MPIImpl::In {reinterpret_cast<char*>(&r_message), sizeof(T)},
                       source,
                       tag);
}


inline bool MPI::isDistributed() const noexcept
{
    return _p_impl->isDistributed();
}


} // namespace cie::mpi


#endif
