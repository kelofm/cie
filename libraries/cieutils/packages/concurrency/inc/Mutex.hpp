#ifndef CIE_CONCURRENCY_MUTEX_HPP
#define CIE_CONCURRENCY_MUTEX_HPP

// --- Utility Includes ---
#include "packages/types/inc/tags.hpp"

// --- STL Includes ---
#include <mutex>


namespace cie::mp {


/// @ingroup cieutils
template <class TTag>
class Mutex
{}; // class Mutex


/// @ingroup cieutils
template <>
class Mutex<tags::SMP>
{
public:
    using Wrapped = std::mutex;

    using ScopedLock = std::scoped_lock<Wrapped>;

public:
    [[nodiscard]] ScopedLock makeScopedLock() {return ScopedLock(_mutex);}

    void lock() {_mutex.lock();}

    void unlock() {_mutex.unlock();}

private:
    Wrapped _mutex;
}; // class Mutex<SMP>


/// @ingroup cieutils
template <>
class Mutex<tags::Serial>
{
public:
    struct Wrapped
    {
        Wrapped& operator=(Wrapped&& r_rhs) = default;

        Wrapped& operator=(Ref<const Wrapped> r_rhs) = delete;

        void lock() {}

        void unlock() {}

        bool try_lock() {return true;}
    }; // struct Mutex

    struct ScopedLock
    {
        ScopedLock& operator=(ScopedLock&& r_rhs) = default;

        ScopedLock& operator=(Ref<const ScopedLock> r_rhs) = delete;
    }; // struct ScopedLock

public:
    [[nodiscard]] ScopedLock makeScopedLock() {return ScopedLock();}

    void lock() {}

    void unlock() {}
}; // class Mutex<Serial>


} // namespace cie::mp


#endif
