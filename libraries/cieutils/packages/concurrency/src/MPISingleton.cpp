// --- External Includes ---
#ifdef CIEUTILS_COMPILE_MPI
    #include "mpi.h"
#endif

// --- Utility Includes ---
#include "packages/concurrency/inc/MPISingleton.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie::mpi {


std::optional<std::shared_ptr<MPI>> MPISingleton::_p_mpi;


MPISingleton::Mutex MPISingleton::_mutex;


Ref<MPI> MPISingleton::get()
{
    CIE_BEGIN_EXCEPTION_TRACING
    constexpr bool enforceDistributed = false;
    constexpr Ptr<MPI_Comm> p_comm = nullptr;
    MPISingleton::ensureInitialized(enforceDistributed, p_comm);
    return *_p_mpi.value();
    CIE_END_EXCEPTION_TRACING
}


Ref<MPI> MPISingleton::get(Ptr<MPI_Comm> p_comm)
{
    CIE_BEGIN_EXCEPTION_TRACING
    constexpr bool enforceDistributed = true;
    MPISingleton::ensureInitialized(enforceDistributed, p_comm);
    return *_p_mpi.value();
    CIE_END_EXCEPTION_TRACING
}


std::shared_ptr<MPI> MPISingleton::share()
{
    CIE_BEGIN_EXCEPTION_TRACING
    constexpr bool enforceDistributed = false;
    constexpr Ptr<MPI_Comm> p_comm = nullptr;
    MPISingleton::ensureInitialized(enforceDistributed, p_comm);
    return _p_mpi.value();
    CIE_END_EXCEPTION_TRACING
}


std::shared_ptr<MPI> MPISingleton::share(Ptr<MPI_Comm> p_comm)
{
    CIE_BEGIN_EXCEPTION_TRACING
    constexpr bool enforceDistributed = true;
    MPISingleton::ensureInitialized(enforceDistributed, p_comm);
    return _p_mpi.value();
    CIE_END_EXCEPTION_TRACING
}


void MPISingleton::ensureInitialized(bool req, Ptr<MPI_Comm> p_comm)
{
    Lock lock(_mutex);

    if (!MPISingleton::isInitialized()) {
        // MPI is not initialized => check what to do
        // 0) if p_comm is non-null, a distributed communicator should be constructed with it
        // 1) otherwise duplicate MPI_COMM_WORLD (if compiled with MPI support)
        //    or construct a serial interface (if compiled without MPI support).
        MPI::MPIKey key;
        if (p_comm) {
            CIE_BEGIN_EXCEPTION_TRACING
            _p_mpi.emplace(std::make_shared<MPI>(p_comm, key));
            CIE_END_EXCEPTION_TRACING
        } else {
            CIE_BEGIN_EXCEPTION_TRACING
            #ifdef CIEUTILS_COMPILE_MPI
            auto zaWarudo = MPI_COMM_WORLD;
            p_comm = &zaWarudo;
            #else
            p_comm = nullptr;
            #endif
            _p_mpi.emplace(std::make_shared<MPI>(p_comm, key));
            CIE_END_EXCEPTION_TRACING
        }
    } else if (req) {
        CIE_CHECK(
            p_comm == _p_mpi.value()->getComm(),
            "MPISingleton is initialized with communicator at " << _p_mpi.value()->getComm()
                << " but is requrested with communicator at " << p_comm
        ) // CIE_CHECK
    } // if req
}


bool MPISingleton::isInitialized()
{
    return _p_mpi.has_value();
}


bool MPISingleton::isDistributed()
{
    return _p_mpi.has_value() && _p_mpi.value()->isDistributed();
}


void MPISingleton::reset()
{
    CIE_BEGIN_EXCEPTION_TRACING
    Lock lock(_mutex);
    if (MPISingleton::isInitialized()) {
        const auto useCount = _p_mpi->use_count();
        CIE_CHECK(
            useCount == 1,
            "The root MPI interface could not be destroyed because its ownership is shared with "
            << useCount - 1 << " other(s)"
        )
        _p_mpi.reset();
    }
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::mpi
