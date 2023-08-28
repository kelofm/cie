// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadUtils.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <unordered_set>
#include <optional>
#include <mutex>
#include <thread>


namespace cie::mp {


namespace {
class ThreadIDGenerator
{
public:
    ThreadIDGenerator()
        : _threadIDs()
    {
        std::mutex mutex;
        const std::size_t numberOfThreads = std::thread::hardware_concurrency();
        for (std::size_t i_thread=0; i_thread<numberOfThreads; ++i_thread) {
            std::scoped_lock<std::mutex> lock(mutex);
            _threadIDs.insert(std::this_thread::get_id());
        } // for i_thread in range(numberOfThreads)
        _threadIDs.insert(std::this_thread::get_id()); // <== main thread
    }

    const std::unordered_set<std::thread::id>& get() const
    {
        return _threadIDs;
    }

private:
    std::unordered_set<std::thread::id> _threadIDs;
}; // class ThreadIDGenerator

class ThreadIDGeneratorSingleton
{
public:
    static const std::unordered_set<std::thread::id>& get()
    {
        if (!_threadIDGenerator.has_value()) {
            _threadIDGenerator = std::make_optional<ThreadIDGenerator>();
        }
        return _threadIDGenerator.value().get();
    }

private:
    static std::optional<ThreadIDGenerator> _threadIDGenerator;
}; // class ThreadIDGeneratorSingleton


std::optional<ThreadIDGenerator> ThreadIDGeneratorSingleton::_threadIDGenerator;
} // unnamed namespace


const std::unordered_set<std::thread::id>& getThreadIDs()
{
    return ThreadIDGeneratorSingleton::get();
}


} // namespace cie::mp
