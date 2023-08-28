// --- Utility Includes ---
#include "packages/logging/inc/Profiler.hpp"
#include "packages/output/inc/FileManager.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/io/inc/json.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <thread>
#include <vector>
#include <tuple>
#include <fstream>
#include <sstream>


namespace cie::utils {


namespace {
template <class TTimeUnit>
std::string getTimeUnit()
{
    CIE_THROW(Exception, "Unknown time unit")
}

template <>
std::string getTimeUnit<std::chrono::milliseconds>()
{
    return "ms";
}

template <>
std::string getTimeUnit<std::chrono::microseconds>()
{
    return "us";
}

template <>

std::string getTimeUnit<std::chrono::nanoseconds>()
{
    return "ns";
}
} // unnamed namespace


template <class T>
Profiler<T>::Item::Item(std::source_location&& r_location)
    : _recursionLevel(0),
      _callCount(0),
      _time(0),
      _location(std::move(r_location))
{
}


template <class T>
typename Profiler<T>::Item& Profiler<T>::Item::operator+=(const Item& r_rhs)
{
    _callCount += r_rhs._callCount;
    _time += r_rhs._time;
    return *this;
}


template <class T>
Profiler<T>::Profiler()
    : Profiler([]() -> std::filesystem::path {
        // Get the time unit string
        utils::FileManager manager;
        auto path = manager.getCurrentDirectory() / ("cie_profiler_output_" + getTimeUnit<T>() + ".json");

        // Get a unique output file for this thread and rank
        manager.makeConcurrentPrivate(path);
        return path;
    }())
{
}


template <class T>
Profiler<T>::Profiler(std::filesystem::path&& r_outputPath)
    : _itemContainerMap(),
      _p_scope(),
      _outputPath(std::move(r_outputPath))
{
    // Reserve thread map
    const auto numberOfThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.reserve(numberOfThreads);
    std::atomic<unsigned> threadCounter = 0;
    for (unsigned i_thread=0; i_thread<numberOfThreads; ++i_thread) {
        threads.emplace_back([i_thread, &threadCounter, this](){
            while (threadCounter < i_thread) {} // <== wait until the previous thread finishes
            _itemContainerMap.emplace(std::this_thread::get_id(), std::list<Item> {});
            ++threadCounter;
        });
    }
    for (auto& r_thread : threads)
        r_thread.join();

    // Insert first item measuring the total runtime
    auto& r_item = this->create(std::source_location());
    _p_scope.reset(new Scope(r_item));
}


template <class T>
typename Profiler<T>::Item& Profiler<T>::create(std::source_location&& r_item)
{
    auto& r_list = _itemContainerMap[std::this_thread::get_id()];
    r_list.emplace_back(std::move(r_item));
    return r_list.back();
}


struct SourceLocationHash
{
    std::size_t operator()(const std::source_location& r_argument) const
    {
        std::string string(r_argument.file_name());
        string.append(std::to_string(r_argument.line()));
        return std::hash<std::string>()(string);
    }
};


struct SourceLocationEquals
{
    bool operator()(const std::source_location& r_lhs,
                    const std::source_location& r_rhs) const
    {
        return (std::string(r_lhs.file_name()) == std::string(r_rhs.file_name())) && (r_lhs.line() == r_rhs.line());
    }
};


template <class T>
Profiler<T>::~Profiler()
{
    // Stop the total timer
    _p_scope.reset();
    const auto totalTime = _itemContainerMap[std::this_thread::get_id()].front()._time;
    _itemContainerMap[std::this_thread::get_id()].pop_front();

    std::unordered_map<
        std::source_location,
        Item,
        SourceLocationHash,
        SourceLocationEquals
    > aggregateMap;

    // Combine maps from all threads
    for (const auto& r_threadMapPair : _itemContainerMap)
    {
        for (const auto& r_item : r_threadMapPair.second)
        {
            auto it = aggregateMap.find(r_item._location);
            if (it == aggregateMap.end())
                it = aggregateMap.emplace(r_item._location, r_item).first;
            else
                it->second += r_item;
        } // for item in vector
    } // for (location, vector) in map

    // Sort items based on their total duration
    using Numeric = typename Duration::rep;
    std::vector<std::tuple<std::source_location,Size,Numeric>> items;
    items.reserve(aggregateMap.size());
    for (const auto& r_pair : aggregateMap)
        items.emplace_back(r_pair.first, r_pair.second._callCount, r_pair.second._time.count());
    std::sort(items.begin(),
              items.end(),
              [](const auto& r_lhs, const auto& r_rhs)
                {return std::get<Numeric>(r_lhs) < std::get<Numeric>(r_rhs);});

    io::JSONObject root;
    std::vector<io::JSONObject> objects;

    // Add metadata
    {
        io::JSONObject object;
        object.add("name", "total");
        object.add("timeUnit", getTimeUnit<T>());
        object.add("time", totalTime.count());

        object.add("mpi", isMPIEnabled());

        {
            io::JSONObject build;
            build.add("compiler", getCompiler());
            build.add("compilerVersion", getCompilerVersion());
            build.add("compilerFlags", getCompileFlags());
            build.add("platform", getCompilePlatform());
            object.add("build", std::move(build));
        }

        root.add("meta", std::move(object));
    }

    // Add all items
    for (const auto& r_item : items)
    {
        objects.emplace_back();
        auto& r_object = objects.back();
        const auto& r_location = std::get<0>(r_item);
        r_object.add("file", std::string(r_location.file_name()));
        r_object.add("line", int(r_location.line()));
        r_object.add("function", std::string(std::get<0>(r_item).function_name()));
        r_object.add("callCount", std::get<1>(r_item));

        std::stringstream stream;
        stream << std::get<2>(r_item);
        r_object.add("time", stream.str());
    }
    root.add("results", objects);

    std::ofstream file(_outputPath);
    file << root << std::endl;
}


template <class T>
Profiler<T>& ProfilerSingleton<T>::get() noexcept
{
    std::scoped_lock<std::mutex> lock(_mutex);
    if (!_profiler.has_value())
        _profiler.emplace();

    return _profiler.value();
}


template <class T>
std::optional<Profiler<T>> ProfilerSingleton<T>::_profiler;


template <class T>
std::mutex ProfilerSingleton<T>::_mutex;


template class Profiler<std::chrono::milliseconds>;
template class ProfilerSingleton<std::chrono::milliseconds>;


template class Profiler<std::chrono::microseconds>;
template class ProfilerSingleton<std::chrono::microseconds>;


template class Profiler<std::chrono::nanoseconds>;
template class ProfilerSingleton<std::chrono::nanoseconds>;


} // namespace cie::utils
