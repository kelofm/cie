#ifndef CIE_UTILS_PROFILER_HPP
#define CIE_UTILS_PROFILER_HPP

// --- STL Includes ---
#include <unordered_map>
#include <filesystem>
#include <source_location>
#include <chrono>
#include <thread>
#include <optional>
#include <list>
#include <mutex>


namespace cie::utils {


/// @ingroup cieutils
template <class TTimeUnit>
class Profiler
{
private:
    using TimeUnit = TTimeUnit;

    using Duration = TimeUnit;

    using Clock = std::chrono::high_resolution_clock;

    class Item
    {
    public:
        Item(std::source_location&& r_location);

    private:
        Item& operator+=(const Item& r_rhs);

    private:
        friend class Profiler;

        unsigned _recursionLevel;

        std::size_t _callCount;

        Duration _time;

        std::source_location _location;
    }; // class Item

public:
    class Scope
    {
    public:
        ~Scope();

    private:
        void start();

        void stop();

        Scope(Item& r_item);

        Scope(Scope&&) = delete;

        Scope(const Scope&) = delete;

        Scope& operator=(Scope&&) = delete;

        Scope& operator=(const Scope&) = delete;

    private:
        friend class Profiler;

        Item& _r_item;

        std::chrono::high_resolution_clock::time_point _begin;
    }; // class Scope

public:
    Profiler();

    Profiler(Profiler&& r_rhs) = default;

    Profiler(std::filesystem::path&& r_outputPath);

    ~Profiler();

    Profiler& operator=(Profiler&& r_rhs) = default;

    [[nodiscard]] Item& create(std::source_location&& r_item);

    [[nodiscard]] Scope profile(Item& r_item);

private:
    Profiler(const Profiler&) = delete;

    Profiler& operator=(const Profiler&) = delete;

private:
    std::unordered_map<std::thread::id,std::list<Item>> _itemContainerMap;

    std::unique_ptr<Scope> _p_scope;

    std::filesystem::path _outputPath;
}; // class Profiler


template <class TTimeUnit>
class ProfilerSingleton
{
public:
    static Profiler<TTimeUnit>& get() noexcept;

private:
    static std::optional<Profiler<TTimeUnit>> _profiler;

    static std::mutex _mutex;
}; // class ProfilerSingleton


} // namespace cie::utils

#include "packages/logging/impl/Profiler_impl.hpp"

#endif
