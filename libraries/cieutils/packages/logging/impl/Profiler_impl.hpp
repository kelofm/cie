#ifndef CIE_UTILS_PROFILER_IMPL_HPP
#define CIE_UTILS_PROFILER_IMPL_HPP

// --- Utility Includes ---
#include "packages/logging/inc/Profiler.hpp"

// --- STL Includes ---
#include <string>
#include <chrono>


namespace cie::utils {


template <class T>
inline Profiler<T>::Scope::Scope(Profiler::Item& r_item)
    : _r_item(r_item),
      _begin(Clock::now())
{
    this->start();
}


template <class T>
inline Profiler<T>::Scope::~Scope()
{
    this->stop();
}


template <class T>
inline void Profiler<T>::Scope::start()
{
    _begin = Clock::now();
    ++_r_item._callCount;
    ++_r_item._recursionLevel;
}


template <class T>
inline void Profiler<T>::Scope::stop()
{
    if (!--_r_item._recursionLevel)
        _r_item._time += std::chrono::duration_cast<Profiler::TimeUnit>(Clock::now() - _begin);
}


template <class T>
inline typename Profiler<T>::Scope Profiler<T>::profile(Item& r_item)
{
    return Scope(r_item);
}


} // namespace cie::utils


#endif
