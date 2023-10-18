#ifndef CIE_MP_THREAD_UTILS_HPP
#define CIE_MP_THREAD_UTILS_HPP

// --- STL Includes ---
#include <unordered_set>
#include <thread>


namespace cie::mp {


/// @addtogroup cieutils
const std::unordered_set<std::thread::id>& getThreadIDs();


} // namespace cie::mp


#endif
