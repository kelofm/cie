#ifndef CIE_UTILS_TAGS_HPP
#define CIE_UTILS_TAGS_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Inlcudes ---
#include <bitset>


namespace cie::tags {


namespace detail {


#define CIE_MAX_NUMBER_OF_TAGS 64

using Flags = std::bitset<CIE_MAX_NUMBER_OF_TAGS>;


template <Size ID, class TSelf>
class Tag
{
public:
    using Self = TSelf;

public:
    static constexpr Size id();

    /// @brief Shadow this in derived classes to specify compatibility.
    /// @note Compatible with all by default.
    static constexpr Flags getCompatibility();

    /// @brief Check compatiblity with another Tag.
    template <class ...TTags>
    static constexpr inline bool isCompatibleWith();
};


} // namespace detail


// --- Execution policy ---
struct Lazy;
struct Eager;

// --- Parallelism ---
struct Serial;
struct SMP;
struct MPI;

// --- Representation ---
struct Binary;
struct Text;


/// @addtogroup cieutils
/// @{


/// @brief Tag to indicate a void property.
struct Null : public detail::Tag<0,Null>
{};


/// @brief Tag to indicate lazy execution.
struct Lazy : public detail::Tag<1,Lazy>
{
    static constexpr detail::Flags getCompatibility();
};


/// @brief Tag to indicate eager execution.
struct Eager : public detail::Tag<2,Eager>
{
    static constexpr detail::Flags getCompatibility();
};


/// @brief Tag to indicate a lack of parallelism.
struct Serial : public detail::Tag<3,Serial>
{
    static constexpr detail::Flags getCompatibility();
};


/// @brief Tag to indicate shared memory parallelism.
struct SMP : public detail::Tag<4,SMP>
{
    static constexpr detail::Flags getCompatibility();
};


///@brief Tag to indicate process parallelism.
struct MPI : public detail::Tag<5,MPI>
{
    static constexpr detail::Flags getCompatibility();
};


/// @brief Tag to indicate binary representation.
struct Binary : public detail::Tag<6,Binary>
{
    static constexpr detail::Flags getCompatibility();
};


/// @brief Tag to indicate text representation.
struct Text : public detail::Tag<7,Text>
{
    static constexpr detail::Flags getCompatibility();
};


/// @}


#undef CIE_MAX_NUMBER_OF_TAGS


} // namespace cie::tags

#include "packages/types/impl/tags_impl.hpp"

#endif
