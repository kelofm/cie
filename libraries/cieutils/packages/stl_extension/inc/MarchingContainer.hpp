#ifndef CIE_UTILS_MARCHING_CONTAINER_HPP
#define CIE_UTILS_MARCHING_CONTAINER_HPP

// --- Internal Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- STL Includes ---
#include <vector>


namespace cie::utils {


/// @brief A container that only holds a specified number of latest values.
/// @addtogroup cieutils
template <class ValueType>
class MarchingContainer
{
public:
    using value_type = ValueType;
    using size_type  = Size;

public:
    MarchingContainer() = delete;
    MarchingContainer( const MarchingContainer<ValueType>& r_rhs ) = default;
    MarchingContainer<ValueType>& operator=( const MarchingContainer<ValueType>& r_rhs ) = default;

    /**
     * Construct a marching container of size identical to the size
     * of the argument container, and initialize it with its values.
     */
    template <class ContainerType>
    requires concepts::Container<ContainerType,ValueType>
    MarchingContainer( const ContainerType& r_initializer );

    /**
     * Construct a marching container of the specified size, initialized
     * with the specified initializer
     */
    MarchingContainer( Size capacity, const ValueType& r_initializer );

    /**
     * Construct a marching container of the specified size,
     * with default initialized values
     */
    MarchingContainer( Size capacity );

    template <class ...Args>
    void emplace_back( Args&&... args );

    void push_back( const ValueType& r_value );

    ValueType& back();
    const ValueType& back() const;

    ValueType& operator[]( Size index );
    const ValueType& operator[]( Size index ) const;

    /** Return the number of times this container was pushed to
     * @note the number of values held will always be <capacity>
     */
    Size size() const;

private:
    bool outOfRangeCheck( Size index );

private:
    using container_type = std::vector<value_type>;

private:
    Size           _capacity;
    Size           _counter;
    container_type _container;
};


} // namespace cie::utils

#include "packages/stl_extension/impl/MarchingContainer_impl.hpp"

#endif