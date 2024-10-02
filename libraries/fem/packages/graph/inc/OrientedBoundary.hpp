#ifndef CIE_FEM_ORIENTED_BOUNDARY_HPP
#define CIE_FEM_ORIENTED_BOUNDARY_HPP

// --- Utility Includes ---
#include "packages/maths/inc/bit.hpp" // getMSBIndex
#include "packages/stl_extension/inc/Hash.hpp" // Hash

// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp" // BoundaryID
#include "packages/graph/inc/OrientedAxes.hpp" // OrientedAxes

// --- STL Includes ---
#include <bitset> // std::bitset, std::hash
#include <iterator> // std::random_access_iterator_tag
#include <cstdint> // std::int8_t
#include <iosfwd> // std::ostream


namespace cie::fem {


/** @brief Class defining a \f$(d-1)\f$-dimensional cube spanning \f$[-1, 1]^{d-1}\f$
 *         in any axis-parallel rotation within \f$d\f$-dimensional space.
 *  @details This class essentially defines a \f$d\f$-dimensional cartesian coordinate system,
 *           and an axis-aligned unit vector in it. The vector is the normal of the face the
 *           instance refers to, while the coordinate system defines which direction the positive
 *           axes are pointing in that plane. This is necessary for comparing ansatz functions on
 *           coincident faces of adjacent elements.
 *
 *           For example, the following case demonstrates the rotated local coordinate systems of
 *           6 elements in 2D in topological space.
 *           \code
 *           + ----------- +                           + ----------- +                           + ----------- +
 *           |             |                           |             |                           |             |
 *           |             |                           |             |                           |             |
 *           |      + > x  | "+x-y:+x"       "-x-y:+x" |  x < +      | "-x-y:-x"       "+x-y:-x" |      + > x  |
 *           |      v      |                           |      v      |                           |      v      |
 *           |      y      |                           |      y      |                           |      y      |
 *           + ----------- +                           + ----------- +                           + ----------- +
 *              "+x-y:+y"                                 "-x-y:+y"                                 "-x-y:+y"
 *
 *              "+x+y:+y"                                 "-x+y:+y"                                 "-x+y:+y"
 *           + ----------- +                           + ----------- +                           + ----------- +
 *           |      y      |                           |      y      |                           |      y      |
 *           |      ^      |                           |      ^      |                           |      ^      |
 *           |      + > x  | "+x+y:+x"       "-x+y:+x" |  x < +      | "-x+y:-x"       "+x+y:-x" |      + > x  |
 *           |             |                           |             |                           |             |
 *           |             |                           |             |                           |             |
 *           + ----------- +                           + ----------- +                           + ----------- +
 *           \endcode
 */
template <unsigned Dimension>
class OrientedBoundary
{
private:
    static constexpr unsigned ComponentBitWidth = 1 + 1 + utils::getMSBIndex(Dimension);

    using Data = std::bitset<(Dimension + 1) * ComponentBitWidth>;

    template <bool Mutable>
    class ValueProxy;

    template <bool Mutable>
    class Iterator;

public:
    static constexpr unsigned dimension = Dimension;

    using value_type = BoundaryID;

    using reference = ValueProxy</*Mutable=*/true>;

    using const_reference = ValueProxy</*Mutable=*/false>;

    using pointer = ValueProxy</*Mutable=*/true>;

    using size_type = std::uint8_t;

    using iterator = Iterator</*Mutable=*/true>;

    using const_iterator = Iterator</*Mutable=*/false>;

public:
    OrientedBoundary() noexcept;

    OrientedBoundary(Ptr<const BoundaryID> itBegin,
                     size_type size,
                     BoundaryID id);

    OrientedBoundary(OrientedAxes<Dimension> axes,
                     BoundaryID id);

    OrientedBoundary(const char axes[2 * Dimension + 1],
                     const char id[3])
    requires (Dimension < 4);

    /// @brief @ref BoundaryID of the current instance in topological space.
    [[nodiscard]] BoundaryID id() const noexcept;

    /// @copydoc OrientedBoundary::id() const
    [[nodiscard]] reference id() noexcept;

    /// @brief @ref BoundaryID of the face in the unrotated system.
    [[nodiscard]] BoundaryID localID() const noexcept;

    [[nodiscard]] bool operator==(OrientedBoundary rhs) const noexcept;

    [[nodiscard]] bool operator!=(OrientedBoundary rhs) const noexcept;

    /// @details Checks the ordering of the ID, and continues with a
    ///          lexicographical comparison of the rotated axes in case of equality.
    [[nodiscard]] bool operator<(OrientedBoundary rhs) const noexcept;

    /// @brief Increment the boundary within the rotated system.
    OrientedBoundary& operator++() noexcept;

    /// @copydoc OrientedBoundary::operator++()
    OrientedBoundary operator++(int) noexcept;

    /// @brief Flip the boundary to represent the opposite face in the same system.
    [[nodiscard]] OrientedBoundary operator-() const noexcept;

    /// @brief Return the size of the array of @ref BoundaryID "boundary IDs" that define the topological space.
    [[nodiscard]] constexpr size_type size() const noexcept
    {return dimension;}

    [[nodiscard]] const_reference operator[](size_type index) const;

    [[nodiscard]] reference operator[](size_type index);

    [[nodiscard]] const_reference at(size_type index) const;

    [[nodiscard]] reference at(size_type index);

    [[nodiscard]] const_iterator begin() const noexcept;

    [[nodiscard]] iterator begin() noexcept;

    [[nodiscard]] const_iterator end() const noexcept;

    [[nodiscard]] iterator end() noexcept;

private:
    /** @brief Data storing the orientation of the local system as well as the @ref BoundaryID in that space.
     *  @details The stored data is an array of bits that consists of <tt>dimension + 1</tt> bit sets of size
     *           @p ComponentBitWidth. The first @p dimension sets define the new orientations of each axis
     *           at the corresponding index, uniquely defining the rotated system. The last set defines which
     *           boundary the current instance refers to in topological space.
     */
    Data _data;

    friend struct utils::Hash<OrientedBoundary>;

    template <bool Mutable>
    class ValueProxy
    {
    public:
        ValueProxy& operator=(BoundaryID rhs)
        requires Mutable;

        operator BoundaryID() const noexcept;

    private:
        using Data = typename OrientedBoundary<Dimension>::Data;

        ValueProxy() = delete;

        ValueProxy(Ref<std::conditional_t<Mutable,Data,const Data>> rData,
                   unsigned index);

        ValueProxy& operator=(ValueProxy rhs) = delete;

        friend class OrientedBoundary<Dimension>;
        Ptr<std::conditional_t<Mutable,Data,const Data>> _pData;

        template <bool M>
        friend class OrientedBoundary<Dimension>::Iterator;

        std::int8_t _index;
    }; // class ValueProxy


    template <bool Mutable>
    class Iterator
    {
    private:
        using Data = typename OrientedBoundary<Dimension>::Data;

    public:
        using value_type = ValueProxy<Mutable>;

        using reference = ValueProxy<Mutable>;

        using pointer = ValueProxy<Mutable>;

        using difference_type = std::int8_t;

        using iterator_category = std::random_access_iterator_tag;

    public:
        Iterator(Iterator&& rhs) noexcept = default;

        Iterator(const Iterator& rhs) noexcept = default;

        Iterator(Iterator</*Mutable=*/true>&& rhs) noexcept
        requires (!Mutable);

        Iterator(const Iterator</*Mutable=*/true>& rhs) noexcept
        requires (!Mutable);

        reference operator*() const;

        Iterator& operator++() noexcept;

        Iterator operator++(int) noexcept;

        Iterator& operator+=(difference_type offset) noexcept;

        Iterator operator+(difference_type offset) noexcept;

        Iterator& operator--() noexcept;

        Iterator operator--(int) noexcept;

        Iterator& operator-=(difference_type offset) noexcept;

        Iterator operator-(difference_type offset) noexcept;

        difference_type operator-(Iterator rhs) noexcept;

        template <bool M>
        bool operator==(Iterator<M> rhs) noexcept;

        template <bool M>
        bool operator!=(Iterator<M> rhs) noexcept;

    private:
        friend class Iterator<!Mutable>;

        friend class OrientedBoundary;

        Iterator() = delete;

        Iterator(Ref<std::conditional_t<Mutable,Data,const Data>> rData,
                 unsigned index) noexcept;

        Ptr<std::conditional_t<Mutable,Data,const Data>> _pData;

        std::int8_t _index;
    }; // class Iterator
}; // class OrientedBoundary


template <unsigned Dimension>
std::ostream& operator<<(std::ostream& rStream, OrientedBoundary<Dimension> boundary);


} // namespace cie::fem


namespace cie::utils {


template <unsigned Dimension>
struct Hash<fem::OrientedBoundary<Dimension>>
{
    auto operator()(fem::OrientedBoundary<Dimension> instance) const noexcept
    {
        return std::hash<typename fem::OrientedBoundary<Dimension>::Data>()(instance._data);
    }
}; // struct hash<OrientedBoundary>



} // namespace cie::utils


#include "packages/graph/impl/OrientedBoundary_impl.hpp"

#endif
