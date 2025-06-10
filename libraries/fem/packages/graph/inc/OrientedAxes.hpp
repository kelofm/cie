#ifndef CIE_FEM_ORIENTED_AXES_HPP
#define CIE_FEM_ORIENTED_AXES_HPP

// --- Utility Includes ---
#include "packages/maths/inc/bit.hpp" // getMSBIndex
#include "packages/stl_extension/inc/Hash.hpp" // Hash

// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp" // BoundaryID
#include "packages/io/inc/GraphML.hpp" // GraphML::Serializer

// --- STL Includes ---
#include <bitset> // std::bitset, std::hash
#include <iterator> // std::random_access_iterator_tag
#include <cstdint> // std::int8_t
#include <iosfwd> // std::ostream


namespace cie::fem {


/** @brief Class defining a \f$d\f$-dimensional cube spanning \f$[-1, 1]^d\f$
 *         in any axis-parallel orientation within \f$d\f$-dimensional space.
 */
template <unsigned Dimension>
class OrientedAxes
{
private:
    static constexpr unsigned ComponentBitWidth = 1 + 1 + utils::getMSBIndex(Dimension);

    using Data = std::bitset<Dimension * ComponentBitWidth>;

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
    OrientedAxes() noexcept;

    OrientedAxes(Ptr<const BoundaryID> itBegin,
                 size_type size);

    OrientedAxes(const char axes[2 * Dimension + 1])
    requires (Dimension < 4);

    [[nodiscard]] bool operator==(OrientedAxes rhs) const noexcept;

    [[nodiscard]] bool operator!=(OrientedAxes rhs) const noexcept;

    /// @details Lexicographical comparison of the rotated axes.
    [[nodiscard]] bool operator<(OrientedAxes rhs) const noexcept;

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
    /** @brief Data storing the orientation of the local system.
     *  @details The stored data is an array of bits that consists of <tt>dimension</tt> bit sets of size
     *           @p ComponentBitWidth. The sets define the new orientations of each axis
     *           at the corresponding index, uniquely defining the rotated system.
     */
    Data _data;

    friend struct utils::Hash<OrientedAxes>;

    template <bool Mutable>
    class ValueProxy
    {
    public:
        ValueProxy& operator=(BoundaryID rhs)
        requires Mutable;

        operator BoundaryID() const noexcept;

    private:
        using Data = typename OrientedAxes<Dimension>::Data;

        ValueProxy() = delete;

        ValueProxy(Ref<std::conditional_t<Mutable,Data,const Data>> rData,
                   unsigned index);

        ValueProxy& operator=(ValueProxy rhs) = delete;

        friend class OrientedAxes<Dimension>;
        Ptr<std::conditional_t<Mutable,Data,const Data>> _pData;

        template <bool M>
        friend class OrientedAxes<Dimension>::Iterator;

        std::int8_t _index;
    }; // class ValueProxy


    template <bool Mutable>
    class Iterator
    {
    private:
        using Data = typename OrientedAxes<Dimension>::Data;

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

        friend class OrientedAxes;

        Iterator() = delete;

        Iterator(Ref<std::conditional_t<Mutable,Data,const Data>> rData,
                 unsigned index) noexcept;

        Ptr<std::conditional_t<Mutable,Data,const Data>> _pData;

        std::int8_t _index;
    }; // class Iterator
}; // class OrientedAxes


template <unsigned D>
std::ostream& operator<<(std::ostream& rStream, OrientedAxes<D> boundary);


template <unsigned Dimension>
struct io::GraphML::Serializer<OrientedAxes<Dimension>>
{
    void header(Ref<XMLElement> rElement) noexcept;

    void operator()(Ref<XMLElement> rElement, Ref<const OrientedAxes<Dimension>> rObject) noexcept;
}; // GraphML::Serializer<OrientedAxes<Dimension>>



} // namespace cie::fem


namespace cie::utils {


template <unsigned Dimension>
struct Hash<fem::OrientedAxes<Dimension>>
{
    auto operator()(fem::OrientedAxes<Dimension> instance) const noexcept
    {
        return std::hash<typename fem::OrientedAxes<Dimension>::Data>()(instance._data);
    }
}; // struct hash<OrientedAxes>


} // namespace cie::utils


#include "packages/graph/impl/OrientedAxes_impl.hpp"

#endif
