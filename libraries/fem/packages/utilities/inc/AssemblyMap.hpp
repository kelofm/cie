#ifndef CIE_FEM_ASSEMBLY_MAP_HPP
#define CIE_FEM_ASSEMBLY_MAP_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <map>


namespace cie::fem {


///@addtogroup fem
///@{

/** @brief Class managing the mapping between local and global DoFs.
 *
 *  @details Depending on whether a degree of freedom (DoF) is represented in
 *           local or global space, the data used to identify it changes. DoFs in
 *           global space are identified by a single @ref Size (@c globalIndex) while
 *           two are used to uniquely identify them in an object's local space. The
 *           first integer in a DoF's local identifier refers to the ID of the object
 *           it belongs to (an @p Element for example) while the second one identifies
 *           it within that object, so @a AssemblyMap is a \f( N^2 \rightarrow N \f)
 *           mapping of \f( \{objectID,localID\} : globalID \f).
 *  @todo Make accessing/inserting multiple DoFs within the same object more efficient.
 */
class AssemblyMap
{
public:
    using Index = Size;

    using Key = StaticArray<Index,2>;

private:
    using Value = Index;

    using Map = std::map<Key,Value>;

public:
    using iterator = Map::iterator;

    using const_iterator = Map::const_iterator;

    using value_type = Map::value_type;

public:
    AssemblyMap() = default;

    AssemblyMap(AssemblyMap&& rRhs) = default;

    explicit AssemblyMap(const AssemblyMap& rRhs) = default;

    AssemblyMap& operator=(AssemblyMap&& rRhs) = default;

    AssemblyMap& operator=(const AssemblyMap& rRhs) = delete;

    /** @brief Convert a DoF's identifier in local space to the global one.
     *
     *  @param rLocalID: identifier of the DoF ibn local space consisting of an object ID and a local index.
     */
    Index at(const Key& rLocalID) const;

    /** @brief Convert a DoF's identifier in local space to the global one.
     *
     *  @param objectID: ID of the object the DoF belongs to.
     *  @param localIndex: DoF's local index within the object.
     */
    Index at(Index objectID, Index localIndex) const;

    /** @brief Convert several DoFs' local identifiers of the same object at once.
     *
     *  @param objectID: ID of the object that owns all provided DoFs.
     *  @param itInputBegin: iterator pointing to the first DoF's local index.
     *  @param itInputEnd iterator pointing past the last DoF's local index.
     *  @param itOutputBegin output iterator specifying where the DoFs' global indices will be written to.
     */
    template <concepts::Iterator<Index> TInputIT, concepts::Iterator<Index> TOutputIT>
    void at(Index objectID,
            TInputIT itInputBegin,
            const TInputIT itInputEnd,
            TOutputIT itOutputBegin) const;

    /** @brief Get the indices of all objects that contain the specified DoF (identified by its ID in global space).
     *
     *  @param globalID: the DoF's identifier in global space.
     *  @param itObjectIDBegin: output iterator pointing to a memory location where object IDs will be written.
     */
    template <concepts::Iterator<Index> TObjectIT>
    void getObjectIDs(Index globalID, TObjectIT itObjectIDBegin) const;

    /** @brief Insert a DoF's local-global identifier pair.
     *
     *  @param rLocalID: Identifier of the DoF in local space consisting of an object ID and a local index.
     *  @param globalID: Identifier of the DoF in global space.
     */
    iterator insert(Key&& rLocalID, Value globalID);

    /** @brief Insert a DoF's local-global identifier pair.
     *
     *  @param rLocalID: Identifier of the DoF in local space consisting of an object ID and a local index.
     *  @param globalID: Identifier of the DoF in global space.
     */
    iterator insert(const Key& rLocalID, Value globalID);

    /** @brief Insert a DoF's local-global identifier pair.
     *
     *  @param objectID: ID of the object the DoF belongs to.
     *  @param localIndex: Index of the DoF withing its containing object.
     *  @param globalID: Identifier of the DoF in global space.
     */
    iterator insert(Index objectID, Index localIndex, Index globalID);

    /** @brief Insert several DoF's local-global identifier pairs sharing an object.
     *
     *  @param objectID: ID of the object all provided DoFs belong to.
     *  @param itLocalBegin: iterator pointing to the first DoF's index in local space.
     *  @param itLocalEnd: iterator past the last DoF's index in local space.
     *  @param itGlobalBegin: iterator pointing to the first DoF's identifier in global space.
     */
    template <concepts::Iterator<Index> TLocalIT, concepts::Iterator<Index> TGlobalIndex>
    iterator insert(Index objectID,
                    TLocalIT itLocalBegin,
                    const TLocalIT itLocalEnd,
                    TGlobalIndex itGlobalBegin);

    /// @brief Number of recorded DoFs.
    Size size() const noexcept;

    const_iterator begin() const noexcept;

    iterator begin() noexcept;

    const_iterator end() const noexcept;

    iterator end() noexcept;

private:
    AssemblyMap(Map&& rMap);

    explicit AssemblyMap(const Map& rMap);

private:
    Map _map;
}; // class AssemblyMap

///@}

} // namespace cie::fem

#include "packages/utilities/impl/AssemblyMap_impl.hpp"

#endif