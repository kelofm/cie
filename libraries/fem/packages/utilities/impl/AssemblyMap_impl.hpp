#ifndef CIE_FEM_ASSEMBLY_MAP_IMPL_HPP
#define CIE_FEM_ASSEMBLY_MAP_IMPL_HPP

// --- FEM Includes ---
#include "packages/utilities/inc/AssemblyMap.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::fem {


inline AssemblyMap::Index AssemblyMap::at(const AssemblyMap::Key& rKey) const
{
    const auto it = _map.find(rKey);
    CIE_OUT_OF_RANGE_CHECK(
        it != _map.end(),
        "No DoF found for local ID: (" << rKey[0] << ", " << rKey[1] << ")"
    );
    return it->second;
}


inline AssemblyMap::Index AssemblyMap::at(AssemblyMap::Index objectID, AssemblyMap::Index localIndex) const
{
    return this->at({objectID, localIndex});
}


template <concepts::Iterator<AssemblyMap::Index> TInputIT, concepts::Iterator<AssemblyMap::Index> TOutputIT>
inline void AssemblyMap::at(AssemblyMap::Index objectID,
                            TInputIT itLocalBegin,
                            const TInputIT itLocalEnd,
                            TOutputIT itGlobalBegin) const
{
    for (; itLocalBegin!=itLocalEnd; ++itLocalBegin, ++itGlobalBegin)
        *itGlobalBegin++ = this->at(objectID, *itLocalBegin);
}


template <concepts::Iterator<AssemblyMap::Index> TObjectIT>
inline void AssemblyMap::getObjectIDs(Index globalID, TObjectIT itObjectIDBegin) const
{
    for (const auto& rPair : _map)
        if (rPair.second == globalID)
            *itObjectIDBegin++ = *rPair.first.begin();
}


template <concepts::Iterator<AssemblyMap::Index> TLocalIT, concepts::Iterator<AssemblyMap::Index> TGlobalIT>
inline AssemblyMap::iterator
AssemblyMap::insert(AssemblyMap::Index objectID,
                    TLocalIT itLocalBegin,
                    const TLocalIT itLocalEnd,
                    TGlobalIT itGlobalBegin)
{
    CIE_OUT_OF_RANGE_CHECK(itLocalBegin != itLocalEnd)
    auto it = this->insert(objectID, itLocalBegin++, itLocalEnd++);
    for (; itLocalBegin!=itLocalEnd; ++itLocalBegin, ++itLocalEnd)
        this->insert({objectID, *itLocalBegin}, *itGlobalBegin);
    return it;
}


} // namespace cie::fem


#endif
