#ifndef CIE_FEM_ASSEMBLY_MAP_IMPL_HPP
#define CIE_FEM_ASSEMBLY_MAP_IMPL_HPP

// --- FEM Includes ---
#include "packages/utilities/inc/AssemblyMap.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::fem {


inline AssemblyMap::Index AssemblyMap::at(const AssemblyMap::Key& r_key) const
{
    const auto it = _map.find(r_key);
    CIE_OUT_OF_RANGE_CHECK(
        it != _map.end(),
        "No DoF found for local ID: (" << r_key[0] << ", " << r_key[1] << ")"
    );
    return it->second;
}


inline AssemblyMap::Index AssemblyMap::at(AssemblyMap::Index objectID, AssemblyMap::Index localIndex) const
{
    return this->at({objectID, localIndex});
}


template <concepts::Iterator<AssemblyMap::Index> TInputIT, concepts::Iterator<AssemblyMap::Index> TOutputIT>
inline void AssemblyMap::at(AssemblyMap::Index objectID,
                            TInputIT it_localBegin,
                            const TInputIT it_localEnd,
                            TOutputIT it_globalBegin) const
{
    for (; it_localBegin!=it_localEnd; ++it_localBegin, ++it_globalBegin)
        *it_globalBegin++ = this->at(objectID, *it_localBegin);
}


template <concepts::Iterator<AssemblyMap::Index> TObjectIT>
inline void AssemblyMap::getObjectIDs(Index globalID, TObjectIT it_objectIDBegin) const
{
    for (const auto& r_pair : _map)
        if (r_pair.second == globalID)
            *it_objectIDBegin++ = *r_pair.first.begin();
}


template <concepts::Iterator<AssemblyMap::Index> TLocalIT, concepts::Iterator<AssemblyMap::Index> TGlobalIT>
inline AssemblyMap::iterator
AssemblyMap::insert(AssemblyMap::Index objectID,
                    TLocalIT it_localBegin,
                    const TLocalIT it_localEnd,
                    TGlobalIT it_globalBegin)
{
    CIE_OUT_OF_RANGE_CHECK(it_localBegin != it_localEnd)
    auto it = this->insert(objectID, it_localBegin++, it_localEnd++);
    for (; it_localBegin!=it_localEnd; ++it_localBegin, ++it_localEnd)
        this->insert({objectID, *it_localBegin}, *it_globalBegin);
    return it;
}


} // namespace cie::fem


#endif
