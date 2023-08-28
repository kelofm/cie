// --- FEM Includes ---
#include "packages/utilities/inc/AssemblyMap.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie::fem {


using iterator = AssemblyMap::iterator;


using const_iterator = AssemblyMap::const_iterator;


iterator AssemblyMap::insert(Key&& r_localID, Value globalID)
{
    const auto pair = _map.emplace(std::move(r_localID), globalID);
    CIE_OUT_OF_RANGE_CHECK(
        pair.second,
        "Attempt to overwrite DoF with local ID (" << r_localID[0] << ", " << r_localID[1] << ")"
        << " and global ID (" << pair.first->second << ") with new global ID (" << globalID << ")";
    )

    return pair.first;
}


iterator AssemblyMap::insert(const Key& r_localID, Value globalID)
{
    return this->insert(Key(r_localID), globalID);
}


iterator AssemblyMap::insert(Index objectID, Index localIndex, Index globalID)
{
    return this->insert(Key {objectID, localIndex}, globalID);
}


Size AssemblyMap::size() const noexcept
{
    return _map.size();
}


const_iterator AssemblyMap::begin() const noexcept
{
    return _map.begin();
}


iterator AssemblyMap::begin() noexcept
{
    return _map.begin();
}


const_iterator AssemblyMap::end() const noexcept
{
    return _map.end();
}


iterator AssemblyMap::end() noexcept
{
    return _map.end();
}


AssemblyMap::AssemblyMap(AssemblyMap::Map&& r_map)
    : _map(std::move(r_map))
{
}


AssemblyMap::AssemblyMap(const AssemblyMap::Map& r_map)
    : _map(r_map)
{
}


} // namespace cie::fem
