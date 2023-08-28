// --- Internal Includes ---
#include "packages/control/inc/ControlKey.hpp"


namespace cie::gl {


ControlKey::ControlKey( ControlKey::key_type key ) :
    _keys()
{
    this->add( key );
}


ControlKey& ControlKey::add( ControlKey::key_type key )
{
    this->_keys.insert( key );
    return *this;
}


ControlKey& ControlKey::remove( ControlKey::key_type key )
{
    this->_keys.erase( key );
    return *this;
}


void ControlKey::clear()
{
    this->_keys.clear();
}


void ControlKey::set( ControlKey::key_type key )
{
    this->_keys.clear();
    this->add( key );
}


ControlKey& ControlKey::replace( ControlKey::key_type oldKey, ControlKey::key_type newKey )
{
    if ( this->_keys.erase(oldKey) != 0 )
        this->add( newKey );

    return *this;
}


ControlKey::operator std::string() const
{
    std::string aggregate;
    aggregate.reserve( this->_keys.size() );

    // TODO: what if there's a key with value 0?
    for ( auto key : this->_keys )
        aggregate.push_back( key );

    return aggregate;
}


const ControlKey::key_container& ControlKey::keys() const
{
    return this->_keys;
}


bool operator==( ControlKey::key_type lhs, const ControlKey& r_rhs )
{
    return r_rhs._keys.contains(lhs);
}


bool operator==( const ControlKey& r_lhs, ControlKey::key_type rhs )
{
    return rhs == r_lhs;
}


bool operator==( const ControlKey& r_lhs, const ControlKey& r_rhs )
{
    auto it_lhs = r_lhs._keys.begin();
    auto it_rhs = r_rhs._keys.begin();

    auto it_lhsEnd = r_lhs._keys.end();
    auto it_rhsEnd = r_rhs._keys.end();

    while ( it_lhs!=it_lhsEnd && it_rhs!=it_rhsEnd )
    {
        if ( (*it_lhs) < (*it_rhs) )
            ++it_lhs;
        else if ( (*it_lhs) > (*it_rhs) )
            ++it_rhs;
        else
            return true;
    }

    return false;
}


} // namespace cie::gl