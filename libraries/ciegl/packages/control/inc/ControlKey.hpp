#ifndef CIE_GL_CONTROL_KEY_HPP
#define CIE_GL_CONTROL_KEY_HPP

// --- Internal Includes ---
#include "packages/control/inc/callback_types.hpp"

// --- STL Includes ---
#include <set>
#include <string>
#include <functional>


namespace cie::gl {


class ControlKey
{
public:
    using key_type      = KeyEnum;
    using key_container = std::set<key_type>;

public:
    ControlKey( key_type key );
    explicit ControlKey() = default;
    ControlKey( const ControlKey& r_rhs ) = default;
    ControlKey& operator=( const ControlKey& r_rhs ) = default;

    ControlKey& add( key_type key );
    ControlKey& remove( key_type key );

    void clear();
    void set( key_type key );
    ControlKey& replace( key_type oldKey, key_type newKey );

    explicit operator std::string() const;

    const key_container& keys() const;

    friend bool operator==( key_type lhs, const ControlKey& r_rhs );
    friend bool operator==( const ControlKey& r_lhs, key_type rhs );
    friend bool operator==( const ControlKey& r_lhs, const ControlKey& r_rhs );

protected:
    key_container _keys;
};


} // namespace cie::gl


// Hash overload for std::unordered_map and std::set
namespace std {


template <>
struct hash<cie::gl::ControlKey>
{
    size_t operator()( const cie::gl::ControlKey& r_key ) const;
};


} // namespace std

#include "packages/control/impl/ControlKey_impl.hpp"

#endif