#ifndef CIE_UTILITIES_CACHE_HPP
#define CIE_UTILITIES_CACHE_HPP

// --- Ínternal Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- STL Includes ---
#include <unordered_map>
#include <functional>


namespace cie::utils {


template <  class InputType,
            class StoredType >
class AbsCache
{
public:
    using input_type            = InputType;
    using stored_type           = StoredType;
    using internal_type         = std::unordered_map<Size,StoredType>;
    using internal_iterator     = typename internal_type::const_iterator;
    using generator_function    = std::function<stored_type(const input_type&)>;

public:
    virtual Size hash( const input_type& input ) const = 0;
    internal_iterator insert(   const input_type& input,
                                generator_function generator,
                                bool force = false );
    internal_iterator insert(   const input_type& input,
                                const StoredType& value,
                                bool force = false );
    const stored_type& operator[]( Size inputID ) const;
    const stored_type& operator[]( const input_type& input ) const;

    bool cached( Size id ) const;
    bool cached( const InputType& input ) const;

    void erase( Size id );
    void clear()                                                { _map.clear(); }
    Size size() const                                           { return _map.size(); }
    internal_iterator begin() const                             { return _map.begin(); }
    internal_iterator end() const                               { return _map.end(); }

protected:
    internal_type   _map;
};


template <  concepts::Container InputContainer,
            concepts::Container OutputContainer>
class ContainerCache : public AbsCache<InputContainer,OutputContainer>
{
public:
    virtual Size hash( const InputContainer& input) const override;
};


} // namespace cie::utils

#include "packages/cache/impl/cache_impl.hpp"

#endif