#ifndef CIE_UTILS_IO_JSON_IMPL_HPP
#define CIE_UTILS_IO_JSON_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::io {


template <concepts::io::SupportedType ValueType>
JSONObject& JSONObject::add( const std::string& r_key,
                             const ValueType& r_value,
                             bool allowOverwrite )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( this->hasKey(r_key) )
    {
        if ( !allowOverwrite )
            CIE_THROW( Exception, "'" + r_key + "' is an existing key" )

        this->operator[]( r_key ).set( r_value );
    }
    else
    {
        this->addDefault<ValueType>( r_key ).set( r_value );
    }

    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::io::SupportedType ValueType>
JSONObject& JSONObject::set( const ValueType& r_value )
{
    CIE_BEGIN_EXCEPTION_TRACING

    JSONObject::SetGet<ValueType>::set( *this, r_value );
    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <concepts::io::SupportedType ValueType>
ValueType JSONObject::as() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    return JSONObject::SetGet<ValueType>::as( *this );

    CIE_END_EXCEPTION_TRACING
}


template <concepts::io::SupportedType ValueType>
bool JSONObject::is() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    return JSONObject::TypeQuery<ValueType>::is( *this );

    CIE_END_EXCEPTION_TRACING
}


template <concepts::io::SupportedType ValueType>
JSONObject JSONObject::addDefault( const std::string& r_key )
{
    CIE_BEGIN_EXCEPTION_TRACING

    return JSONObject::TypeQuery<ValueType>::addDefault(
        *this,
        r_key
    );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::io


#endif