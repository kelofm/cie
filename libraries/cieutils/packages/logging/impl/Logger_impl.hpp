#ifndef CIE_UTILS_LOGGER_IMPL_HPP
#define CIE_UTILS_LOGGER_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <tuple>
#include <sstream>


namespace cie::utils {


template <class ...Args>
Logger& Logger::logs( Args&&... args )
{
    auto argTuple = std::make_tuple( std::forward<Args>(args)... );

    std::stringstream stream;
    std::apply(
        [&]( auto&&... arguments ) { ((stream << arguments),...); },
        argTuple
    );

    this->log( stream.str() );

    return *this;
}


template <class ExceptionType>
Logger& Logger::error( const std::string& r_message )
{
    log( "ERROR: " + r_message );
    CIE_THROW( ExceptionType, r_message );
    return *this;
}


} // namespace cie::utils


#endif