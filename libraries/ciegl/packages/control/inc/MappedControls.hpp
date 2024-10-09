#ifndef CIE_GL_MAPPED_CONTROLS_HPP
#define CIE_GL_MAPPED_CONTROLS_HPP

// --- Internal Includes ---
#include "packages/control/inc/ControlKey.hpp"

// --- STL Includes ---
#include <unordered_map>
#include <deque>
#include <string>
#include <functional>
#include <istream>


namespace cie::gl {


class MappedControls
{
public:
    using control_function       = std::function<void()>;
    using map_type               = std::unordered_map<KeyEnum,control_function>;
    using configuration_map_type = std::unordered_map<ControlKey,control_function>;
    using configuration_contents = std::deque<std::pair<std::string,ControlKey>>;

public:
    MappedControls();

    void configure( configuration_map_type& r_configMap );

    /** @details Load control configuration from a json file with the following layout:
     *  @code
     *  "<control_name>" : ["key_caharacter_1", ..., "key_character_n"]
     *  .
     *  .
     *  .
     *  @endcode
     */
    void configure( std::istream& r_configurationFile );

    void callControl( KeyEnum key );

protected:
    virtual configuration_map_type makeConfigurationMap( configuration_contents& r_configContents ) = 0;

protected:
    map_type _controlMap;
};


} // namespace cie::gl

#include "packages/control/impl/MappedControls_impl.hpp"

#endif
