#ifndef CIE_GL_MAPPED_CONTROLS_IMPL_HPP
#define CIE_GL_MAPPED_CONTROLS_IMPL_HPP

namespace cie::gl {


inline void
MappedControls::callControl( KeyEnum key )
{
    auto it_map = this->_controlMap.find( key );
    
    if ( it_map != this->_controlMap.end() )
        if ( it_map->second )
            it_map->second();
}


} // namespace cie::gl

#endif