#ifndef CIE_GL_CONTROL_KEY_IMPL_HPP
#define CIE_GL_CONTROL_KEY_IMPL_HPP

namespace std {


inline size_t
hash<cie::gl::ControlKey>::operator()( const cie::gl::ControlKey& r_key ) const
{
    return hash<string>()( string(r_key) );
}


} // namespace std

#endif