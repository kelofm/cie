#ifndef CIE_GL_SCENE_IMPL_HPP
#define CIE_GL_SCENE_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::gl {


template <class CameraType, class ...Args>
inline CameraPtr
Scene::makeCamera( Args&&... args )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto p_camera = CameraPtr(
        new CameraType( std::forward<Args>(args)... )
    );

    this->addCamera( p_camera );

    return p_camera;

    CIE_END_EXCEPTION_TRACING
}


template <class ...Args>
void
Scene::loadTexture( const std::string& r_name,
                    Args&&... args )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto& rp_texture = this->findTexture( r_name );

    rp_texture->write( std::forward<Args>(args)... );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl


#endif