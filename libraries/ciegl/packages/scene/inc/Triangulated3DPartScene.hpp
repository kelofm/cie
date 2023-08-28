#ifndef CIE_GL_TRIANGULATED_3D_PART_SCENE_HPP
#define CIE_GL_TRIANGULATED_3D_PART_SCENE_HPP

// --- Internal Includes ---
#include "packages/scene/inc/PartScene.hpp"


namespace cie::gl {


class Triangulated3DPartScene : public PartScene
{
public:
    // Constructor for the AbsWindow::makeScene
    Triangulated3DPartScene( const std::string& r_name,
                             Triangulated3DPartScene::part_container&& r_parts,
                             utils::Logger& r_logger = utils::LoggerSingleton::get() );

protected:
    virtual bool update_impl() override;

    /// Constructor for derived classes
    Triangulated3DPartScene( const std::string& r_name,
                             part_container&& r_parts,
                             CameraPtr p_camera,
                             Shader::SharedPointer p_vertexShader,
                             Shader::SharedPointer p_geometryShader,
                             Shader::SharedPointer p_fragmentShader,
                             utils::Logger& r_logger,
                             VertexBuffer::SharedPointer p_vertexBuffer,
                             ElementBuffer::SharedPointer p_elementBuffer );
};


} // namespace cie::gl


#endif