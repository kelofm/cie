// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/scene/inc/Triangulated3DPartScene.hpp"
#include "packages/camera/inc/Camera.hpp"
#include "packages/camera/inc/PerspectiveProjection.hpp"
#include "packages/shaders/inc/Shader.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <filesystem>


namespace cie::gl {


const std::filesystem::path SHADER_DIR = getDataPath() / "ciegl/shaders/Triangulated3DPartScene";


Triangulated3DPartScene::Triangulated3DPartScene( const std::string& r_name,
                                                  Triangulated3DPartScene::part_container&& r_parts,
                                                  CameraPtr p_camera,
                                                  Shader::SharedPointer p_vertexShader,
                                                  Shader::SharedPointer p_geometryShader,
                                                  Shader::SharedPointer p_fragmentShader,
                                                  utils::Logger& r_logger,
                                                  VertexBuffer::SharedPointer p_vertexBuffer,
                                                  ElementBuffer::SharedPointer p_elementBuffer ) :
    PartScene( r_name,
               std::forward<Triangulated3DPartScene::part_container>(r_parts),
               p_vertexShader,
               p_geometryShader,
               p_fragmentShader,
               r_logger,
               p_vertexBuffer,
               p_elementBuffer )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Register camera or construct one if it was not provided
    if ( p_camera )
        this->addCamera( p_camera );
    else
        p_camera = this->makeCamera<Camera<PerspectiveProjection>>();

    // Bind required uniforms
    this->bindUniform( "transformation", p_camera->transformationMatrix() );
    this->bindUniform( "cameraPosition", p_camera->position() );

    CIE_END_EXCEPTION_TRACING
}


Triangulated3DPartScene::Triangulated3DPartScene( const std::string& r_name,
                                                  Triangulated3DPartScene::part_container&& r_parts,
                                                  utils::Logger& r_logger ) :
    Triangulated3DPartScene( r_name,
                             std::forward<Triangulated3DPartScene::part_container>(r_parts),
                             nullptr,
                             makeVertexShader(SHADER_DIR / "vertexShader.json",
                                              SHADER_DIR / "vertexShader.glsl" ),
                             makeGeometryShader(SHADER_DIR / "geometryShader.json",
                                                SHADER_DIR / "geometryShader.glsl" ),
                             makeFragmentShader(SHADER_DIR / "fragmentShader.json",
                                                SHADER_DIR / "fragmentShader.glsl" ),
                             r_logger,
                             nullptr,
                             nullptr )
{
}


bool Triangulated3DPartScene::update_impl()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( this->_partsRequireUpdate )
        this->updateParts();

    GLint64 numberOfIndices;
    glGetBufferParameteri64v( GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &numberOfIndices );
    numberOfIndices /= sizeof( GLuint );
    glDrawElements( GL_TRIANGLES, numberOfIndices, GL_UNSIGNED_INT, 0 );

    return true;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl