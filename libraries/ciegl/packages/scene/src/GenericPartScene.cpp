// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Incldues ---
#include "packages/scene/inc/GenericPartScene.hpp"
#include "packages/shaders/inc/Shader.hpp"


namespace cie::gl {


GenericPartScene::GenericPartScene( const std::string& r_name,
                                    Size drawMode,
                                    GenericPartScene::part_container&& r_parts,
                                    Shader::SharedPointer p_vertexShader,
                                    Shader::SharedPointer p_geometryShader,
                                    Shader::SharedPointer p_fragmentShader,
                                    utils::Logger& r_logger,
                                    VertexBuffer::SharedPointer p_vertexBuffer,
                                    ElementBuffer::SharedPointer p_elementBuffer ) :
    PartScene( r_name,
               std::forward<GenericPartScene::part_container>(r_parts),
               p_vertexShader,
               p_geometryShader,
               p_fragmentShader,
               r_logger,
               p_vertexBuffer,
               p_elementBuffer ),
    _drawMode( drawMode )
{
}


GenericPartScene::GenericPartScene( const std::string& r_name,
                                    Size drawMode,
                                    const std::filesystem::path& r_shaderDirectory,
                                    utils::Logger& r_logger ) :
    GenericPartScene( r_name,
                      drawMode,
                      {},
                      makeVertexShader(r_shaderDirectory / "vertexShader.json",
                                       r_shaderDirectory / "vertexShader.glsl" ),
                      makeGeometryShader(r_shaderDirectory / "geometryShader.json",
                                         r_shaderDirectory / "geometryShader.glsl" ),
                      makeFragmentShader(r_shaderDirectory / "fragmentShader.json",
                                         r_shaderDirectory / "fragmentShader.glsl" ),
                      r_logger,
                      nullptr,
                      nullptr )
{
}


void GenericPartScene::setDrawMode( Size drawMode )
{
    this->_drawMode = drawMode;
}


bool GenericPartScene::update_impl()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( this->_partsRequireUpdate )
        this->updateParts();

    GLint64 numberOfIndices;
    glGetBufferParameteri64v( GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &numberOfIndices );
    numberOfIndices /= sizeof( GLuint );
    glDrawElements( this->_drawMode, numberOfIndices, GL_UNSIGNED_INT, 0 );

    return true;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl