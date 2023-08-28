#ifndef CIE_GL_SCENE_GENERIC_PART_SCENE_HPP
#define CIE_GL_SCENE_GENERIC_PART_SCENE_HPP

// --- Internal Includes ---
#include "packages/scene/inc/PartScene.hpp"

// --- STL Includes ---
#include <filesystem>


namespace cie::gl {


/// Convenience Scene class
class GenericPartScene : public PartScene
{
public:
    GenericPartScene( const std::string& r_name,
                      Size drawMode,
                      GenericPartScene::part_container&& r_parts,
                      Shader::SharedPointer p_vertexShader,
                      Shader::SharedPointer p_geometryShader,
                      Shader::SharedPointer p_fragmentShader,
                      utils::Logger& r_logger,
                      VertexBuffer::SharedPointer p_vertexBuffer,
                      ElementBuffer::SharedPointer p_elementBuffer );

    GenericPartScene( const std::string& r_name,
                      Size drawMode,
                      const std::filesystem::path& r_shaderDirectory,
                      utils::Logger& r_logger = utils::LoggerSingleton::get() );

    void setDrawMode( Size drawMode );

protected:
    virtual bool update_impl() override;

protected:
    Size _drawMode;
};


} // namespace cie::gl


#endif