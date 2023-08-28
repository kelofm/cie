#ifndef CIE_GL__PART_SCENE_HPP
#define CIE_GL__PART_SCENE_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- Internal Includes ---
#include "packages/scene/inc/Scene.hpp"
#include "packages/file/inc/Part.hpp"

// --- STL Includes ---
#include <list>
#include <memory>


namespace cie::gl {


/*! Interface Scene for drawing Part objects
 * @note: shaders require the following uniforms to be defined:
 *  - 'transformation'  : mat4 for the camera transformation matrix
 *  - 'cameraPosition'  : vec3 for the camera position
 */
class PartScene : public Scene
{
public:
    using part_container = std::list<PartPtr>;

public:
    PartScene( const std::string& r_name,
               part_container&& r_parts,
               Shader::SharedPointer p_vertexShader,
               Shader::SharedPointer p_geometryShader,
               Shader::SharedPointer p_fragmentShader,
               utils::Logger& r_logger,
               VertexBuffer::SharedPointer p_vertexBuffer,
               ElementBuffer::SharedPointer p_elementBuffer );

    void addPart( PartPtr p_part );
    void removePart( PartPtr p_part );

    const part_container& parts() const;

    /// Rewrite buffers if the parts were modified
    virtual void updateParts();

    virtual bool update_impl() override = 0;

private:
    PartScene() = delete;
    PartScene( const PartScene& r_rhs ) = delete;
    PartScene( PartScene&& r_rhs ) = delete;
    PartScene& operator=( const PartScene& r_rhs ) = delete;

protected:
    part_container _parts;
    bool           _partsRequireUpdate;
};


} // namespace cie::gl

#endif