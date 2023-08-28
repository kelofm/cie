// --- External Includes ---
#include "glm/glm.hpp"

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/texture/inc/Texture.hpp"
#include "packages/context/inc/GLFWContext.hpp"
#include "packages/utility/inc/paths.hpp"
#include "packages/shaders/inc/Shader.hpp"

// --- STL Includes ---
#include <chrono>
#include <thread>


namespace cie::gl {


class TextureTestScene : public Scene
{
public:
    template <class ...Args>
    TextureTestScene( Args&&... args ) : Scene( std::forward<Args>(args)... ) {}

private:
    bool update_impl() override
    {
        GLint64 numberOfIndices;
        glGetBufferParameteri64v( GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &numberOfIndices );
        numberOfIndices /= sizeof( GLuint );
        glDrawElements( GL_TRIANGLES, numberOfIndices, GL_UNSIGNED_INT, 0 );

        return true;
    }
};



CIE_TEST_CASE( "Texture", "[texture]" )
{
    CIE_TEST_CASE_INIT( "Texture" )

    auto p_context = gl::GLFWContextSingleton::get();
    CIE_CHECK_POINTER(p_context.lock())
    auto p_window  = p_context.lock()->newWindow();

    auto shaderPath = GL_SHADER_SOURCE_PATH / "defaultWithTexture";

    auto p_vertexShader = makeVertexShader(
        shaderPath / "vertexShader.json",
        shaderPath / "vertexShader.glsl"
    );

    auto p_geometryShader = makeGeometryShader(
        shaderPath / "geometryShader.json",
        shaderPath / "geometryShader.glsl"
    );

    auto p_fragmentShader = makeFragmentShader(
        shaderPath / "fragmentShader.json",
        shaderPath / "fragmentShader.glsl"
    );

    auto p_scene = p_window->makeScene<TextureTestScene>(
        "TextureTestScene",
        p_vertexShader,
        p_geometryShader,
        p_fragmentShader
    );

    glm::mat4 transformationMatrix( 1.0f );
    Image image( GL_TEXTURE_SOURCE_PATH / "cage.png");

    p_scene->bindUniform( "transformation", transformationMatrix );
    p_scene->loadTexture( "image", image );

    const VertexBuffer::data_type tmp = 0.5;
    VertexBuffer::data_container_type components {
        -tmp, -tmp, -tmp, 0.0, 0.0,
        tmp, -tmp, -tmp, 1.0, 0.0,
        -tmp, tmp, -tmp, 0.0, 1.0,
        tmp, tmp, -tmp, 1.0, 1.0
    };

    ElementBuffer::data_container_type triangles {
        0, 1, 3,
        0, 3, 2
    };

    p_scene->bufferManager()->template writeToBoundBuffer<VertexBuffer>(
        components.begin(),
        components.size()
    );
    p_scene->bufferManager()->template writeToBoundBuffer<ElementBuffer>(
        triangles.begin(),
        triangles.size()
    );

    p_window->update();
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));
}


} // namespace cie::gl