#define _USE_MATH_DEFINES

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/scene/inc/Scene.hpp"
#include "packages/context/inc/GLFWContext.hpp"
#include "packages/context/inc/GLFWWindow.hpp"
#include "packages/context/inc/GLFWMonitor.hpp"
#include "packages/shaders/inc/Shader.hpp"
#include "packages/buffer/inc/BufferManager.hpp"
#include "packages/camera/inc/PerspectiveCamera.hpp"
#include "packages/utility/inc/paths.hpp"

// --- STL Includes ---
#include <chrono>
#include <thread>
#include <math.h>


namespace cie::gl {


class TestScene : public Scene
{
public:
    template <class ...Args>
    TestScene(Args&&... args) : Scene(std::forward<Args>(args)...)
    {
        glEnable(GL_DEPTH_TEST);
    }

private:
    bool update_impl() override
    {
        GLint64 numberOfIndices;
        glGetBufferParameteri64v(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &numberOfIndices);
        numberOfIndices /= sizeof(GLuint);
        glDrawElements(GL_TRIANGLES, numberOfIndices, GL_UNSIGNED_INT, 0);

        return true;
    }
};


CIE_TEST_CASE("Scene", "[scene]")
{
    CIE_TEST_CASE_INIT("Scene")

    // Context
    auto p_context = gl::GLFWContextSingleton::get();

    // Window
    AbsWindow::SharedPointer p_window;
    CIE_CHECK_POINTER(p_context.lock())
    CIE_TEST_REQUIRE_NOTHROW(p_window = p_context.lock()->newWindow());

    // Shaders
    const std::filesystem::path shaderDir = GL_SHADER_SOURCE_PATH;

    auto shaderPaths = [&shaderDir](const std::filesystem::path& r_shaderName)
    {
        return std::pair<std::filesystem::path,std::filesystem::path>(
            (shaderDir / r_shaderName) += ".json",
            (shaderDir / r_shaderName) += ".glsl"
       );
    };

    auto p_vertexShader = makeVertexShader(
        shaderPaths("default/vertexShader").first,
        shaderPaths("default/vertexShader").second
   );

    auto p_geometryShader = makeGeometryShader(
        shaderPaths("default/geometryShader").first,
        shaderPaths("default/geometryShader").second
   );

    auto p_fragmentShader = makeFragmentShader(
        shaderPaths("default/fragmentShader").first,
        shaderPaths("default/fragmentShader").second
   );

    // Scene
    ScenePtr p_scene;
    BufferManager::SharedPointer p_bufferManager;

    {
        CIE_TEST_CASE_INIT("create scene")
        auto localBlock = p_context.lock()->logger().newBlock("basic scene");

        CIE_TEST_REQUIRE_NOTHROW(
            p_scene = p_window->makeScene<TestScene>(
                "TestScene",
                p_vertexShader,
                p_geometryShader,
                p_fragmentShader)
        );

        VertexBuffer::data_container_type components
        {
            0.0, 0.5, 0.0,
            0.5, -0.5, 0.0,
            -0.5, -0.5, 0.0
        };

        ElementBuffer::data_container_type triangles
        {
            0, 1, 2
        };

        CIE_TEST_REQUIRE_NOTHROW(p_bufferManager = p_scene->bufferManager());
        p_bufferManager->template writeToBoundBuffer<VertexBuffer>(
            components.begin(),
            components.size()
        );
        p_bufferManager->template writeToBoundBuffer<ElementBuffer>(
            triangles.begin(),
            triangles.size()
        );

        CIE_TEST_CHECK_NOTHROW(p_window->update());

        CIE_TEST_CHECK_NOTHROW(p_window->removeScene(p_scene));
    }

    {
        CIE_TEST_CASE_INIT("multiple scenes")
        auto localBlock = p_context.lock()->logger().newBlock("two scenes");

        ScenePtr p_secondScene;

        CIE_TEST_CHECK_NOTHROW(
            p_scene = p_window->makeScene<TestScene>(
                "TestScene0",
                p_vertexShader,
                p_geometryShader,
                p_fragmentShader
           )
       );

        CIE_TEST_CHECK_NOTHROW(
            p_secondScene = p_window->makeScene<TestScene>(
                "TestScene1",
                p_vertexShader,
                p_geometryShader,
                p_fragmentShader
           )
       );

        VertexBuffer::data_container_type components0
        {
            0.0, 0.5, 0.0,
            0.5, -0.5, 0.0,
            0.0, -0.5, 0.0
        };

        ElementBuffer::data_container_type triangles0
        {
            0, 1, 2
        };

        VertexBuffer::data_container_type components1
        {
            0.0, 0.5, 0.0,
            0.5, 0.5, 0.0,
            0.5, -0.5, 0.0
        };

        ElementBuffer::data_container_type triangles1
        {
            0, 1, 2
        };

        p_scene->bufferManager()->template writeToBoundBuffer<VertexBuffer>(
            components0.begin(),
            components0.size()
        );
        p_scene->bufferManager()->template writeToBoundBuffer<ElementBuffer>(
            triangles0.begin(),
            triangles0.size()
        );

        p_secondScene->bufferManager()->template writeToBoundBuffer<VertexBuffer>(
            components1.begin(),
            components1.size()
        );
        p_secondScene->bufferManager()->template writeToBoundBuffer<ElementBuffer>(
            triangles1.begin(),
            triangles1.size()
        );

        CIE_TEST_CHECK_NOTHROW(p_window->update());

        CIE_TEST_CHECK_NOTHROW(p_window->removeScene(p_scene));
        CIE_TEST_CHECK_NOTHROW(p_window->removeScene(p_secondScene));
    }



    // Scene with uniforms
    p_vertexShader = makeVertexShader(
        shaderPaths("trianglesWithNormals/vertexShader").first,
        shaderPaths("trianglesWithNormals/vertexShader").second
   );

    p_geometryShader = makeGeometryShader(
        shaderPaths("trianglesWithNormals/geometryShader").first,
        shaderPaths("trianglesWithNormals/geometryShader").second
   );

    p_fragmentShader = makeFragmentShader(
        shaderPaths("trianglesWithNormals/fragmentShader").first,
        shaderPaths("trianglesWithNormals/fragmentShader").second
   );

    auto p_wireframeVertexShader = makeVertexShader(
        shaderPaths("wireframeWithNormal/vertexShader").first,
        shaderPaths("wireframeWithNormal/vertexShader").second
   );

    auto p_wireframeGeometryShader = makeGeometryShader(
        shaderPaths("wireframeWithNormal/geometryShader").first,
        shaderPaths("wireframeWithNormal/geometryShader").second
   );

    auto p_wireframeFragmentShader = makeFragmentShader(
        shaderPaths("wireframeWithNormal/fragmentShader").first,
        shaderPaths("wireframeWithNormal/fragmentShader").second
   );


    {
        CIE_TEST_CASE_INIT("scene with uniforms")
        auto localBlock = p_context.lock()->logger().newBlock("scene with uniforms");

        // Primary scene for triangles
        CIE_TEST_REQUIRE_NOTHROW(
            p_scene = p_window->makeScene<TestScene>(
                "Scene_triangles",
                p_vertexShader,
                p_geometryShader,
                p_fragmentShader)
       );

        // Secondary scene for the wireframe
        ScenePtr p_secondScene;
        CIE_TEST_REQUIRE_NOTHROW(
            p_secondScene = p_window->makeScene<TestScene>(
                "Scene_wireframe"
                ,p_wireframeVertexShader
                ,p_wireframeGeometryShader
                ,p_wireframeFragmentShader
                ,p_context.lock()->logger()
                ,p_scene->bufferManager()->getBoundBuffer<VertexBuffer>().lock()
                ,p_scene->bufferManager()->getBoundBuffer<ElementBuffer>().lock()
               )
       );

        CIE_TEST_CHECK(p_scene->bufferManager()->template getBoundBuffer<VertexBuffer>().lock() == p_secondScene->bufferManager()->template getBoundBuffer<VertexBuffer>().lock());
        CIE_TEST_CHECK(p_scene->bufferManager()->template getBoundBuffer<ElementBuffer>().lock() == p_secondScene->bufferManager()->template getBoundBuffer<ElementBuffer>().lock());

        const float a = 0.5;

        typename VertexBuffer::data_container_type vertexData
        {
            // position ----------- normal
            a, 0.0, 0.0,    1.0, 0.0, 0.0,  // 1
            a, a, 0.0,      1.0, 0.0, 0.0,  // 3
            a, a, a,        1.0, 0.0, 0.0,  // 7
            a, a, a,        1.0, 0.0, 0.0,  // 7
            a, 0.0, a,      1.0, 0.0, 0.0,  // 5
            a, 0.0, 0.0,    1.0, 0.0, 0.0,  // 1

            0.0, 0.0, 0.0,  -1.0, 0.0, 0.0, // 0
            0.0, 0.0, a,    -1.0, 0.0, 0.0, // 4
            0.0, a, a,      -1.0, 0.0, 0.0, // 6
            0.0, a, a,      -1.0, 0.0, 0.0, // 6
            0.0, a, 0.0,    -1.0, 0.0, 0.0, // 2
            0.0, 0.0, 0.0,  -1.0, 0.0, 0.0, // 0

            0.0, a, 0.0,    0.0, 1.0, 0.0,  // 2
            0.0, a, a,      0.0, 1.0, 0.0,  // 6
            a, a, a,        0.0, 1.0, 0.0,  // 7
            a, a, a,        0.0, 1.0, 0.0,  // 7
            a, a, 0.0,      0.0, 1.0, 0.0,  // 3
            0.0, a, 0.0,    0.0, 1.0, 0.0,  // 2

            0.0, 0.0, 0.0,  0.0, -1.0, 0.0, // 0
            a, 0.0, 0.0,    0.0, -1.0, 0.0, // 1
            a, 0.0, a,      0.0, -1.0, 0.0, // 5
            a, 0.0, a,      0.0, -1.0, 0.0, // 5
            0.0, 0.0, a,    0.0, -1.0, 0.0, // 4
            0.0, 0.0, 0.0,  0.0, -1.0, 0.0, // 0

            0.0, 0.0, a,    0.0, 0.0, 1.0,  // 4
            a, 0.0, a,      0.0, 0.0, 1.0,  // 5
            a, a, a,        0.0, 0.0, 1.0,  // 7
            a, a, a,        0.0, 0.0, 1.0,  // 7
            0.0, a, a,      0.0, 0.0, 1.0,  // 6
            0.0, 0.0, a,    0.0, 0.0, 1.0,  // 4

            0.0, 0.0, 0.0,  0.0, 0.0, -1.0, // 0
            0.0, a, 0.0,    0.0, 0.0, -1.0, // 2
            a, a, 0.0,      0.0, 0.0, -1.0, // 3
            a, a, 0.0,      0.0, 0.0, -1.0, // 3
            a, 0.0, 0.0,    0.0, 0.0, -1.0, // 1
            0.0, 0.0, 0.0,  0.0, 0.0, -1.0  // 0
        };

        // Cube triangles
        ElementBuffer::data_container_type triangles
        {
            0, 1, 2,        3, 4, 5,
            6, 7, 8,        9, 10, 11,
            12, 13, 14,     15, 16, 17,
            18, 19, 20,     21, 22, 23,
            24, 25, 26,     27, 28, 29,
            30, 31, 32,     33, 34, 35
        };

        CIE_TEST_REQUIRE_NOTHROW(p_bufferManager = p_scene->bufferManager());
        p_bufferManager->template writeToBoundBuffer<VertexBuffer>(
            vertexData.begin(),
            vertexData.size()
        );
        p_bufferManager->template writeToBoundBuffer<ElementBuffer>(
            triangles.begin(),
            triangles.size()
        );

        using TCamera = PerspectiveCamera;
        auto p_camera = p_scene->makeCamera<TCamera>();
        p_secondScene->addCamera(p_camera);

        p_camera->setAspectRatio(p_window->getSize().first / double(p_window->getSize().second));
        p_camera->setPosition({a/2.0, a/2.0, 2.0*a + p_camera->clippingPlanes().first});
        p_camera->setFieldOfView(60.0 * M_PI / 180.0);

        TCamera::vector_type center { a/2.0, a/2.0, a/2.0 };

        CIE_TEST_REQUIRE_NOTHROW(p_scene->bindUniform("transformation", p_camera->transformationMatrix()));
        CIE_TEST_REQUIRE_NOTHROW(p_scene->bindUniform("cameraPosition", p_camera->position()));

        GLfloat wireframeScale = 1e-4;
        CIE_TEST_REQUIRE_NOTHROW(p_secondScene->bindUniform("transformation", p_camera->transformationMatrix()));
        CIE_TEST_REQUIRE_NOTHROW(p_secondScene->bindUniform("wireframeScale", wireframeScale));

        const int steps = 360;
        std::chrono::microseconds delay(100);

        // Translate
        for (int i=0; i<steps; ++i)
        {
            p_camera->translate({ i<steps/2 ? 0.01 : -0.01 ,0.0,0.0});
            CIE_TEST_CHECK_NOTHROW(p_window->update());
            std::this_thread::sleep_for(delay);
        }

        // Rotate
        for (int i=0; i<steps; ++i)
        {
            p_camera->rotate(M_PI / 180.0,
                              { 0.0, 0.0, 1.0 },
                              center);
            CIE_TEST_CHECK_NOTHROW(p_window->update());
            std::this_thread::sleep_for(delay);
        }

        // Rotate
        for (int i=0; i<steps; ++i)
        {
            p_camera->rotate((i<steps/2 ? 1.0 : -1.0) * M_PI / 180.0,
                              { 1.0, 0.0, 0.0 },
                              center);
            CIE_TEST_CHECK_NOTHROW(p_window->update());
            std::this_thread::sleep_for(delay);
        }

        // Zoom
        for (int i=0; i<steps; ++i)
        {
            p_camera->zoom(i<steps/2 ? 1.01 : 0.99);
            CIE_TEST_CHECK_NOTHROW(p_window->update());
            std::this_thread::sleep_for(delay);
        }

        CIE_TEST_CHECK_NOTHROW(p_window->removeScene(p_scene));
        CIE_TEST_CHECK_NOTHROW(p_window->removeScene(p_secondScene));
    }
}


} // namespace cie::gl