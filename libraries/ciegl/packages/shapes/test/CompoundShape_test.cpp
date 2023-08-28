// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/context/inc/GLFWWindow.hpp"
#include "packages/context/inc/GLFWContext.hpp"
#include "packages/shapes/inc/Arrow.hpp"
#include "packages/shapes/inc/CompoundShape.hpp"
#include "packages/camera/inc/OrthographicCamera.hpp"
#include "packages/buffer/inc/Vertex3.hpp"
#include "packages/file/inc/GenericPart.hpp"
#include "packages/scene/inc/GenericPartScene.hpp"
#include "packages/utility/inc/paths.hpp"

// --- STL Includes ---
#include <chrono>
#include <thread>


namespace cie::gl {


CIE_TEST_CASE( "CompoundShape", "[shapes]" )
{
    CIE_TEST_CASE_INIT( "CompoundShape" )

    auto p_context = GLFWContextSingleton::get();
    CIE_CHECK_POINTER(p_context.lock())
    auto p_window  = p_context.lock()->newWindow();

    auto p_scene   = p_window->makeScene<GenericPartScene>(
        "ArrowScene",
        GL_TRIANGLES,
        GL_SHADER_SOURCE_PATH / "Triangulated3DPartScene"
    );

    auto p_camera = p_scene->makeCamera<OrthographicCamera>();

    p_camera->setPose( {1.0, 1.0, 0.0},
                       {-1.0, -1.0, 0.0},
                       {0.0, 0.0, 1.0} );
    p_camera->setFieldOfView( 90.0 * 3.14159265 / 180.0 );
    p_camera->setClippingPlanes( 1e-3, 1e2 );

    p_scene->bindUniform( "transformation", p_camera->transformationMatrix() );
    p_scene->bindUniform( "cameraPosition", p_camera->position() );

    using VertexType   = Vertex3;
    using ArrowType    = Arrow<VertexType>;
    using CompoundType = CompoundShape<VertexType>;

    auto p_part = std::make_shared<GenericPart>( 3, 3, 3 );

    auto p_arrow0 = CompoundType::shape_ptr( new ArrowType(
        p_part->attributes(),
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0}
    ) );

    auto p_arrow1 = CompoundType::shape_ptr( new ArrowType(
        p_part->attributes(),
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0}
    ) );

    auto compoundShape = CompoundType( {p_arrow0, p_arrow1} );

    p_part->indices() = compoundShape.indices();

    p_scene->addPart( p_part );
    p_window->update();
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));
}


} // namespace cie::gl