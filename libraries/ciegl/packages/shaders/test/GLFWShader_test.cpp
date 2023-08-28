// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/shaders/inc/Shader.hpp"
#include "packages/context/inc/GLFWContext.hpp"
#include "packages/context/inc/GLFWWindow.hpp"
#include "packages/context/inc/GLFWMonitor.hpp"
#include "packages/utility/inc/paths.hpp"


namespace cie::gl {


CIE_TEST_CASE( "Shader", "[shader]" )
{
    CIE_TEST_CASE_INIT( "Shader" )

    // Context
    auto p_context = GLFWContextSingleton::get();

    // Window
    AbsWindow::SharedPointer p_window;
    CIE_CHECK_POINTER(p_context.lock())
    CIE_TEST_REQUIRE_NOTHROW( p_window = p_context.lock()->newWindow() );

    // Shader setup
    const std::filesystem::path shaderDir   = GL_SHADER_SOURCE_PATH / "default";

    auto shaderPaths = [&shaderDir]( const std::filesystem::path& r_shaderName )
    {
        return std::pair<std::filesystem::path,std::filesystem::path>(
            (shaderDir / r_shaderName) += ".json",
            (shaderDir / r_shaderName) += ".glsl"
        );
    };

    {
        CIE_TEST_CASE_INIT( "default vertex shader" )

        Shader::SharedPointer p_shader;
        CIE_TEST_CHECK_NOTHROW( p_shader = makeVertexShader(shaderPaths("vertexShader").first,
                                                            shaderPaths("vertexShader").second ) );

        CIE_TEST_REQUIRE( !p_shader->attributes().empty() );
        CIE_TEST_CHECK( p_shader->attributes().size() == 1 );

        const auto& r_attribute = p_shader->attributes()[0];
        CIE_TEST_CHECK( r_attribute.name() == "position" );
        CIE_TEST_CHECK( r_attribute.size() == 3 );
        CIE_TEST_CHECK( r_attribute.stride() == 0 );
        CIE_TEST_CHECK( r_attribute.offset() == 0 );

        CIE_TEST_CHECK( p_shader->uniforms().empty() );
        CIE_TEST_CHECK( p_shader->textures().empty() );
        CIE_TEST_CHECK( p_shader->outputs().empty() );
    }

    {
        CIE_TEST_CASE_INIT( "default geometry shader" )

        Shader::SharedPointer p_shader;
        CIE_TEST_CHECK_NOTHROW( p_shader = makeGeometryShader(shaderPaths("geometryShader").first,
                                                              shaderPaths("geometryShader").second ) );

        CIE_TEST_CHECK( p_shader->attributes().empty() );
        CIE_TEST_CHECK( p_shader->uniforms().empty() );
        CIE_TEST_CHECK( p_shader->textures().empty() );
        CIE_TEST_CHECK( p_shader->outputs().empty() );
    }

    {
        CIE_TEST_CASE_INIT( "default fragment shader" )

        Shader::SharedPointer p_shader;
        CIE_TEST_CHECK_NOTHROW( p_shader = makeFragmentShader(shaderPaths("fragmentShader").first,
                                                              shaderPaths("fragmentShader").second ) );

        CIE_TEST_REQUIRE( !p_shader->outputs().empty() );
        CIE_TEST_CHECK( p_shader->outputs().size() == 1 );

        const auto& r_output = p_shader->outputs()[0];
        CIE_TEST_CHECK( r_output.name() == "color" );

        CIE_TEST_CHECK( p_shader->attributes().empty() );
        CIE_TEST_CHECK( p_shader->uniforms().empty() );
        CIE_TEST_CHECK( p_shader->textures().empty() );
    }
}


} // namespace cie::gl