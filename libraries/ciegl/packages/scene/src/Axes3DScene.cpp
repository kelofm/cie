// --- External Includes ---
#include "glm/glm.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/ranges/inc/TransformView.hpp"

// --- Internal Includes ---
#include "packages/scene/inc/Axes3DScene.hpp"
#include "packages/file/inc/GenericPart.hpp"
#include "cmake_variables.hpp"


namespace cie::gl {


Axes3DScene::Axes3DScene(const std::string& r_name,
                         AbsCamera::SharedPointer p_camera,
                         utils::Logger& r_logger) :
    GenericPartScene(r_name,
                     GL_TRIANGLES,
                     getDataPath() / "ciegl/shaders/coloredAxes",
                     r_logger),
    _p_axes( nullptr ),
    _box( {{0.0, 0.2}, {0.0, 0.2}} ),
    _transformationMatrix()
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER( p_camera )

    // Set camera and bind uniforms
    this->addCamera( p_camera );

    // Get number of vertex attributes
    Size vertexAttributeSize = makeVertex<ColoredVertex3>(
        ColoredVertex3::attribute_container::SharedPointer(new ColoredVertex3::attribute_container)
    )->numberOfAttributes();

    auto p_part = GenericPartPtr(
        new GenericPart( 3, vertexAttributeSize, 3 )
    );

    this->_p_axes = Axes3DScene::axes_ptr( new Axes(
        p_part->attributes(),
        3
    ) );

    p_part->indices() = this->_p_axes->indices();

    // Add the part containing the arrows
    this->addPart( p_part );

    // Bind uniforms to internal camera
    this->bindUniform( "transformation", this->_transformationMatrix );
    this->bindUniform( "cameraPosition", p_camera->position() );

    CIE_END_EXCEPTION_TRACING
}


bool Axes3DScene::update_impl()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto p_camera     = this->camera( 0 );
    auto p_axes       = this->_p_axes;
    const auto& r_box = this->_box;

    Axes::vector_type position {
        ( r_box[0].second + r_box[0].first ) / 2.0,
        ( r_box[1].second + r_box[1].first ) / 2.0,
        0.0
    };

    p_axes->setPose( this->_p_axes->position(),
                     p_camera->direction(),
                     p_camera->up() );

    double aspectRatio = this->camera()->aspectRatio();
    double scaleX      = 1.0 / ( ( r_box[0].second - r_box[0].first ) / 2.0 ) * aspectRatio;
    double scaleY      = 1.0 / ( ( r_box[1].second - r_box[1].first ) / 2.0 );

    position[0] /= aspectRatio;

    this->_transformationMatrix = glm::ortho(
        -position[0] * scaleX,
        (1.0 - position[0]) * scaleX,
        -position[1] * scaleY,
        (1.0 - position[1]) * scaleY
    )
    *
    glm::lookAt(
        this->_p_axes->position(),
        this->_p_axes->position() + this->_p_axes->direction(),
        this->_p_axes->up()
    );

    return GenericPartScene::update_impl();

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl