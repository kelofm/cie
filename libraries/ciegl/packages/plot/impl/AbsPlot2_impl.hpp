#ifndef CIE_GL_ABS_PLOT_2_IMPL_HPP
#define CIE_GL_ABS_PLOT_2_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/plot/inc/Plot2.hpp"
#include "packages/shaders/inc/Shader.hpp"
#include "packages/utility/inc/GLError.hpp"
#include "packages/camera/inc/camera_utilities.hpp"
#include "packages/utility/inc/paths.hpp"


namespace cie::gl {


/* --- AbsPlot2 --- */

template <class Tag, class ...Args>
requires Tag::value
AbsPlot2<Tag,Args...>::AbsPlot2( AbsWindow::SharedPointer p_window ) :
    AbsPlot( p_window ),
    _vertices(),
    _p_indices( new typename AbsPlot2<Tag,Args...>::index_container ),
    _p_camera( nullptr ),
    _p_scene( nullptr ),
    _p_controls( nullptr )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->initializeScene();

    CIE_END_EXCEPTION_TRACING
}


template <class Tag, class ...Args>
requires Tag::value
AbsPlot2<Tag,Args...>::AbsPlot2() :
    AbsPlot(),
    _vertices(),
    _p_indices( new typename AbsPlot2<Tag,Args...>::index_container ),
    _p_camera( nullptr ),
    _p_scene( nullptr ),
    _p_controls( nullptr )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->initializeScene();

    CIE_END_EXCEPTION_TRACING
}


template <class Tag, class ...Args>
requires Tag::value
void AbsPlot2<Tag,Args...>::fit( bool keepAspectRatio )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Get data domain
    double xMin = std::numeric_limits<double>().max();
    double xMax = std::numeric_limits<double>().min();
    double yMin = std::numeric_limits<double>().max();
    double yMax = std::numeric_limits<double>().min();

    for ( const auto& rp_vertex : this->_vertices )
    {
        if ( rp_vertex->x() < xMin )
            xMin = rp_vertex->x();
        if ( xMax < rp_vertex->x() )
            xMax = rp_vertex->x();

        if ( rp_vertex->y() < yMin )
            yMin = rp_vertex->y();
        if ( yMax < rp_vertex->y() )
            yMax = rp_vertex->y();
    }

    // Adjust aspect ratio
    if ( keepAspectRatio )
    {
        auto maxRange = xMax - xMin;

        if ( yMax - yMin > maxRange )
        {
            maxRange = yMax - yMin;

            auto midPoint = (xMax + xMin) / 2;
            xMax = midPoint + maxRange / 2;
            xMin = midPoint - maxRange / 2;
        }
        else
        {
            auto midPoint = (yMax + yMin) / 2;
            yMax = midPoint + maxRange / 2;
            yMin = midPoint - maxRange / 2;
        }

    }

    // Scale data
    auto makeNormalizeFunctor = []( double minValue, double maxValue )
    {
        CIE_CHECK( minValue < maxValue, "Invalid domain" )
        double offset = minValue;
        double scale  = maxValue - minValue;

        return [offset,scale]( typename AbsPlot2<Tag,Args...>::vertex_type::value_type& r_component ) -> void
        {
            r_component -= offset;
            r_component /= scale;
        };
    };

    auto normalizeX = makeNormalizeFunctor( xMin, xMax );
    auto normalizeY = makeNormalizeFunctor( yMin, yMax );

    this->transform( normalizeX, normalizeY );

    this->_p_scene->setUpdateFlag();

    // Adjust camera
    using BoxType = geo::Box<3>;
    BoxType box( BoxType::Point { 0.0, 0.0, -0.5 },
                 BoxType::Point { 1.0, 1.0, 1.0 } );

    this->_p_camera->setPosition( { 0.5, 0.5, 2.0 } );
    this->_p_camera->setClippingPlanes( 0.1, 4.0 );

    fitViewToBox(*this->_p_camera, box);

    // Adjust controls
    this->_p_controls->setMovementScale( 0.1 );

    CIE_END_EXCEPTION_TRACING
}


template <class Tag, class ...Args>
requires Tag::value
void AbsPlot2<Tag,Args...>::initializeScene()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto p_scene = this->_p_window->makeScene<typename AbsPlot2<Tag,Args...>::Plot2Scene>( "Plot2Scene",
                                                                                           this->_p_attributes,
                                                                                           this->_p_indices );

    this->_p_scene = p_scene;

    this->_p_camera = p_scene->getCamera();
    this->_p_camera->setPose( { 0.0, 0.0, 1.0 },
                              { 0.0, 0.0, -1.0 },
                              { 0.0, 1.0, 0.0 } );
    this->_p_camera->setClippingPlanes( 0.5, 1.5 );

    this->_p_camera->setAspectRatio(_p_window->getAspectRatio());

    this->_p_controls = typename AbsPlot2<Tag,Args...>::controls_ptr(
        new typename AbsPlot2<Tag,Args...>::controls_type(true)
    );
    this->_p_controls->bind( this->_p_window, this->_p_camera );

    CIE_END_EXCEPTION_TRACING
}


template <class Tag, class ...Args>
requires Tag::value
void AbsPlot2<Tag,Args...>::transform( typename AbsPlot2<Tag,Args...>::transform_functor xTransform,
                                       typename AbsPlot2<Tag,Args...>::transform_functor yTransform )
{
    CIE_BEGIN_EXCEPTION_TRACING

    for ( auto& rp_vertex : this->_vertices )
    {
        xTransform( rp_vertex->x() );
        yTransform( rp_vertex->y() );
    }

    CIE_END_EXCEPTION_TRACING
}


/* --- AbsPlot2<Tag,Args...>::Plot2Scene --- */

template <class Tag, class ...Args>
requires Tag::value
const std::filesystem::path
AbsPlot2<Tag,Args...>::Plot2Scene::_shaderDirectory = GL_SHADER_INSTALL_PATH / "plot2";


template <class Tag, class ...Args>
requires Tag::value
AbsPlot2<Tag,Args...>::Plot2Scene::Plot2Scene(const std::string& r_name,
                                              AttributeContainer::SharedPointer p_attributes,
                                              typename AbsPlot2<Tag,Args...>::index_container_ptr p_indices) :
    Scene(
        r_name,
        makeVertexShader(this->_shaderDirectory / "vertexShader.json",
                         this->_shaderDirectory / "vertexShader.glsl"),
        nullptr,
        makeFragmentShader(this->_shaderDirectory / "fragmentShader.json",
                           this->_shaderDirectory / "fragmentShader.glsl"),
        utils::LoggerSingleton::get()
    ),
    _p_attributes( p_attributes ),
    _p_indices( p_indices ),
    _updateFlag( true ),
    _drawMode( GL_LINE_STRIP )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto p_camera = this->makeCamera<typename AbsPlot2<Tag,Args...>::camera_type>();
    this->bindUniform( "transformation", p_camera->transformationMatrix() );

    CIE_END_EXCEPTION_TRACING
}


template <class Tag, class ...Args>
requires Tag::value
CameraPtr AbsPlot2<Tag,Args...>::Plot2Scene::getCamera()
{
    return *this->_cameras.begin();
}


template <class Tag, class ...Args>
requires Tag::value
void AbsPlot2<Tag,Args...>::Plot2Scene::setUpdateFlag()
{
    this->_updateFlag = true;
}


template <class Tag, class ...Args>
requires Tag::value
void AbsPlot2<Tag,Args...>::Plot2Scene::setDrawMode( GLenum drawMode )
{
    this->_drawMode = drawMode;
}


template <class Tag, class ...Args>
requires Tag::value
bool AbsPlot2<Tag,Args...>::Plot2Scene::update_impl()
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto& r_bufferManager = *this->_p_bufferManager;

    // Write attributes into the buffer if requested
    if ( this->_updateFlag )
    {
        r_bufferManager.template writeToBoundBuffer<VertexBuffer>(
            this->_p_attributes->begin(),
            this->_p_attributes->size()
        );
        r_bufferManager.template writeToBoundBuffer<ElementBuffer>(
            this->_p_indices->begin(),
            this->_p_indices->size()
        );
        this->_updateFlag = false;

        return true;
    }

    checkGLErrors( *this, "Error writing to buffer" );

    GLint64 numberOfIndices;
    glGetBufferParameteri64v( GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &numberOfIndices );
    numberOfIndices /= sizeof( ElementBuffer::data_type );

    glDrawElements( this->_drawMode, numberOfIndices, GL_UNSIGNED_INT, 0 );

    return true;

    CIE_END_EXCEPTION_TRACING
}




} // namespace cie::gl

#endif