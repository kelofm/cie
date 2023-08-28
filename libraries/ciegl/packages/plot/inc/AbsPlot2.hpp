#ifndef CIE_GL_ABS_PLOT_2_HPP
#define CIE_GL_ABS_PLOT_2_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- Internal Includes ---
#include "packages/plot/inc/AbsPlot.hpp"
#include "packages/buffer/inc/Vertex2.hpp"
#include "packages/scene/inc/Scene.hpp"

#include "packages/camera/inc/OrthographicCamera.hpp"
#include "packages/control/inc/PanZoomCameraControls.hpp"

// --- STL Includes ---
#include <vector>
#include <memory>


namespace cie::gl {


template <class Tag, class ...Args>
requires Tag::value
class AbsPlot2 : public AbsPlot
{
public:
    using vertex_type         = Vertex2;
    using vertex_container    = std::vector<vertex_type::SharedPointer>;

    using index_container     = ElementBuffer::data_container_type;
    using index_container_ptr = std::shared_ptr<index_container>;

    using camera_type         = OrthographicCamera;
    using controls_type       = PanZoomCameraControls;
    using controls_ptr        = std::shared_ptr<controls_type>;

protected:
    using transform_functor   = std::function<void(typename AbsPlot2<Tag,Args...>::vertex_type::value_type&)>;

public:
    AbsPlot2( AbsWindow::SharedPointer p_window );
    AbsPlot2();

    void fit( bool keepAspectRatio = false ) override;

private:
    void initializeScene();

    void transform( transform_functor xTransform,
                    transform_functor yTransform );

protected:
    /// Internal Scene class used by AbsPlot2
    class Plot2Scene : public Scene
    {
    public:
        Plot2Scene(const std::string& r_name,
                   AttributeContainer::SharedPointer p_attributes,
                   index_container_ptr p_indices);

        AbsCamera::SharedPointer getCamera();

        /// Sets a flag that forces all attributes to be rebuffered on the next update
        void setUpdateFlag();

        /// Set draw mode
        void setDrawMode( GLenum drawMode );

    private:
        bool update_impl() override;

    private:
        AttributeContainer::SharedPointer  _p_attributes;

        index_container_ptr                _p_indices;

        index_container                    _indices;

        bool                               _updateFlag;

        GLenum                             _drawMode;

        static const std::filesystem::path _shaderDirectory;
    };

protected:
    vertex_container            _vertices;

    index_container_ptr         _p_indices;

    AbsCamera::SharedPointer    _p_camera;

    std::shared_ptr<Plot2Scene> _p_scene;

    controls_ptr                _p_controls;
};


} // namespace cie::gl

#include "packages/plot/impl/AbsPlot2_impl.hpp"

#endif