#ifndef CIE_GL_SCENE_AXES_3D_SCENE_HPP
#define CIE_GL_SCENE_AXES_3D_SCENE_HPP

// --- Internal Includes ---
#include "packages/scene/inc/GenericPartScene.hpp"
#include "packages/camera/inc/AbsCamera.hpp"
#include "packages/shapes/inc/Axes.hpp"
#include "packages/traits/inc/GLTraits.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <vector>
#include <array>
#include <utility>


namespace cie::gl {


class Axes3DScene :
    public GenericPartScene,
    public GLTraits
{
public:
    using axes_ptr   = std::shared_ptr<Axes>;
    using screen_box = StaticArray<std::pair<double,double>,2>;

public:
    Axes3DScene(const std::string& r_name,
                AbsCamera::SharedPointer p_camera,
                utils::Logger& r_logger = utils::LoggerSingleton::get());

protected:
    virtual bool update_impl() override;

protected:
    axes_ptr                        _p_axes;
    screen_box                      _box;
    AbsCamera::internal_matrix_type _transformationMatrix;
};





} // namespace cie::gl


#endif