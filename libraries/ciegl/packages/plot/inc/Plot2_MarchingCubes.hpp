#ifndef CIE_GL_PLOT_2_MARCHING_CUBES_HPP
#define CIE_GL_PLOT_2_MARCHING_CUBES_HPP

// --- Mesh Includes ---
#include "packages/marchingprimitives/inc/MarchingCubes.hpp"
#include "packages/marchingprimitives/inc/concepts.hpp"

// --- Internal Includes ---
#include "packages/plot/inc/Plot2.hpp"


namespace cie::gl {


template <concepts::MarchingPrimitives TObject>
requires concepts::Cube<typename TObject::primitive_type>
struct Plot2<TObject> : public AbsPlot2<std::true_type>
{
    Plot2( TObject& r_marchingCubes );
};


template <concepts::MarchingPrimitives TObject>
requires concepts::Cube<typename TObject::primitive_type>
void plot2( TObject& r_marchingCubes );


} // namespace cie::gl

#include "packages/plot/impl/Plot2_MarchingCubes_impl.hpp"

#endif