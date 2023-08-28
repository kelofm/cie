// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- GEO Includes ---
#include "packages/primitives/inc/Sphere.hpp"

// --- Internal Includes ---
#include "packages/plot/inc/Plot2_MarchingCubes.hpp"

// --- STL Includes ---
#include <chrono>
#include <thread>


namespace cie::gl {


CIE_TEST_CASE( "Plot2 - MarchingCubes", "[plot]" )
{
    CIE_TEST_CASE_INIT( "Plot2 - MarchingCubes" )

    // Construct marching cubes
    const Size Dimension = 2;
    using CoordinateType = Double;
    using TargetType     = geo::Object<Dimension,Bool,CoordinateType>;
    using MarchingCubes  = geo::StructuredMarchingCubes<TargetType>;

    auto p_target = MarchingCubes::target_ptr(
        new geo::boolean::Sphere<Dimension,CoordinateType>( {0.0,0.0}, 1.0 )
    );

    MarchingCubes::Domain domain {{ {-1.0,1.0}, {-1.0,1.0} }};
    MarchingCubes::Resolution numberOfPoints { 25, 25 };

    MarchingCubes marchingCubes( p_target,
                                 domain,
                                 numberOfPoints,
                                 nullptr );

    Plot2<MarchingCubes> plot( marchingCubes );
    plot.update();
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    //plot2( marchingCubes );
}


} // namespace cie::gl