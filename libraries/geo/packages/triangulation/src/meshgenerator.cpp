// --- Internal Includes ---
#include "packages/triangulation/inc/meshgenerator.hpp"
#include "packages/triangulation/inc/meshgenerator_helper.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <cmath>
#include <numeric>
#include <algorithm>

namespace cie::geo {

Triangulation triangulate(
    std::span<const Vertex2D> polygon,
    TriangulationParameters parameters) {
        CIE_BEGIN_EXCEPTION_TRACING
            std::vector<IndexVector> regions{IndexVector(polygon.size())};
            std::iota(regions[0].begin(), regions[0].end(), 0);
            return triangulate(polygon, regions, parameters);
        CIE_END_EXCEPTION_TRACING
}

Triangulation triangulate(
    std::span<const Vertex2D> vertices,
    const std::vector<IndexVector>& polygonRegions,
    TriangulationParameters parameters) {
        CIE_BEGIN_EXCEPTION_TRACING
            using namespace meshgeneratorhelper;

            meshgeneratorhelper::prepareForTriangulating(
                vertices,
                polygonRegions,
                parameters);

            Triangulation triangulation;
            triangulation.first.insert(
                triangulation.first.end(),
                vertices.begin(),
                vertices.end());

            auto regions = polygonRegions; // create copy

            while( checkBreakingCriteria(regions)) {
                for(size_t iRegion = 0; iRegion<regions.size( ); ++iRegion) {
                    while(chopTriangle(
                        triangulation,
                        regions[iRegion],
                        parameters.goodChopRatio));
                    if(!attemptDivision(
                        triangulation.first,
                        regions, iRegion, parameters.divisionAngle, parameters.edgeLength)) {
                            chopTriangle( triangulation, regions[iRegion], 0.0 );
                    }
                }
            } // while checkBreakCriteria

            return triangulation;
        CIE_END_EXCEPTION_TRACING
}


}// namespace cie::geo