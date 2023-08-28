#ifndef CIE_GL_PLOT_PLOT_2_HPP
#define CIE_GL_PLOT_PLOT_2_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- Internal Includes ---
#include "packages/plot/inc/AbsPlot2.hpp"


namespace cie::gl {


/// Default Plot2 that will cause a compilation error
template <class ...Args>
class Plot2 : public AbsPlot2<std::false_type,Args...>
{};


/// Standard x-y plot
template < concepts::NumericContainer XContainer,
           concepts::NumericContainer YContainer >
struct Plot2<XContainer,YContainer> : public AbsPlot2<std::true_type>
{
    Plot2( const XContainer& r_XContainer,
           const YContainer& r_YContainer );
};


/// Convenience function for deduced template arguments
template < concepts::NumericContainer XContainer,
           concepts::NumericContainer YContainer >
void plot2( const XContainer& r_XContainer,
            const YContainer& r_YContainer );


} // namespace cie::gl

#include "packages/plot/impl/Plot2_impl.hpp"

#endif