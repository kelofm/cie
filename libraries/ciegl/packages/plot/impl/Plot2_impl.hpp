#ifndef CIE_GL_PLOT_PLOT_2_IMPL_HPP
#define CIE_GL_PLOT_PLOT_2_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/stl_extension/inc/resize.hpp"

// --- STL Includes ---
#include <limits>


namespace cie::gl {


template < concepts::NumericContainer XContainer,
           concepts::NumericContainer YContainer >
Plot2<XContainer,YContainer>::Plot2(const XContainer& r_xContainer,
                                     const YContainer& r_yContainer) :
    AbsPlot2<std::true_type>()
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Check data sizes
    CIE_CHECK(
        r_xContainer.size() == r_yContainer.size(),
        "Container size mismatch! "
            + std::to_string(r_xContainer.size()) + " != "
            + std::to_string(r_yContainer.size())
    )

    CIE_CHECK(
        r_xContainer.size() > 0,
        "Cannot plot empty data"
    )

    using VertexType = typename Plot2<XContainer,YContainer>::vertex_type;

    // Initialize plot
    this->_vertices.reserve(r_xContainer.size());
    this->_p_indices->reserve(r_xContainer.size());
    ElementBuffer::data_type indexCounter = 0;

    {
        auto p_vertex = makeVertex<VertexType>(
            this->_p_attributes,
            r_xContainer[0],
            r_yContainer[0]
        );
        this->_p_attributes->reserve(r_xContainer.size() * p_vertex->numberOfAttributes());
        this->_vertices.push_back(p_vertex);
        this->_p_indices->push_back(indexCounter++);
    }

    // Load data
    for (; indexCounter<r_xContainer.size(); ++indexCounter)
    {
        this->_vertices.push_back(makeVertex<VertexType>(
            this->_p_attributes,
            r_xContainer[indexCounter],
            r_yContainer[indexCounter]
        ));

        this->_p_indices->push_back(indexCounter);
    }

    // Adjust camera
    this->fit();

    CIE_END_EXCEPTION_TRACING
}


template < concepts::NumericContainer XContainer,
           concepts::NumericContainer YContainer >
void plot2(const XContainer& r_XContainer,
            const YContainer& r_YContainer)
{
    Plot2<XContainer,YContainer>(r_XContainer, r_YContainer).show();
}


} // namespace cie::gl


#endif