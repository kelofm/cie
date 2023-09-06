#ifndef CIE_GL_ABS_VERTEX_IMPL_HPP
#define CIE_GL_ABS_VERTEX_IMPL_HPP

// --- Graphics Includes ---
#include "packages/buffer/inc/AbsVertex.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie::gl {


inline AbsVertex::attribute_range AbsVertex::attribute( Size index )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK( index < this->_r_attributeOffsets.size() )

    auto it_begin = this->_p_attributes->begin() + this->_offset + this->_r_attributeOffsets[index];

    Size attributeSize = 0;

    if ( index != this->_r_attributeOffsets.size() - 1 )
        attributeSize = this->_r_attributeOffsets[index+1] - this->_r_attributeOffsets[index];
    else
        attributeSize = this->numberOfAttributes() - this->_r_attributeOffsets.back();

    CIE_OUT_OF_RANGE_CHECK( this->_offset + attributeSize <= this->_p_attributes->size() )

    return AbsVertex::attribute_range( it_begin, it_begin + attributeSize );

    CIE_END_EXCEPTION_TRACING
}


inline AbsVertex::attribute_const_range AbsVertex::attribute( Size index ) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK( index < this->_r_attributeOffsets.size() )

    auto it_begin = this->_p_attributes->begin() + this->_offset + this->_r_attributeOffsets[index];

    Size attributeSize = 0;

    if ( index != this->_r_attributeOffsets.size() - 1 )
        attributeSize = this->_r_attributeOffsets[index+1] - this->_r_attributeOffsets[index];
    else
        attributeSize = this->numberOfAttributes() - this->_r_attributeOffsets.back();

    CIE_OUT_OF_RANGE_CHECK( this->_offset + attributeSize <= this->_p_attributes->size() )

    return AbsVertex::attribute_const_range( it_begin, it_begin + attributeSize );

    CIE_END_EXCEPTION_TRACING
}


inline AbsVertex::value_type& AbsVertex::at( Size attributeIndex,
                                             Size componentIndex )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK( attributeIndex < this->_r_attributeOffsets.size() )

    return this->at(
        this->_r_attributeOffsets[attributeIndex] + componentIndex
    );

    CIE_END_EXCEPTION_TRACING
}


inline AbsVertex::value_type AbsVertex::at(Size attributeIndex,
                                           Size componentIndex) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_OUT_OF_RANGE_CHECK( attributeIndex < this->_r_attributeOffsets.size() )

    return this->at(
        this->_r_attributeOffsets[attributeIndex] + componentIndex
    );

    CIE_END_EXCEPTION_TRACING
}


inline AbsVertex::value_type& AbsVertex::at( Size index )
{
    CIE_OUT_OF_RANGE_CHECK( index < this->numberOfAttributes() )
    Size indexInContainer = this->_offset + index;
    CIE_OUT_OF_RANGE_CHECK( indexInContainer < this->_p_attributes->size() )

    return this->_p_attributes->at( indexInContainer );
}


inline AbsVertex::value_type AbsVertex::at(Size index) const
{
    CIE_OUT_OF_RANGE_CHECK( index < this->numberOfAttributes() )
    Size indexInContainer = this->_offset + index;
    CIE_OUT_OF_RANGE_CHECK( indexInContainer < this->_p_attributes->size() )

    return this->_p_attributes->at( indexInContainer );
}


template <class VertexType, class ...Args>
requires std::derived_from<VertexType,AbsVertex>
std::shared_ptr<VertexType> makeVertex(AttributeContainer::SharedPointer p_attributeContainer,
                                       Args&&... args)
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::shared_ptr<VertexType> p_vertex( new VertexType(
        p_attributeContainer,
        p_attributeContainer->size(),
        std::forward<Args>(args)...
    ));

    return p_vertex;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl


#endif
