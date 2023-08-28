#ifndef CIE_GL_TEXTURE_TEXTURE_IMPL_HPP
#define CIE_GL_TEXTURE_TEXTURE_IMPL_HPP

// --- External Includes ---
#include "glad/glad.h"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/traits/inc/types.hpp"
#include "packages/utility/inc/GLError.hpp"

namespace cie::gl {


template <class ValueType>
void Texture::write( const ValueType* p_data,
                     const std::vector<Size>& r_dimensions,
                     Size numberOfChannels )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Deduce format
    GLint format = 0;
    switch ( numberOfChannels )
    {
        case 1: format = GL_RED; break;
        case 2: format = GL_RG; break;
        case 3: format = GL_RGB; break;
        case 4: format = GL_RGBA; break;
        default:
            CIE_THROW( Exception, "Invalid number of channels (" + std::to_string(numberOfChannels) + ")" )
    }

    // Disable aligned reading (otherwise segfaults when reading exactly allocated data)
    // https://stackoverflow.com/questions/7380773/glteximage2d-segfault-related-to-width-height
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);

    // Load to GPU
    std::vector<GLint> dimensions(r_dimensions.size(), 0);

    switch( r_dimensions.size() )
    {
        case 1:
            glTexImage1D( this->_properties.target,
                          0,
                          this->_properties.format,
                          r_dimensions[0],
                          0,
                          format,
                          getGLTypeEnum<ValueType>(),
                          p_data );

            glGetTexLevelParameteriv( this->_properties.target,
                                      0,
                                      GL_TEXTURE_WIDTH,
                                      &dimensions[0] );
            CIE_CHECK(
                Size(dimensions[0]) == r_dimensions[0],
                "Loaded texture width mismatch"
            )
            break;
        case 2:
            glTexImage2D( this->_properties.target,
                          0,
                          this->_properties.format,
                          r_dimensions[0], r_dimensions[1],
                          0,
                          format,
                          getGLTypeEnum<ValueType>(),
                          p_data );
            glGetTexLevelParameteriv( this->_properties.target,
                                      0,
                                      GL_TEXTURE_WIDTH,
                                      &dimensions[0] );
            glGetTexLevelParameteriv( this->_properties.target,
                                      0,
                                      GL_TEXTURE_HEIGHT,
                                      &dimensions[1] );
            CIE_CHECK(
                Size(dimensions[0]) == r_dimensions[0],
                "Loaded texture width mismatch"
            )
            CIE_CHECK(
                Size(dimensions[1]) == r_dimensions[1],
                "Loaded texture height mismatch"
            )
            break;
        case 3:
            glTexImage3D( this->_properties.target,
                          0,
                          this->_properties.format,
                          r_dimensions[0], r_dimensions[1], r_dimensions[2],
                          0,
                          format,
                          getGLTypeEnum<ValueType>(),
                          p_data );
            glGetTexLevelParameteriv( this->_properties.target,
                                      0,
                                      GL_TEXTURE_WIDTH,
                                      &dimensions[0] );
            glGetTexLevelParameteriv( this->_properties.target,
                                      0,
                                      GL_TEXTURE_HEIGHT,
                                      &dimensions[1] );
            glGetTexLevelParameteriv( this->_properties.target,
                                      0,
                                      GL_TEXTURE_DEPTH,
                                      &dimensions[2] );
            CIE_CHECK(
                Size(dimensions[0]) == r_dimensions[0],
                "Loaded texture width mismatch"
            )
            CIE_CHECK(
                Size(dimensions[1]) == r_dimensions[1],
                "Loaded texture height mismatch"
            )
            CIE_CHECK(
                Size(dimensions[2]) == r_dimensions[2],
                "Loaded texture depth mismatch"
            )
            break;
        default:
            CIE_THROW( Exception, "Invalid number of texture dimensions (" + std::to_string(r_dimensions.size()) + ")" )
    }

    checkGLErrors( "Error writing data to texture" );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl

#endif