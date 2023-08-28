// --- External Includes ---
#include "glad/glad.h"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/texture/inc/Texture.hpp"
#include "packages/utility/inc/GLError.hpp"


namespace cie::gl {


Texture::Texture( const Shader::Texture& r_shaderTexture ) :
    Shader::Texture( r_shaderTexture )
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Deduce texture type
    switch (this->_dimension)
    {
        case 1: this->_properties.target = GL_TEXTURE_1D; break;
        case 2: this->_properties.target = GL_TEXTURE_2D; break;
        case 3: this->_properties.target = GL_TEXTURE_3D; break;
        default:
            CIE_THROW( Exception, "Invalid texture dimension (" + std::to_string(this->_dimension) + ")" )
    }

    // Deduce texture format
    switch (this->_channels)
    {
        case 1: this->_properties.format = GL_RED; break;
        case 2: this->_properties.format = GL_RG; break;
        case 3: this->_properties.format = GL_RGB; break;
        case 4: this->_properties.format = GL_RGBA; break;
        default:
            CIE_THROW( Exception, "Invalid number of channels (" + std::to_string(this->_channels) + ")" )
    }

    // Generate GL texture and get its ID
    {
        GLuint id;
        glGenTextures( 1, &id );
        this->setID( id );
        checkGLErrors( "Failed to generate texture" );
    }

    glBindTexture( this->_properties.target, this->getID() );
    checkGLErrors( "Failed to bind texture | ID" + std::to_string(this->getID()) );

    glTexParameteri( this->_properties.target, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( this->_properties.target, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    checkGLErrors( "Failed to set texture properties" );

    CIE_END_EXCEPTION_TRACING
}


void Texture::write( const Image& r_image )
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->write( r_image.data(),
                 { r_image.width(), r_image.height() },
                 r_image.numberOfChannels() );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl