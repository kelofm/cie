#ifndef CIE_GL_TEXTURE_TEXTURE_HPP
#define CIE_GL_TEXTURE_TEXTURE_HPP

// --- External Includes ---
#include "glad/glad.h"

// --- Utility Includes ---
#include "packages/types/inc/IDObject.hpp"

// --- Internal Includes ---
#include "packages/shaders/inc/Shader.hpp"
#include "packages/texture/inc/Image.hpp"

// --- STL Includes ---
#include <vector>
#include <memory>


namespace cie::gl {


/// Interface for loading textures onto the GPU
/**
 * Texture format is deduced from the number of channels:
 *  1:  GL_RED
 *  2:  GL_RG
 *  3:  GL_RGB
 *  4:  GL_RGBA
 */
class Texture :
    public Shader::Texture,
    public utils::IDObject<GLuint>
{
public:
    struct Properties
    {
        GLenum target;
        GLint  format;
        GLenum valueType;
    };

public:
    Texture( const Shader::Texture& r_shaderTexture );

    /// Load data from image into the texture on the GPU
    void write( const Image& r_image );

    /// Load data from an arbitrary* source into the texture on the GPU
    /**
     * @param p_data contiguous data to be copied
     * @param r_dimensions [width [,height] [,depth]]
     * @param numberOfChannels number of components per 'pixel'
     * @note r_dimensions must be of size 1-3
     * @note channels determine the data format:
     *  1:  GL_RED
     *  2:  GL_RG
     *  3:  GL_RGB
     *  4:  GL_RGBA
     */
    template <class ValueType>
    void write( const ValueType* p_data,
                const std::vector<Size>& r_dimensions,
                Size numberOfChannels );

    const Properties& properties() const;

private:
    Properties _properties;
};


using TexturePtr = std::shared_ptr<Texture>;


} // namespace cie::gl

#include "packages/texture/impl/Texture_impl.hpp"

#endif