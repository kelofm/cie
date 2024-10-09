#ifndef CIE_GL_TEXTURE_IMAGE
#define CIE_GL_TEXTURE_IMAGE

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <filesystem>
#include <memory>


namespace cie::gl {


/// An interface class for the stb_image libraries
class Image
{
public:
    using value_type = unsigned char;

public:
    /// Initialize all to 0
    Image();

    /// Default constructor -> Image::load
    /**
     * @param r_filePath path to an image file
     */
    Image(const std::filesystem::path& r_filePath);

    /// Alloc constructor
    Image( Size width, Size height, Size numberOfChannels );

    Image( Image&& r_rhs ) ;

    /// Default constructor -> copy assignment operator
    Image( const Image& r_rhs );

    Image& operator=(Image&& r_rhs);

    /// Allocate and copy everything
    Image& operator=( const Image& r_rhs );

    ~Image();

    /// Load image from file
    /**
     * Supported image formats:
     * - jpeg (or jpg), png, tga, bmp, psd, gif, pic, pnm
     * see @ref stbi_load
     *
     * @param r_filePath path to an image file
     * @param flip flip y axis and move the origin to the bottom left
     * @note only the first frame will be loaded from gifs
     * @note extensions must be lowercase, otherwise they won't be recognised
     */
    void load( const std::filesystem::path& r_filePath);

    /// @brief Write image to file
    /** @details Method is deduced from the extension. Supported image formats:
     *            - png, bmp, tga, jpg
     *              (all lower case)
     *
     *  @param r_filePath path to file to be written to
     *  @param quality quality used for JPEG compression [1;100]
     *  @note other extensions (even uppercase) will throw an exception
     */
    void write( const std::filesystem::path& r_filePath, Size quailty=100 ) const;

    /// Resize image
    /**
     * @note no resizing is performed if the new width and height are identical
     * to the old ones
     * @note alpha channels are processed with default @ref stbi_resize behaviour
     */
    void resize( Size width, Size height );

    /// Flip the image vertically
    void verticalFlip();

    /// Delete and deallocate data, reset sizes
    void clear();

    /// Get image width in pixels
    Size width() const noexcept;

    /// Get image height in pixels
    Size height() const noexcept;

    /// Get number of channels in the image
    /**
     *  1: gray
     *  2: gray, alpha
     *  3: RGB
     *  4: RGBA
     */
    Size numberOfChannels() const noexcept;

    /// Get total number of components in the image
    Size size() const noexcept;

    const value_type* data() const noexcept;

    value_type* data() noexcept;

private:
    /// Attempt to allocate memory and check whether it was successful
    value_type* allocate( Size width, Size height, Size numberOfChannels ) const;

    /// Dellocate memory
    void deallocate( value_type*& rp_memory ) const;

    struct STBInitializer
    {
        STBInitializer();
    };

private:
    value_type* _p_data;

    Size        _width;

    Size        _height;

    Size        _numberOfChannels;

    // Sets the static flip flag in stb_image
    static STBInitializer _initializerDummy;
};


using ImagePtr = std::shared_ptr<Image>;


} // namespace cie::gl

#include "packages/texture/impl/Image_impl.hpp"

#endif