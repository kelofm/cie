// --- External Includes ---
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb/stb_image_resize.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/texture/inc/Image.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::gl {


Image::STBInitializer Image::_initializerDummy;


Image::Image() :
    _p_data(nullptr),
    _width(0),
    _height(0),
    _numberOfChannels(0)
{
}


Image::Image(const std::filesystem::path& r_filePath) :
    Image()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->load(r_filePath);

    CIE_END_EXCEPTION_TRACING
}


Image::Image(Size width, Size height, Size numberOfChannels) :
    _p_data(nullptr),
    _width(width),
    _height(height),
    _numberOfChannels(numberOfChannels)
{
    this->_p_data = this->allocate(width, height, numberOfChannels);
}


Image::Image(Image&& r_rhs)
    : Image()
{
    CIE_BEGIN_EXCEPTION_TRACING

    *this = r_rhs;

    CIE_END_EXCEPTION_TRACING
}


Image::Image(const Image& r_rhs) :
    Image()
{
    CIE_BEGIN_EXCEPTION_TRACING

    *this = r_rhs;

    CIE_END_EXCEPTION_TRACING
}


Image& Image::operator=(Image&& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _p_data = r_rhs._p_data;
    _width = r_rhs._width;
    _height = r_rhs._height;
    _numberOfChannels = r_rhs._numberOfChannels;

    r_rhs._p_data = nullptr;
    r_rhs._width = 0;
    r_rhs._height = 0;
    r_rhs._numberOfChannels = 0;

    return *this;

    CIE_END_EXCEPTION_TRACING
}


Image& Image::operator=(const Image& r_rhs)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_p_data = this->allocate(r_rhs.width(),
                                    r_rhs.height(),
                                    r_rhs.numberOfChannels());

    this->_width            = r_rhs.width();
    this->_height           = r_rhs.height();
    this->_numberOfChannels = r_rhs.numberOfChannels();

    std::copy(r_rhs._p_data,
               r_rhs._p_data + r_rhs.size(),
               this->_p_data);

    return *this;

    CIE_END_EXCEPTION_TRACING
}


Image::~Image()
{
    this->deallocate(this->_p_data);
}


void Image::load(const std::filesystem::path& r_filePath)
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->deallocate(this->_p_data);

    auto extension = r_filePath.extension().string();
    if (extension != ".jpeg"
         && extension != ".jpg"
         && extension != ".png"
         && extension != ".tga"
         && extension != ".bmp"
         && extension != ".psd"
         && extension != ".gif"
         && extension != ".hdr"
         && extension != ".pic"
         && extension != ".pnm")
        CIE_THROW(Exception, "Unsupported input image format: " + extension)

    int width=0, height=0, numberOfChannels=0;

    this->_p_data = stbi_load(
        r_filePath.string().c_str(),
        &width, &height, &numberOfChannels,
        0
   );

    if (!this->_p_data)
        CIE_THROW(Exception, "Failed to load " + r_filePath.string())

    //if (flip)
    //    stbi_set_flip_vertically_on_load(0);

    this->_width            = Size(width);
    this->_height           = Size(height);
    this->_numberOfChannels = Size(numberOfChannels);

    CIE_END_EXCEPTION_TRACING
}


void Image::write(const std::filesystem::path& r_filePath, Size quality) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK(
        1 <= quality && quality <= 100,
        "Invalid output image quality (" + std::to_string(quality) + ")"
   )

    auto fileName       = r_filePath.string();
    auto extension      = r_filePath.extension().string();
    int operationResult = 0;

    if (extension == ".jpg")
        operationResult = stbi_write_jpg(
            fileName.c_str(),
            static_cast<int>(this->_width),
            static_cast<int>(this->_height),
            static_cast<int>(this->_numberOfChannels),
            (void*) this->_p_data,
            static_cast<int>(quality)
       );

    else if (extension == ".png")
        operationResult = stbi_write_png(
            fileName.c_str(),
            static_cast<int>(this->_width),
            static_cast<int>(this->_height),
            static_cast<int>(this->_numberOfChannels),
            (void*) this->_p_data,
            0
       );

    else if (extension == ".bmp")
        operationResult = stbi_write_bmp(
            fileName.c_str(),
            static_cast<int>(this->_width),
            static_cast<int>(this->_height),
            static_cast<int>(this->_numberOfChannels),
            (void*) this->_p_data
       );

    else if (extension == ".tga")
        operationResult = stbi_write_tga(
            fileName.c_str(),
            static_cast<int>(this->_width),
            static_cast<int>(this->_height),
            static_cast<int>(this->_numberOfChannels),
            (void*) this->_p_data
       );

    else
        CIE_THROW(Exception, "Unsupported output image format: " + fileName)

    if (!operationResult)
        CIE_THROW(Exception, "Error writing image " + fileName)

    CIE_END_EXCEPTION_TRACING
}


void Image::resize(Size width, Size height)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (width == this->_width && height == this->_height)
        return;

    if (width == 0 || height == 0)
        CIE_THROW(Exception, "Attempt to resize image to zero size (" + std::to_string(width) + "x" + std::to_string(height) + ")")

    Image::value_type* newData = this->allocate(width, height, this->_numberOfChannels);

    int operationResult = stbir_resize_uint8(
        this->_p_data,
        static_cast<int>(this->_width),
        static_cast<int>(this->_height),
        0,
        newData,
        static_cast<int>(width),
        static_cast<int>(height),
        0,
        static_cast<int>(this->_numberOfChannels)
   );

    if (!operationResult)
        CIE_THROW(
            Exception,
            "Failed to resize image from ("
            + std::to_string(this->_width) + "x" + std::to_string(this->_height) + ") to ("
            + std::to_string(width) + "x" + std::to_string(height) + ")"
       )

    this->deallocate(this->_p_data);

    this->_p_data = newData;
    this->_width  = width;
    this->_height = height;

    CIE_END_EXCEPTION_TRACING
}


void Image::verticalFlip()
{
    CIE_BEGIN_EXCEPTION_TRACING

    Size tripleWidth = 3 * this->_width;
    Size halfHeight = this->_height/2;

    for (Size rowIndex=0; rowIndex < halfHeight; ++rowIndex)
        std::swap_ranges(
            this->_p_data + tripleWidth * rowIndex,
            this->_p_data + tripleWidth * (rowIndex+1),
            this->_p_data + tripleWidth * (this->_height - rowIndex - 1)
       );

    CIE_END_EXCEPTION_TRACING
}


void Image::clear()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->_width = 0;
    this->_height = 0;
    this->_numberOfChannels = 0;

    this->deallocate(this->_p_data);

    CIE_END_EXCEPTION_TRACING
}


Image::value_type* Image::allocate(Size width, Size height, Size numberOfChannels) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    using value_type   = Image::value_type;
    Size numberOfBytes = width * height * numberOfChannels * sizeof(value_type);

    value_type* p_memory = (value_type*) malloc(numberOfBytes);

    if (!p_memory)
        CIE_THROW(MemoryException, "Failed to allocate memory of size " + std::to_string(numberOfBytes))

    return p_memory;

    CIE_END_EXCEPTION_TRACING
}


void Image::deallocate(Image::value_type*& rp_memory) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (rp_memory)
        free(rp_memory);

    rp_memory = nullptr;

    CIE_END_EXCEPTION_TRACING
}


Image::STBInitializer::STBInitializer()
{
    stbi_set_flip_vertically_on_load(1);
}


} // namespace cie::gl