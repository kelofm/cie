#ifndef CIE_GL_IMAGE_IMPL_HPP
#define CIE_GL_IMAGE_IMPL_HPP

// --- Graphics Includes ---
#include "packages/texture/inc/Image.hpp"


namespace cie::gl {


inline Size Image::width() const noexcept
{
    return this->_width;
}


inline Size Image::height() const noexcept
{
    return this->_height;
}


inline Size Image::numberOfChannels() const noexcept
{
    return this->_numberOfChannels;
}


inline Size Image::size() const noexcept
{
    return this->_width * this->_height * this->_numberOfChannels;
}


inline const Image::value_type* Image::data() const noexcept
{
    return this->_p_data;
}


inline Image::value_type* Image::data() noexcept
{
    return this->_p_data;
}


} // namespace cie::gl


#endif
