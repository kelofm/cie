#ifndef CIE_GL_SHADER_IMPL_HPP
#define CIE_GL_SHADER_IMPL_HPP

// --- GL Includes ---
#include "packages/shaders/inc/Shader.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::gl {


template <class TOStream>
TOStream& Shader::print(TOStream& r_stream) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    r_stream << "Shader loaded from '" << this->sourcePath()
    << "' and configured from '" << this->configurationPath() << "'";

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl


#endif