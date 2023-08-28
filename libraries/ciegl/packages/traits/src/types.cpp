// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/traits/inc/types.hpp"

// --- STL Includes ---
#include <map>


namespace cie::gl {


const std::map<std::string,GLenum> typeMap {

    // Integral types
    { "char", GL_BYTE },
    { "unsigned char", GL_UNSIGNED_BYTE },
    { "short", GL_SHORT },
    { "unsigned short", GL_UNSIGNED_SHORT },
    { "int", GL_INT },
    { "unsigned int", GL_UNSIGNED_INT },
    { "half", GL_HALF_FLOAT },
    { "float", GL_FLOAT },
    { "double", GL_DOUBLE },

    // Vector types (incomplete)
    { "vec2", GL_FLOAT_VEC2 },
    { "vec3", GL_FLOAT_VEC3 },
    { "vec4", GL_FLOAT_VEC4 },

    // Matrix types (incomplete)
    { "mat2", GL_FLOAT_MAT2 },
    { "mat3", GL_FLOAT_MAT3 },
    { "mat4", GL_FLOAT_MAT4 }
};


GLenum getGLTypeEnum( const std::string& r_typeName )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto it = typeMap.find( r_typeName );

    if ( it == typeMap.end() )
        CIE_THROW( Exception, "Unrecognized type name: " + r_typeName )

    return it->second;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl