#ifndef CIE_GL_VERTEX_ARRAY_OBJECT_HPP
#define CIE_GL_VERTEX_ARRAY_OBJECT_HPP

// --- External Includes ---
#include "glad/glad.h"

// --- Utility Includes ---
#include "packages/types/inc/IDObject.hpp"
#include "packages/logging/inc/Loggee.hpp"


namespace cie::gl {


class VertexArrayObject : public utils::IDObject<GLuint>, public utils::Loggee
{
public:
    VertexArrayObject();

    VertexArrayObject(utils::Logger& r_logger, const std::string& r_name = "VertexArrayObject");

    VertexArrayObject(VertexArrayObject&& r_rhs) = default;

    VertexArrayObject(const VertexArrayObject& r_rhs) = delete;

    VertexArrayObject& operator=(VertexArrayObject&& r_rhs) = delete;

    VertexArrayObject& operator=(const VertexArrayObject& r_rhs) = delete;

    ~VertexArrayObject();
}; // class VertexArrayObject


} // namespace cie::gl


#endif