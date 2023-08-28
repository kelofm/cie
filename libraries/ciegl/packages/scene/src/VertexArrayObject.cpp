// --- Internal Includes ---
#include "packages/scene/inc/VertexArrayObject.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"

// --- STL Includes ---
#include <limits>

namespace cie::gl {


VertexArrayObject::VertexArrayObject(utils::Logger& r_logger, const std::string& r_name)
    : utils::IDObject<GLuint>(std::numeric_limits<GLuint>::max()),
      utils::Loggee(r_logger, r_name)
{
    GLuint id;
    glGenVertexArrays(1, &id);
    this->setID(id);
    this->logID("Create VAO", id);
}


VertexArrayObject::~VertexArrayObject()
{
    GLuint id = this->getID();
    this->logID("Destroy VAO", id);
    glDeleteVertexArrays(1, &id);
}


} // namespace cie::gl