// --- Internal Includes ---
#include "packages/scene//inc/AbsProgram.hpp"

// --- GL Includes ---
#include "packages/utility/inc/GLError.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/types/inc/IDObject.hpp"

// --- STL Includes ---
#include <limits>
#include <optional>
#include <mutex>


namespace cie::gl {


struct AbsProgram::PrivateState
{
    PrivateState()
        : isLinked(false),
          isActive(false),
          _mutex()
    {}

    ~PrivateState() = default;

    bool isLinked;

    bool isActive;

    void unsetActiveProgramID()
    {
        std::scoped_lock<std::mutex> lock(_mutex);
        _activeProgramID.reset();
    }

    void setActiveProgramID(Size id)
    {
        std::scoped_lock<std::mutex> lock(_mutex);
        _activeProgramID = id;
    }

    std::optional<Size> getActiveProgramID()
    {
        std::scoped_lock<std::mutex> lock(_mutex);
        return _activeProgramID;
    }

private:
    std::mutex _mutex;

    static std::optional<Size> _activeProgramID;
}; // struct AbsProgram::PrivateState


std::optional<Size> AbsProgram::PrivateState::_activeProgramID = std::optional<Size>();


AbsProgram::ActiveContextToken::ActiveContextToken(AbsProgram& r_program, bool log)
    : _r_program(r_program),
      _log(log)
{
}


AbsProgram::ActiveContextToken::~ActiveContextToken()
{
    _r_program.deactivate(_log);
}


AbsProgram::ActiveContextToken AbsProgram::activate(bool log)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (!this->isActive())
    {
        if (log)
            this->logID("Set as active program", this->getID());
        glUseProgram(this->getID());

        if (log)
            this->logID("Bind Vertex Array Object", this->getID());

        _p_state->setActiveProgramID(this->getID());

        checkGLErrors(*this, "Activation failed!\n" + this->getErrorString());
    }

    return ActiveContextToken(*this, log);

    CIE_END_EXCEPTION_TRACING
}


bool AbsProgram::isActive()
{
    CIE_BEGIN_EXCEPTION_TRACING

    return _p_state->getActiveProgramID() == this->getID();

    CIE_END_EXCEPTION_TRACING
}


void AbsProgram::deactivate(bool log)
{
    CIE_BEGIN_EXCEPTION_TRACING

    if (this->isActive())
    {
        if (log)
            this->logID("Deactivate", this->getID());
        _p_state->unsetActiveProgramID();
    }

    CIE_END_EXCEPTION_TRACING
}


AbsProgram::AbsProgram() : AbsProgram(utils::LoggerSingleton::get())
{
}


AbsProgram::AbsProgram(utils::Logger& r_logger, const std::string& r_name)
    : utils::Loggee(r_logger, r_name),
      utils::IDObject<Size>(std::numeric_limits<Size>().max()),
      _p_bufferManager(new BufferManager(r_logger)),
      _p_state()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->setID(glCreateProgram());
    this->logID("Created a new program", this->getID());

    CIE_END_EXCEPTION_TRACING
}


AbsProgram::~AbsProgram()
{
    this->logID("Delete", this->getID());
    glDeleteProgram(this->getID());

    //this->logID("Delete Vertex Array Object", this->_vao->getID());
}


void AbsProgram::link()
{
    CIE_BEGIN_EXCEPTION_TRACING

    glLinkProgram(this->getID());

    GLint linkStatus = GL_FALSE;
    glGetProgramiv(this->getID(), GL_LINK_STATUS, &linkStatus);
    if (linkStatus != GL_TRUE)
        this->logID("Failed to link program!\n" + this->getErrorString(),
                    this->getID(),
                    LOG_TYPE_ERROR);

    CIE_END_EXCEPTION_TRACING
}


BufferManager& AbsProgram::getBufferManager()
{
    return *_p_bufferManager;
}


const BufferManager& AbsProgram::getBufferManager() const
{
    return *_p_bufferManager;
}


std::string AbsProgram::getErrorString() const
{
    CIE_BEGIN_EXCEPTION_TRACING

    constexpr Size bufferSize = 512;
    std::string buffer(bufferSize + 1, '\0');
    glGetProgramInfoLog(this->getID(),
                        bufferSize,
                        NULL,
                        &buffer[0]);

    buffer.resize(buffer.find('\0') - 1);
    return buffer;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl