#ifndef CIE_GL_ABS_PROGRAM_HPP
#define CIE_GL_ABS_PROGRAM_HPP

// --- External Includes ---
#include "glad/glad.h"

// --- GL Includes ---
#include "packages/buffer/inc/BufferManager.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/Loggee.hpp"
#include "packages/types/inc/IDObject.hpp"


namespace cie::gl {


class AbsProgram : public utils::Loggee, public utils::IDObject<Size>
{
public:
    /// A class that disables the program that constructed it upon destruction
    class ActiveContextToken
    {
    public:
        friend class AbsProgram;

    public:
        ~ActiveContextToken();

    protected:
        ActiveContextToken(AbsProgram& r_program, bool log);

    private:
        ActiveContextToken(ActiveContextToken&& r_rhs) = delete;

        ActiveContextToken(const ActiveContextToken& r_rhs) = delete;

        ActiveContextToken& operator=(ActiveContextToken&& r_rhs) = delete;

        ActiveContextToken& operator=(const ActiveContextToken& r_rhs) = delete;

    private:
        AbsProgram& _r_program;

        bool _log;
    }; // class ActiveContextToken

public:
    AbsProgram(AbsProgram&& r_rhs) = default;

    virtual ~AbsProgram();

    virtual void update() = 0;

    [[nodiscard]] ActiveContextToken activate(bool log = false);

    bool isActive();

protected:
    AbsProgram();

    AbsProgram(utils::Logger& r_logger, const std::string& r_name = "Program");

    void link();

    /// This doesn't do anything except log a request to disable the program
    virtual void deactivate(bool log = false);

    BufferManager& getBufferManager();

    const BufferManager& getBufferManager() const;

    std::string getErrorString() const;

private:
    AbsProgram(const AbsProgram& r_rhs) = delete;

    AbsProgram& operator=(AbsProgram&& r_rhs) = delete;

    AbsProgram& operator=(const AbsProgram& r_rhs) = delete;

private:
    BufferManager::UniquePointer _p_bufferManager;

    struct PrivateState;
    std::unique_ptr<PrivateState> _p_state;
}; // class AbsProgram


} // namespace cie::gl


#endif