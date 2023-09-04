#ifndef CIE_GL_SHADER_DEFINITIONS_HPP
#define CIE_GL_SHADER_DEFINITIONS_HPP

// --- STL Includes ---
#include <string>
#include <optional>
#include <unordered_map>
#include <mutex>


namespace cie::gl {


class ShaderDefinitions
{
public:
    struct ShaderDefinition
    {
        std::string configuration;
        std::string source;
    };

public:
    static const ShaderDefinition& get(const std::string& r_key);

    static void clear();

private:
    static std::optional<std::unordered_map<
        std::string,
        ShaderDefinition
    >> _map;

    static std::mutex _mutex;
};


} // namespace cie::gl


#endif
