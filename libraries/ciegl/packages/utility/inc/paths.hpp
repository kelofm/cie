#ifndef CIE_GL_UTILITY_PATHS_HPP
#define CIE_GL_UTILITY_PATHS_HPP

// --- Internal Includes ---
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <filesystem>


namespace cie::gl {


// --- SOURCE PATHS --- //

const std::filesystem::path GL_DATA_SOURCE_PATH          = getSourcePath() / "libraries/ciegl/data";
const std::filesystem::path GL_CONFIGURATION_SOURCE_PATH = GL_DATA_SOURCE_PATH / "configuration";
const std::filesystem::path GL_GEOMETRY_SOURCE_PATH      = GL_DATA_SOURCE_PATH / "geometry";
const std::filesystem::path GL_SHADER_SOURCE_PATH        = GL_DATA_SOURCE_PATH / "shaders";
const std::filesystem::path GL_TEXTURE_SOURCE_PATH       = GL_DATA_SOURCE_PATH / "textures";


// --- INSTALL PATHS --- //

const std::filesystem::path GL_DATA_INSTALL_PATH          = getDataPath();
const std::filesystem::path GL_CONFIGURATION_INSTALL_PATH = GL_DATA_INSTALL_PATH / "configuration";
const std::filesystem::path GL_GEOMETRY_INSTALL_PATH      = GL_DATA_INSTALL_PATH / "geometry";
const std::filesystem::path GL_SHADER_INSTALL_PATH        = GL_DATA_INSTALL_PATH / "shaders";
const std::filesystem::path GL_TEXTURE_INSTALL_PATH       = GL_DATA_INSTALL_PATH / "textures";


} // namespace cie::gl

#endif
