#ifndef CIEUTILS_CMAKE_VARIABLES_HPP
#define CIEUTILS_CMAKE_VARIABLES_HPP

// --- STL Includes ---
#include <filesystem>
#include <string>


namespace cie {


std::filesystem::path getInstallPath();

std::filesystem::path getBinaryPath();

std::filesystem::path getSourcePath();

std::filesystem::path getDataPath();

std::filesystem::path getOutputPath();

bool isMPIEnabled();

std::string getCompiler();

std::string getCompilerVersion();

std::string getCompileFlags();

std::string getCompilePlatform();


} // namespace cie

#endif