#ifndef CIE_UTILS_FILE_MANAGER_HPP
#define CIE_UTILS_FILE_MANAGER_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <filesystem>


namespace cie::utils {


/// @addtogroup cieutils
class FileManager
{
public:
    using Path = std::filesystem::path;

public:
    /// @brief Add information about the current thread ID to the input path.
    static void makeThreadPrivate(Ref<Path> r_path);

    /// @brief Add information about the MPI rank to the input path.
    void makeMPIPrivate(Ref<Path> r_path) const;

    /// @brief Add information about the current thread ID and MPI rank to the input path.
    void makeConcurrentPrivate(Ref<Path> r_path) const;

    /// @brief Return the global CiE output directory.
    static Path getOutputDirectory();

    /// @brief Return the current working directory.
    static Path getCurrentDirectory();
}; // class FileManager


} // namespace cie::utils


#endif
