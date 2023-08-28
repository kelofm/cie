#ifndef CIE_UTILS_FILE_OUTPUT_STREAM_HPP
#define CIE_UTILS_FILE_OUTPUT_STREAM_HPP

// --- Internal Includes ---
#include "packages/logging/inc/OutputStream.hpp"

// --- STL Includes ---
#include <filesystem>


namespace cie {


class FileOutputStream : public OutputStream
{
public:
    using Path = std::filesystem::path;
    using OpenMode = std::ios_base::openmode;

public:
    FileOutputStream(const Path& r_filePath,
                     OpenMode openMode = std::ios::out);

    FileOutputStream(FileOutputStream&& r_rhs) = default;

    FileOutputStream(const FileOutputStream& r_rhs) = delete;

    FileOutputStream& operator=(FileOutputStream&& r_rhs) = default;

    FileOutputStream& operator=(const FileOutputStream& r_rhs) = delete;
};


} // namespace cie

#endif