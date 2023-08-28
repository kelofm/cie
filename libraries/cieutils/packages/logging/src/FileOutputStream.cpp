// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/logging/inc/FileOutputStream.hpp"

// --- STL Inlcudes ---
#include <fstream>


namespace cie {


std::ofstream* makeNewFileStreamWithParentPath(const FileOutputStream::Path& r_path,
                                               FileOutputStream::OpenMode openMode)
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::filesystem::create_directories(r_path.parent_path());
    return new std::ofstream(r_path.string(), openMode);

    CIE_END_EXCEPTION_TRACING
}


FileOutputStream::FileOutputStream(const FileOutputStream::Path& r_filePath,
                                   FileOutputStream::OpenMode openMode)
    : OutputStream(makeNewFileStreamWithParentPath(r_filePath, openMode), true)
{
}


} // namespace cie