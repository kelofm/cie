// --- Utility Includes ---
#include "packages/output/inc/FileManager.hpp"
#include "packages/concurrency/inc/MPISingleton.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <thread>
#include <string>
#include <sstream>


namespace cie::utils {


namespace {
inline static void appendPathName(Ref<FileManager::Path> r_path, RightRef<std::string> r_extra)
{
    if (!r_extra.empty())
    {
        const auto suffix = r_path.extension();
        r_path.replace_extension();
        r_path = (r_path.string() + std::move(r_extra));
        r_path.replace_extension(suffix);
    }
}


inline static std::string getThreadSuffix()
{
    std::stringstream stream;
    stream << "_thread_" << std::this_thread::get_id();
    return stream.str();
}


inline static std::string getMPISuffix()
{
    std::stringstream stream;

    if (mpi::MPISingleton::isDistributed())
        stream << "_rank_" << mpi::MPISingleton::get().getRankID();

    return stream.str();
}
} // unnamed namespace


void FileManager::makeThreadPrivate(Ref<Path> r_path)
{
    appendPathName(r_path, getThreadSuffix());
}


void FileManager::makeMPIPrivate(Ref<Path> r_path) const
{
    appendPathName(r_path, getMPISuffix());
}


void FileManager::makeConcurrentPrivate(Ref<Path> r_path) const
{
    std::string suffix = getThreadSuffix();
    suffix.append(getMPISuffix());
    appendPathName(r_path, std::move(suffix));
}


FileManager::Path FileManager::getOutputDirectory()
{
    return getOutputPath();
}


FileManager::Path FileManager::getCurrentDirectory()
{
    return std::filesystem::current_path();
}


} // namespace cie::utils
