#!/bin/bash

scriptName="$(basename ${BASH_SOURCE[0]})"
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
projectName="$(basename $scriptDir)"
projectNameUpper=$(printf '%s\n' "$projectName" | awk '{ print toupper($0) }')

print_help() {
    echo "$scriptName - Configure, build, and install $projectName"
    echo "Usage: $scriptName [OPTION [ARGUMENT]]"
    echo "-h                : print this help and exit"
    echo "-p                : package after building"
    echo "-t build-type     : build type [Debug, Release, RelWithDebInfo] (default: Release)"
    echo "-b build-dir      : build directory"
    echo "-i install-dir    : install directory (python site package by default)"
    echo "-o [opts]         : options/arguments to pass on to CMake. Semicolon (;) delimited, or defined repeatedly."
}

if ! command -v python3 &>/dev/null; then
    echo "Error: $scriptName requires python3, but could not find it on the system"
    exit 1
fi

# Function for getting the module paths associated with the current interpreter
get_site_packages_dir() {
    echo $(python3 -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
}

# Default arguments
package=0
buildType="Release"
buildDir="$scriptDir/build"
installDir=$(get_site_packages_dir)
generator="Unix Makefiles"
cCacheFlag=""
cmakeArguments=""                       # <== semicolon-separated list of options to pass to CMake
cmakeCxxFlags=""                        # <== value to set for CMAKE_CXX_FLAGS
cxx="g++"

while getopts ":h p t: b: i: o:" arg; do
    case "$arg" in
        h)  # Print help and exit without doing anything
            print_help
            exit 0
            ;;
        p)  # Package after building
            package=1
            ;;
        t)  # Set build type
            buildType="$OPTARG"
            [[ "${buildType}" = "Debug" || "${buildType}" = "RelWithDebInfo" || "${buildType}" = "Release" ]] || (print_help && echo "Invalid build type: $buildType" && exit 1)
            ;;
        b)  # Set build directory
            buildDir="$OPTARG"
            ;;
        i)  # Set install directory
            installDir="$OPTARG"
            ;;
        o)  # Append CMake arguments
            if [[ "$OPTARG" == -DCMAKE_CXX_FLAGS* ]]; then
                cmakeCxxFlags="${cmakeCxxFlags}${OPTARG#*-DCMAKE_CXX_FLAGS=} "
            else
                cmakeArguments="$cmakeArguments;$OPTARG"
            fi
            ;;
        \?) # Unrecognized argument
            print_help
            echo "Error: unrecognized argument: -$OPTARG"
            exit 1
            ;;
    esac
done

case "$(uname -s)" in
    Linux*)
        #export ACPP_TARGETS="hip:gfx1030"
        #toolchainRoot="/opt/hipSYCL/ROCm"
        #toolchainBin="$toolchainRoot/bin"
        #export cc="$toolchainBin/acpp"
        #export cxx="$toolchainBin/acpp"
        ;;
    Darwin*)
        # Set clang from homebrew
        if ! command -v brew &> /dev/null; then
            echo "Error: $scriptName requires Homebrew"
            exit 1
        fi

        foundPackage=""
        get_homebrew_package() {
            foundPackage=""
            package_versions="$(brew search --formula "/$1@[0-9]+/")"
            for package_version in $(echo $package_versions | tr ' ' '\n' | sort -r | tr '\n' ' '); do
                if brew list "$package_version" >/dev/null 2>&1; then
                    foundPackage="$package_version"
                    echo "using '$package_version' for dependency '$1'"
                    return 0
                fi
            done

            echo "Error: no installed version of '$1' was found."
            echo "Consider running 'brew install $1'."
            exit 1
        }

        get_homebrew_package llvm
        toolchainRoot="$(brew --prefix $foundPackage)"
        toolchainBin="${toolchainRoot}/bin"
        toolchainLib="${toolchainRoot}/lib"
        toolchainInclude="${toolchainRoot}/include"
        cxx="$toolchainBin/clang++"
        ;;
    \?)
        echo "Error: unsupported OS $(uname -s)"
        exit 1
esac

# Create or clear the build directory
if ! [ -d "$buildDir" ]; then
    mkdir -p "$buildDir"
else
    rm -f "$buildDir/cmake_install.cmake"
    rm -f "$buildDir/CMakeCache.txt"
    rm -rf "$buildDir/CMakeFiles"
fi

# Check optional dependency - ninja
if command -v ninja &>/dev/null; then
    generator="Ninja"
fi

# Check optional dependency - ccache
if command -v ccache &>/dev/null; then
    cCacheFlag="-DCMAKE_CXX_COMPILER_LAUNCHER:STRING=ccache"
fi

# Generate with CMake
if ! cmake                                                  \
    "-H$scriptDir"                                          \
    "-B$buildDir"                                           \
    "-DCMAKE_INSTALL_PREFIX:STRING=$installDir"             \
    "-G${generator}"                                        \
    "-DCMAKE_BUILD_TYPE:STRING=$buildType"                  \
    "-DCMAKE_CXX_COMPILER:STRING=$cxx"                      \
    "-DCMAKE_CXX_FLAGS=${cmakeCxxFlags}"                    \
    "-DCMAKE_COLOR_DIAGNOSTICS:BOOL=ON"                     \
    "-D${projectNameUpper}_BUILD_SHARED_LIBRARY:BOOL=ON"    \
    "-D${projectNameUpper}_BUILD_PYTHON_MODULE:BOOL=ON"     \
    "-D${projectNameUpper}_BUILD_TESTS:BOOL=ON"             \
    "$cCacheFlag"                                           \
    "${cmakeArguments[@]}"                                  \
    ; then
    exit 1
fi

physicalCores=$(grep "^cpu\\scores" /proc/cpuinfo | uniq |  awk '{print $4}')

# Build and install
if ! cmake --build "$buildDir" --config "$buildType" --target install -j $physicalCores; then
    exit 1
fi

# Package
if [ $package -eq 1 ]; then
    cd "$buildDir"
    if [ "$generator" = "Ninja" ]; then
        echo ninja package
        ninja package
        ninja package_source
    else
        echo make package
        make package
        make package_source
    fi
fi

# Success
exit 0
