#!/bin/bash

scriptName="$(basename ${BASH_SOURCE[0]})"
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
projectName="$(basename $scriptDir)"
projectNameUpper=$(printf '%s\n' "$projectName" | awk '{ print toupper($0) }')

print_help() {
    echo "$scriptName - Configure, build, and install $projectName"
    echo "Usage: $scriptName [OPTION [ARGUMENT]]"
    echo "-a                : Add project to the list of compiled ones. Options are [cieutils, linalg, geo, fem, ciegl, bad_apple, benchmarks]."
    echo "-b build-dir      : Build directory."
    echo "-c compiler-name  : Compiler family [gcc, clang, intel, acpp] (Default: gcc)."
    echo "-h                : Print this help and exit."
    echo "-i install-dir    : Install directory (python site package by default)."
    echo "-o [opts]         : Options/arguments to pass on to CMake. Semicolon (;) delimited, or defined repeatedly."
    echo "-p                : Package after building."
    echo "-t build-type     : Build type [Debug, Release, RelWithDebInfo] (default: Release)."
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
cxx="g++"
projects=()                             # <== list of projects to build
cmakeArguments=""                       # <== semicolon-separated list of options to pass to CMake
cmakeCxxFlags=""                        # <== value to set for CMAKE_CXX_FLAGS
enableSYCL="OFF"                        # <== set to true if the requested compiler is acpp

while getopts ":a: h p t: b: c: i: o:" arg; do
    case "$arg" in
        a)  # Add a project to the list of compiled ones.
            [[ "$OPTARG" = "cieutils" || "$OPTARG" = "linalg" || "$OPTARG" = "geo" || "$OPTARG" = "fem" || "$OPTARG" = "ciegl" || "$OPTARG" = "bad_apple" || "$OPTARG" = "benchmarks" ]] || (print_help && echo "Invalid project name: '$buildType'" && exit 1)
            projectNameUpper=$(printf '%s\n' "$OPTARG" | awk '{ print toupper($0) }')
            projects+=("-DCIE_ENABLE_$projectNameUpper=ON")
            ;;
        h)  # Print help and exit without doing anything.
            print_help
            exit 0
            ;;
        p)  # Package after building.
            package=1
            ;;
        t)  # Set build type.
            buildType="$OPTARG"
            [[ "${buildType}" = "Debug" || "${buildType}" = "RelWithDebInfo" || "${buildType}" = "Release" ]] || (print_help && echo "Invalid build type: $buildType" && exit 1)
            ;;
        b)  # Set build directory.
            buildDir="$OPTARG"
            ;;
        c)  # Select compilers.
            compilerFamily="${OPTARG}"
            if [ "$compilerFamily" = "gcc" ]; then
                export cc="$(which gcc)"
                export cxx="$(which g++)"
            elif [ "$compilerFamily" = "clang" ]; then
                export cc="$(which clang)"
                export cxx="$(which clang++)"
            elif [ "$compilerFamily" = "intel" ]; then
                if [ -f "/opt/intel/oneapi/setvars.sh" ] ; then
                    source "/opt/intel/oneapi/setvars.sh" intel64
                fi
                export cc="$(which icx)"
                export cxx="$(which icpx)"
            elif [ "$compilerFamily" = "acpp" ]; then
                export cc="$(which acpp)"
                export cxx="$(which acpp)"
                enableSYCL="ON"
            else
                echo "Error: unsupported compiler family: $compilerFamily"
                exit 1
            fi
            ;;
        i)  # Set install directory.
            installDir="$OPTARG"
            ;;
        o)  # Append CMake arguments.
            if [[ "$OPTARG" == -DCMAKE_CXX_FLAGS* ]]; then
                cmakeCxxFlags="${cmakeCxxFlags}${OPTARG#*-DCMAKE_CXX_FLAGS=} "
            else
                cmakeArguments="$cmakeArguments;$OPTARG"
            fi
            ;;
        \?) # Unrecognized argument.
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
    "-DCIE_ENABLE_SYCL=${enableSYCL}"                       \
    "$cCacheFlag"                                           \
    "${cmakeArguments[@]}"                                  \
    "${projects[@]}"                                        \
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
