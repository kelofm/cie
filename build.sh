#!/bin/bash

scriptName="$(basename ${BASH_SOURCE[0]})"
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
projectName="$(basename $scriptDir)"
projectNameUpper=$(printf '%s\n' "$projectName" | awk '{ print toupper($0) }')

print_help() {
    projectNames=""
    for item in $(ls "$scriptDir/libraries"); do
        if [ -d "$scriptDir/libraries/$item" ]; then
            if [ "$projectNames" = "" ]; then
                projectNames=$item
            else
                projectNames="$projectNames, $item"
            fi
        fi
    done
    for item in $(ls "$scriptDir/executables"); do
        if [ -d "$scriptDir/executables/$item" ]; then
            if [ "$projectNames" = "" ]; then
                projectNames=$item
            else
                projectNames="$projectNames, $item"
            fi
        fi
    done

    echo "$scriptName - Configure, build, and install $projectName."
    echo "Usage: $scriptName [OPTION [ARGUMENT]]"
    echo "-a project-name   : Add a project to the list of compiled ones. Options are [$projectNames]."
    echo "-b build-dir      : Build directory."
    echo "-c compiler-name  : Compiler family [gcc, clang, intel, acpp] (Default: gcc)."
    echo "-h                : Print this help and exit."
    echo "-i install-dir    : Install directory (python site package by default)."
    echo "-j job-count      : Number of jobs to launch the compilation with (default is the number of physical cores)."
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
installDir="$HOME/.local"
generator="Unix Makefiles"
cCacheFlag=""
cxx="g++"
projects=()                             # <== list of projects to build
cmakeArguments=""                       # <== semicolon-separated list of options to pass to CMake
cmakeCxxFlags=""                        # <== value to set for CMAKE_CXX_FLAGS
cmakeLinkerFlags=""                     # <== value to set for CMAKE_{EXE,SHARED,MODULE}_LINKER_FLAGS
openmpRoot=""                           # <== root of the OpenMP runtime (acpp only)
cmakeToolchainFile=""                   # <== toolchain file to use (macOS + system compilers only)
enableSYCL="OFF"                        # <== set to true if the requested compiler is acpp
jobCount=""                             # <== number of processes to launch the compilation with

case "$(uname -s)" in
    Linux*)
        jobCount=$(grep "^cpu\\scores" /proc/cpuinfo | uniq |  awk '{print $4}')
        ;;
    Darwin*)
        jobCount=$(sysctl -n machdep.cpu.thread_count)
        # The system compiler on macOS cannot locate the C++ standard library,
        # so the Homebrew LLVM toolchain file is required for it.
        cmakeToolchainFile="-DCMAKE_TOOLCHAIN_FILE:STRING=${scriptDir}/cmake/homebrew_llvm_toolchain.cmake"
        ;;
    \?)
        echo "Error: unsupported OS $(uname -s)"
        exit 1
esac

# Parse command line arguments.
while getopts ":a: h p t: b: c: i: j: o:" arg; do
    case "$arg" in
        a)  # Add a project to the list of compiled ones.
            [[ "$OPTARG" = "cieutils" || "$OPTARG" = "linalg" || "$OPTARG" = "geo" || "$OPTARG" = "fem" || "$OPTARG" = "ciegl" || "$OPTARG" = "bad_apple" || "$OPTARG" = "benchmarks" || "$OPTARG" = "fem_examples" || "$OPTARG" = "geo_examples" ]] || (print_help && echo "Invalid project name: '$OPTARG'" && exit 1)
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
                if [ "$cmakeCxxFlags" = "" ]; then
                    cmakeCxxFlags="-fsycl"
                else
                    cmakeCxxFlags="${cmakeCxxFlags} -fsycl "
                fi
                enableSYCL="ON"
            elif [ "$compilerFamily" = "acpp" ]; then
                export cc="$(which acpp)"
                export cxx="$(which acpp)"
                enableSYCL="ON"
                # The toolchain file would replace the acpp driver with clang++.
                cmakeToolchainFile=""
                # The acpp driver links against the OpenMP runtime of the LLVM it
                # was built with. That runtime is not necessarily on the default
                # linker search path, so locate it and tell both the linker and
                # CMake where to find it.
                openmpRoot="$("$scriptDir/cmake/find_openmp_runtime.sh" "$cxx")"
                if [ $? -ne 0 ]; then
                    echo "Warning: could not locate the OpenMP runtime for '$cxx'." >&2
                    echo "         Linking will likely fail with \"library 'omp' not found\"." >&2
                    echo "         Install an OpenMP runtime (e.g. 'brew install libomp'," >&2
                    echo "         'apt install libomp-<version>-dev') or point the compiler" >&2
                    echo "         at an existing installation." >&2
                fi
                if [ "$openmpRoot" != "" ]; then
                    cmakeLinkerFlags="${cmakeLinkerFlags}-L${openmpRoot}/lib -Wl,-rpath,${openmpRoot}/lib "
                fi
            else
                echo "Error: unsupported compiler family: $compilerFamily"
                exit 1
            fi
            ;;
        i)  # Set install directory.
            installDir="$OPTARG"
            ;;
        j)  # Set the number of compiler instances.
            jobCount="$OPTARG"
            ;;
        o)  # Append CMake arguments.
            if [[ "$OPTARG" == -DCMAKE_CXX_FLAGS* ]]; then
                cmakeCxxFlags="${cmakeCxxFlags}${OPTARG#*-DCMAKE_CXX_FLAGS=} "
            elif [[ "$OPTARG" =~ ^-DCMAKE_(EXE|SHARED|MODULE)_LINKER_FLAGS=(.*)$ ]]; then
                cmakeLinkerFlags="${cmakeLinkerFlags}${BASH_REMATCH[2]} "
            else
                if [ "$cmakeArguments" = "" ]; then
                    cmakeArguments="$OPTARG"
                else
                    cmakeArguments="$cmakeArguments;$OPTARG"
                fi
            fi
            ;;
        \?) # Unrecognized argument.
            print_help
            echo "Error: unrecognized argument: -$OPTARG"
            exit 1
            ;;
    esac
done

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

# Collect compiler family specific CMake arguments
extraCMakeArguments=()
if [ "$cmakeLinkerFlags" != "" ]; then
    extraCMakeArguments+=(
        "-DCMAKE_EXE_LINKER_FLAGS:STRING=${cmakeLinkerFlags}"
        "-DCMAKE_SHARED_LINKER_FLAGS:STRING=${cmakeLinkerFlags}"
        "-DCMAKE_MODULE_LINKER_FLAGS:STRING=${cmakeLinkerFlags}")
fi
if [ "$openmpRoot" != "" ]; then
    extraCMakeArguments+=("-DOpenMP_ROOT:STRING=${openmpRoot}")
fi
if [ "$cmakeToolchainFile" != "" ]; then
    extraCMakeArguments+=("$cmakeToolchainFile")
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
    "-DCIE_ENABLE_HDF5=ON"                                  \
    "$cCacheFlag"                                           \
    "${extraCMakeArguments[@]}"                             \
    "${cmakeArguments[@]}"                                  \
    "${projects[@]}"                                        \
    ; then
    exit 1
fi

# Build and install
if ! cmake --build "$buildDir" --config "$buildType" --target install -j $jobCount; then
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
