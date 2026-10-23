#!/usr/bin/env bash

# Locate the OpenMP runtime and print the prefix that contains it, i.e. the
# directory that holds the 'lib' subdirectory with libomp/libgomp in it.
# Nothing is printed when the runtime cannot be located.
#
# This is needed whenever the compiler cannot link '-fopenmp' on its own:
#   - On macOS Homebrew ships libomp as a formula separate from LLVM, so
#     '-lomp' is not on the default linker search path.
#   - On Linux the runtime usually lives in the LLVM installation, where the
#     compiler picks it up by itself.
#
# Exits with 1 if the compiler needs the runtime but it cannot be located.
#
# Usage: find_openmp_runtime.sh [compiler]

set -u

scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
compiler="${1:-}"

# The prefix of a runtime library is the parent of the directory holding it.
prefixOf() {
    [ -n "${1:-}" ] || return 1
    dirname "$(dirname "$1")"
}

probeSource="$(mktemp -t cie_openmp_probe).cpp"
probeBinary="${probeSource%.cpp}"
trap 'rm -f "$probeSource" "$probeBinary"' EXIT
printf 'int main() { return 0; }\n' > "$probeSource"

# If the compiler links an OpenMP program unaided, the runtime is already on
# the linker search path and there is nothing to do.
if [ -n "$compiler" ] && "$compiler" -fopenmp "$probeSource" -o "$probeBinary" >/dev/null 2>&1; then
    exit 0
fi

runtimeLibrary=""

# 1. Platform specific locations.
case "$(uname -s)" in
    Darwin)
        # Homebrew ships libomp as a formula separate from LLVM.
        runtimePrefix="$("$scriptDir/find_homebrew_dependency.sh" libomp 2>/dev/null)" || runtimePrefix=""
        if [ -n "$runtimePrefix" ]; then
            runtimeLibrary="$(ls "$runtimePrefix"/lib/libomp.* 2>/dev/null | head -1)"
        fi
        ;;
    Linux)
        # The LLVM installation the compiler belongs to.
        for llvmConfig in $(command -v llvm-config llvm-config-21 llvm-config-20 llvm-config-19 llvm-config-18 2>/dev/null); do
            libDir="$("$llvmConfig" --libdir 2>/dev/null)"
            if [ -n "$libDir" ] && ls "$libDir"/libomp.* >/dev/null 2>&1; then
                runtimeLibrary="$(ls "$libDir"/libomp.* | head -1)"
                break
            fi
        done
        # The dynamic linker cache.
        if [ -z "$runtimeLibrary" ]; then
            runtimeLibrary="$(ldconfig -p 2>/dev/null | grep -oE '/[^ ]*libomp\.so[^ ]*' | head -1)"
        fi
        ;;
esac

# 2. Ask the AdaptiveCpp driver which OpenMP runtime it was configured with.
if [ -z "$runtimeLibrary" ] && [ -n "$compiler" ]; then
    runtimeLibrary="$("$compiler" --acpp-dryrun "$probeSource" 2>/dev/null \
        | tr ' ' '\n' \
        | grep -E '/lib(omp|gomp|iomp5)[.a-zA-Z0-9]*$' \
        | head -1)"
fi

runtimePrefix="$(prefixOf "$runtimeLibrary")" || runtimePrefix=""
if [ -n "$runtimePrefix" ] && [ -d "$runtimePrefix/lib" ]; then
    printf '%s\n' "$runtimePrefix"
    exit 0
fi

# The compiler needs the runtime but it could not be located.
exit 1
