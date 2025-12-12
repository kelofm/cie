# Table of Contents

- [Building from Source](#building-from-source)
    - [CMake](#building-from-source-cmake)
    - [Shell Script](#building-from-source-shell-script)
    - [SYCL Support](#building-from-source-sycl-support)
        - [Arch](#building-from-source-sycl-support-arch)
    - [Sanitizers](#building-from-source-sanitizers)
        - [GCC](#building-from-source-sanitizers-gcc)
        - [Clang](#building-from-source-sanitizers-clang)

<a name="building-from-source"></a>

# Building from Source

<a name="building-from-source-cmake"></a>

## CMake

`CiE` consists of several components that include libraries (`/libraries`) and executables (`/executables`). A root `CMake` configuration is available, but each component can be configured separately as well, if necessary.

<a name="building-from-source-shell-script"></a>

## Shell Script

Alternatively, a build script (`/build.sh`) is provided for a bit more convenience.
```
$> ./build.sh -h

build.sh - Configure, build, and install cie.
Usage: build.sh [OPTION [ARGUMENT]]
-a project-name   : Add a project to the list of compiled ones. Options are [ciegl, cieutils, fem, geo, linalg, bad_apple, benchmarks].
-b build-dir      : Build directory.
-c compiler-name  : Compiler family [gcc, clang, intel, acpp] (Default: gcc).
-h                : Print this help and exit.
-i install-dir    : Install directory (python site package by default).
-o [opts]         : Options/arguments to pass on to CMake. Semicolon (;) delimited, or defined repeatedly.
-p                : Package after building.
-t build-type     : Build type [Debug, Release, RelWithDebInfo] (default: Release).
```

<a name="building-from-source-sycl-support"></a>

## SYCL Support

At the moment, `SYCL` features are only supported when compiling with `AdaptiveCpp`, which can be set by passing `-c acpp` to the build script.

```bash
./build.sh -a cieutils -c acpp
```

<a name="building-from-source-sycl-support-arch"></a>

### Arch

The `AdaptiveCpp` package is only available in the *AUR* and is highly unstable. It might be counter-intuitive but building it from source is the more reliable alternative, though it may need a specific `LLVM` version (e.g.: `llvm-19` at the moment, available in the *AUR*).

<a name="building-from-source-sanitizers"></a>

## Sanitizers

Building with sanitizers involves manipulating the compiler flags at `CMake` configuration time, as well as setting `LD_PRELOAD` before running built executables.

<a name="building-from-source-sanitizers-gcc"></a>

### GCC

Configuration
```bash
./build.sh                                  \
    -c gcc                                  \
    -a cieutils                             \
    -o -DCMAKE_CXX_FLAGS=-fsanitize=address \
    -o -DCMAKE_CXX_FLAGS=-shared-libasan
```

Runtime environment
```bash
export LD_PRELOAD=$(gcc -print-file-name=libasan.so)
```

<a name="building-from-source-sanitizers-clang"></a>

### Clang

*Note: compiling with `LLVM` sanitizers requires using `clang`'s standard library, and requires all `C++` dependencies be built with `clang`'s standard libraries.*

Configuration
```bash
./build.sh                                  \
    -c gcc                                  \
    -a cieutils                             \
    -o -DCMAKE_CXX_FLAGS=-fsanitize=address \
    -o -DCMAKE_CXX_FLAGS=-shared-libasan    \
    -o -DCMAKE_CXX_FLAGS=-stdlib=libc++     \
    -o -DCMAKE_EXE_LINKER_FLAGS=-lc++
```

Runtime environment
```bash
export LD_PRELOAD=$(clang -print-file-name=libclang_rt.asan-x86_64.so)
```
