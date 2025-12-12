# Table of Contents

- [Building from Source](#building-from-source)
    - [SYCL Support](#building-from-source-sycl-support)
        - [Arch](#building-from-source-sycl-support-arch)
    - [Sanitizers](#building-from-source-sanitizers)
        - [GCC](#building-from-source-sanitizers-gcc)
        - [Clang](#building-from-source-sanitizers-clang)

<a name="building-from-source"></a>

# Building from Source

<a name="building-from-source-sycl-support"></a>

## SYCL Support

At the moment, `SYCL` features are only supported when compiling with `AdaptiveCpp`, which can be set by passing `-c acpp` to the build script.

```bash
build.sh -a cieutils -c acpp
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
build.sh                                    \
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
build.sh                                    \
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
