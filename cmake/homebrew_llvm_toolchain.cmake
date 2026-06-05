# Find dependencies installed from Homebrew.
execute_process(COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/cmake/find_homebrew_dependency.sh" llvm OUTPUT_VARIABLE llvm_directory)
execute_process(COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/cmake/find_homebrew_dependency.sh" libomp OUTPUT_VARIABLE OpenMP_ROOT)

# Find dependencies compiled locally.
#set(nlohmann_json_schema_validator_DIR "$ENV{HOME}/.local/lib")

# Select the compiler and point env vars to it.
set(CMAKE_C_COMPILER "${llvm_directory}/bin/clang")
set(CMAKE_CXX_COMPILER "${llvm_directory}/bin/clang++")

# Set MacOS-related crap.
set(CMAKE_CXX_FLAGS_INIT "-stdlib=libc++ -I${llvm_directory}/include/c++/v1")
set(CMAKE_EXE_LINKER_FLAGS_INIT "-L${llvm_directory}/lib -Wl,-rpath,${llvm_directory}/lib -L${OpenMP_ROOT}/lib")
set(CMAKE_STATIC_LINKER_FLAGS_INIT "-L${llvm_directory}/lib -Wl,-rpath,${llvm_directory}/lib -L${OpenMP_ROOT}/lib")
set(CMAKE_SHARED_LINKER_FLAGS_INIT "-L${llvm_directory}/lib -Wl,-rpath,${llvm_directory}/lib -L${OpenMP_ROOT}/lib")
