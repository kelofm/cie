cmake_minimum_required(VERSION 3.15.0)
project(benchmarks LANGUAGES CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)

# Declare all CMake settings
set(${PROJECT_NAME_UPPER}_ENABLE_PROFILING                  OFF CACHE BOOL "Enable profiling macros.")
set(${PROJECT_NAME_UPPER}_ENABLE_OUT_OF_RANGE_CHECKS        OFF CACHE BOOL "Enable macros for checking out-of-range errors.")
set(${PROJECT_NAME_UPPER}_ENABLE_DIVISION_BY_ZERO_CHECKS    OFF CACHE BOOL "Enable macros for checking division by zero errors.")
# Language server
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Collect dependencies
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}/cie/lib/cmake")
find_package(benchmark REQUIRED)
find_package(fem CONFIG REQUIRED)
find_package(geo CONFIG REQUIRED)
find_package(cieutils CONFIG REQUIRED)
find_package(linalg CONFIG REQUIRED)

# Build main executable
file(GLOB_RECURSE sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(${PROJECT_NAME} ${sources})

# Compile definitions
if(${PROJECT_NAME_UPPER}_ENABLE_PROFILING)
    add_compile_definitions(CIE_ENABLE_PROFILING=1)
endif()

if(${PROJECT_NAME_UPPER}_ENABLE_OUT_OF_RANGE_CHECKS)
    add_compile_definitions(CIE_ENABLE_OUT_OF_RANGE_CHECKS=1)
endif()

if(${PROJECT_NAME_UPPER}_ENABLE_DIVISION_BY_ZERO_CHECKS)
    add_compile_definitions(CIE_ENABLE_DIVISION_BY_ZERO_CHECKS=1)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-Wall -Wpedantic -Wextra -Wno-dangling-reference -Werror)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(-Wall -Wpedantic -Wextra -Werror)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    add_compile_options(-Wall -Wpedantic -Wextra -Werror)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(FATAL_ERROR "MSVC is not supported yet")
endif ()

target_link_libraries(${PROJECT_NAME} PRIVATE cie::fem benchmark::benchmark)
target_include_directories(${PROJECT_NAME} PRIVATE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")

set_target_properties(${PROJECT_NAME} PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
                      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
                      INSTALL_RPATH "$ORIGIN/../lib;$ORIGIN/../bin")
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION "cie/bin")

# Install scripts
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/scripts/"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/cie/scripts/"
        USE_SOURCE_PERMISSIONS)

# Package
include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_VENDOR "")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/readme")

set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /build
    *.gitignore)
include(CPack)

# Generate CMake config
include(CMakePackageConfigHelpers)
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in"
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION "bin/cmake/${PROJECT_NAME}"
                              NO_SET_AND_CHECK_MACRO
                              NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
                                 COMPATIBILITY AnyNewerVersion)

# Install targets
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        RUNTIME DESTINATION "cie/bin"
        INCLUDES DESTINATION "cie/include")
install(EXPORT ${PROJECT_NAME}Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE cie::
        DESTINATION "cie/bin/cmake/${PROJECT_NAME}")

# Install CMake config
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION "cie/bin/cmake/${PROJECT_NAME}")

# Export
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(TARGETS ${PROJECT_NAME}
       NAMESPACE cie::
       FILE ${PROJECT_NAME}Targets.cmake)
export(PACKAGE ${PROJECT_NAME})
