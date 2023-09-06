cmake_minimum_required(VERSION 3.15.0)
project(ciegl LANGUAGES C CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)

# Declare all CMake settings
set(${PROJECT_NAME_UPPER}_BUILD_SHARED_LIBRARY              ON CACHE BOOL "build shared libraries instead of static ones")
set(${PROJECT_NAME_UPPER}_BUILD_PYTHON_MODULE               ON CACHE BOOL "build python bindings for ${PROJECT_NAME}")
set(${PROJECT_NAME_UPPER}_BUILD_TESTS                       ON CACHE BOOL "build testrunner for ${PROJECT_NAME}")
set(${PROJECT_NAME_UPPER}_ENABLE_MPI                        OFF CACHE BOOL "Build with MPI support")
set(${PROJECT_NAME_UPPER}_ENABLE_PROFILING                  OFF CACHE BOOL "Enable profiling macros.")
set(${PROJECT_NAME_UPPER}_ENABLE_OUT_OF_RANGE_CHECKS        OFF CACHE BOOL "Enable macros for checking out-of-range errors.")
set(${PROJECT_NAME_UPPER}_ENABLE_DIVISION_BY_ZERO_CHECKS    OFF CACHE BOOL "Enable macros for checking division by zero errors.")

# Language server
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Collect dependencies
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}/cie/lib")
find_package(cieutils CONFIG REQUIRED)
find_package(linalg CONFIG REQUIRED)
find_package(geo CONFIG REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)

# Build GLAD
add_library(glad STATIC "${CMAKE_CURRENT_SOURCE_DIR}/external/glad/src/glad.c")
target_include_directories(glad
                           PRIVATE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/glad/inc>"
                           INTERFACE "$<INSTALL_INTERFACE:cie/external/glad/inc>"
                                     "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/glad/inc>")
set_target_properties(glad PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_compile_definitions(glad PUBLIC GLAD_GLAPI_EXPORT PRIVATE GLAD_GLAPI_EXPORT_BUILD)

# Collect shader sources
set(FILL_MAP_IN_CMAKE "")
set(shader_types "vertexShader;geometryShader;fragmentShader")
file(GLOB shaderDirectory "${CMAKE_CURRENT_SOURCE_DIR}/data/shaders/*")
foreach(directory ${shaderDirectory})
    if(IS_DIRECTORY "${directory}")
        foreach(shader_type ${shader_types})
            if(EXISTS "${directory}/${shader_type}.glsl")
                if(EXISTS "${directory}/${shader_type}.json")
                    get_filename_component(shader_key "${directory}" NAME)
                    file(READ "${directory}/${shader_type}.glsl" shader_source)
                    file(READ "${directory}/${shader_type}.json" shader_configuration)
                    string(APPEND
                           FILL_MAP_IN_CMAKE
                           "    r_map.emplace(R\"(${shader_key}::${shader_type})\", ShaderDefinitions::ShaderDefinition {\nR\"(${shader_configuration})\"\n,\nR\"(\n${shader_source})\"}\n);\n")
                else()
                    message(FATAL_ERROR "Shader '${directory}/shader_type' is missing its configuration json")
                endif()
            endif()
        endforeach()
    endif()
endforeach()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/${PROJECT_NAME}_shaders.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/src/${PROJECT_NAME}_shaders.cpp")

# Withdraw write permission from the generated file
file(CHMOD
     "${CMAKE_CURRENT_SOURCE_DIR}/src/${PROJECT_NAME}_shaders.cpp"
     FILE_PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

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

add_compile_options(-Wall -Wpedantic -Wextra -Werror)

# Build main library
file(GLOB_RECURSE sources "${CMAKE_CURRENT_SOURCE_DIR}/packages/*/src/*.cpp")
if (${${PROJECT_NAME_UPPER}_BUILD_SHARED_LIBRARY})
    add_library(${PROJECT_NAME}_shaders SHARED "${CMAKE_CURRENT_SOURCE_DIR}/src/${PROJECT_NAME}_shaders.cpp")
    add_library(${PROJECT_NAME} SHARED ${sources})
    set_property(TARGET ${PROJECT_NAME}_shaders PROPERTY POSITION_INDEPENDENT_CODE 1)
    set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE 1)
else()
    add_library(${PROJECT_NAME} STATIC ${sources})
    add_library(${PROJECT_NAME}_shaders STATIC "${CMAKE_CURRENT_SOURCE_DIR}/src/${PROJECT_NAME}_shaders.cpp")
endif()

# GLM removed its legacy cmake target, but then readded it? idk
if(TARGET glm)
    set(libglm "glm")
elseif(TARGET glm::glm)
    set(libglm "glm::glm")
else()
    message(FATAL_ERROR "Neither glm nor glm::glm is a valid target")
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC cie::geo glfw OpenGL::GL ${libglm} glad ${PROJECT_NAME}_shaders)
target_include_directories(${PROJECT_NAME}
                           PRIVATE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                           INTERFACE "$<INSTALL_INTERFACE:cie/include>"
                                     "$<INSTALL_INTERFACE:cie/include/external/stb>"
                                     "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                                     "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/stb>")
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/stb>")
target_include_directories(${PROJECT_NAME}_shaders
                           PRIVATE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>"
                           INTERFACE "$<INSTALL_INTERFACE:cie/include/inc>"
                                     "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>")
foreach(target "${PROJECT_NAME};${PROJECT_NAME}_shaders")
    set_target_properties(${PROJECT_NAME} PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                          ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
                          INSTALL_RPATH "$ORIGIN")
endforeach()

# Build test executable
if (${PROJECT_NAME_UPPER}_BUILD_TESTS)
    file(GLOB_RECURSE test_sources "${CMAKE_CURRENT_SOURCE_DIR}/packages/*/test/*.cpp")
    list(INSERT test_sources 0 "${CMAKE_CURRENT_SOURCE_DIR}/test/${PROJECT_NAME}_testrunner.cpp")

    add_executable(${PROJECT_NAME}_testrunner ${test_sources})
    target_link_libraries(${PROJECT_NAME}_testrunner PRIVATE ${PROJECT_NAME})
    set_target_properties(${PROJECT_NAME}_testrunner PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

# Build python module
if (${PROJECT_NAME_UPPER}_BUILD_PYTHON_MODULE)
    # Each subdirectory in modules is assumed to be a python module
    # and is installed as-is
    file(GLOB children "${CMAKE_CURRENT_SOURCE_DIR}/python/modules/*")
    foreach(child ${children})
        if (IS_DIRECTORY "${child}")
            get_filename_component(module_name "${child}" NAME)
            install(DIRECTORY "${child}"
                    DESTINATION "cie")
        endif()
    endforeach()
endif()

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
                              INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
                              NO_SET_AND_CHECK_MACRO
                              NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
                                 COMPATIBILITY AnyNewerVersion)

# Install targets
install(TARGETS glad ${PROJECT_NAME}_shaders ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        LIBRARY DESTINATION "cie/lib"
        ARCHIVE DESTINATION "cie/lib"
        RUNTIME DESTINATION "cie/bin"
        INCLUDES DESTINATION "cie/include")
install(EXPORT ${PROJECT_NAME}Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE cie::
        DESTINATION "cie/lib/cmake/${PROJECT_NAME}")

# Install headers
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/packages"
        DESTINATION "cie/include"
        FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/inc"
        DESTINATION "cie/include"
        FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/external"
        DESTINATION "cie/include"
        FILES_MATCHING PATTERN "*.h*")

# Install data
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data/" DESTINATION "${CMAKE_INSTALL_PREFIX}/cie/data/")

# Install CMake config
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION "cie/lib/cmake/${PROJECT_NAME}")

# Export
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(TARGETS glad ${PROJECT_NAME}_shaders ${PROJECT_NAME}
       NAMESPACE cie::
       FILE ${PROJECT_NAME}Targets.cmake)
export(PACKAGE ${PROJECT_NAME})
