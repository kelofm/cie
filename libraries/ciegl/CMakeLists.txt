cmake_minimum_required(VERSION 3.15.0)
project(ciegl LANGUAGES C CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)

# Declare all CMake settings
set(${PROJECT_NAME_UPPER}_BUILD_SHARED_LIBRARY ON CACHE BOOL "build shared libraries instead of static ones")
set(${PROJECT_NAME_UPPER}_BUILD_TESTS OFF CACHE BOOL "build testrunner for ${PROJECT_NAME}")

# Language server
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Collect dependencies
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}/cie/lib")
find_package(cieutils REQUIRED)
find_package(linalg REQUIRED)
find_package(geo REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)

# Build GLAD
add_library(glad SHARED "${CMAKE_CURRENT_SOURCE_DIR}/external/glad/src/glad.c")
target_include_directories(glad
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/glad/inc>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/external/glad/inc>")
set_target_properties(glad PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
target_compile_definitions(glad PUBLIC GLAD_GLAPI_EXPORT PRIVATE GLAD_GLAPI_EXPORT_BUILD)
install(TARGETS glad EXPORT ${PROJECT_NAME}Targets LIBRARY DESTINATION "cie/lib")

# Build main library
file(GLOB_RECURSE sources "${CMAKE_CURRENT_SOURCE_DIR}/packages/*/src/*.cpp")
if (${${PROJECT_NAME_UPPER}_BUILD_SHARED_LIBRARY})
    add_library(${PROJECT_NAME} SHARED ${sources})
    set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE 1)
else()
    add_library(${PROJECT_NAME} STATIC ${sources})
endif()

# GLM removed its legacy cmake target, but then readded it? idk
if(TARGET glm)
    set(libglm "glm")
elseif(TARGET glm::glm)
    set(libglm "glm::glm")
else()
    message(FATAL_ERROR "Neither glm nor glm::glm is a valid target")
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC cie::geo glfw OpenGL::GL ${libglm} glad)
target_include_directories(${PROJECT_NAME}
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_include_directories(${PROJECT_NAME}
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/stb>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/external/stb>")
set_target_properties(${PROJECT_NAME} PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                      INSTALL_RPATH "$ORIGIN")
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets LIBRARY DESTINATION "cie/lib")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/packages" DESTINATION "cie/include" FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/external" DESTINATION "cie/include" FILES_MATCHING PATTERN "*.h*")

# Build test executable
if (${PROJECT_NAME_UPPER}_BUILD_TESTS)
    file(GLOB_RECURSE test_sources "${CMAKE_CURRENT_SOURCE_DIR}/packages/*/test/*.cpp")
    list(INSERT test_sources 0 "${CMAKE_CURRENT_SOURCE_DIR}/test/${PROJECT_NAME}_testrunner.cpp")

    add_executable(${PROJECT_NAME}_testrunner ${test_sources})
    target_link_libraries(${PROJECT_NAME}_testrunner PRIVATE ${PROJECT_NAME})
    set_target_properties(${PROJECT_NAME}_testrunner PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

# Build python module
if (${PROJECT_NAME_UPPER}_BUILD_PYTHON_MODULE)
    # Each subdirectory in modules is assumed to be a python module
    # and is installed as-is
    file(GLOB children "${CMAKE_CURRENT_SOURCE_DIR}/python/modules/*")
    foreach(child ${children})
        if (IS_DIRECTORY "${child}")
            get_filename_component(module_name "${child}" NAME)
            install(DIRECTORY "${child}" DESTINATION "cie")
        endif()
    endforeach()
endif()

# Install data
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data/" DESTINATION "${CMAKE_INSTALL_PREFIX}/cie/data/")

# Export
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(TARGETS glad ${PROJECT_NAME} NAMESPACE cie:: FILE ${PROJECT_NAME}Targets.cmake)
export(PACKAGE ${PROJECT_NAME})
install(EXPORT ${PROJECT_NAME}Targets FILE ${PROJECT_NAME}Targets.cmake DESTINATION "cie/lib")
install(EXPORT ${PROJECT_NAME}Targets FILE ${PROJECT_NAME}Targets.cmake DESTINATION "cie/lib/cmake/${PROJECT_NAME}")

# Package
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

include(CMakePackageConfigHelpers)
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION "cie/lib/cmake/${PROJECT_NAME}")
