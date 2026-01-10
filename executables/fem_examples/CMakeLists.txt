cmake_minimum_required(VERSION 3.15.0)
project(fem_examples LANGUAGES CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)

# Language server
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Collect dependencies
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}/cie/lib/cmake")
if (NOT TARGET fem)
    find_package(fem CONFIG REQUIRED)
endif()

# Build executables.
file(GLOB drivers "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
set(targets "")
foreach(driver ${drivers})
    # Extract target name and define executable.
    get_filename_component(driver_name ${driver} NAME_WE)
    add_executable(${driver_name} "${driver}")
    target_include_directories(${driver_name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    list(APPEND targets ${driver_name})

    # Configure the compiler.
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(${driver_name} PRIVATE
                               -Wall
                               -Wpedantic
                               -Wextra
                               -Wno-dangling-reference
                               -Wno-maybe-uninitialized
                               -Wno-uninitialized
                               -Wno-array-bounds           #< spurious warning for 1-item copies
                               -Wno-stringop-overflow      #< spurious warning for 1-item copies
                               -Werror)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_compile_options(${driver_name} PRIVATE
                               -Wall
                               -Wpedantic
                               -Wextra
                               -Werror)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        target_compile_options(${driver_name} PRIVATE
                               -Wall
                               -Wpedantic
                               -Wextra
                               -Werror)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        message(FATAL_ERROR "MSVC is not supported yet")
    endif ()

    # Link dependencies.
    target_link_libraries(${driver_name} PRIVATE cie::fem)

    # Installation.
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set_target_properties(
            ${driver_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            INSTALL_RPATH "@loader_path/../lib")
    else()
        set_target_properties(
            ${driver_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            INSTALL_RPATH "$ORIGIN/../lib")
    endif()

    install(TARGETS ${driver_name} RUNTIME DESTINATION "cie/bin")

    # Package.
    include(InstallRequiredSystemLibraries)
    set(CPACK_PACKAGE_VENDOR "")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "")
    set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
    set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
    set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/readme.md")

    set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
    set(CPACK_SOURCE_IGNORE_FILES
        /.git
        /build
        *.gitignore)
    include(CPack)
endforeach()

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
install(TARGETS ${targets}
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
export(TARGETS ${targets}
       NAMESPACE cie::
       FILE ${PROJECT_NAME}Targets.cmake)
export(PACKAGE ${PROJECT_NAME})
