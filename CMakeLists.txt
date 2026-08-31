cmake_minimum_required(VERSION 3.20)
project(cie LANGUAGES CXX VERSION 0.0.0)

set(CMAKE_INSTALL_MESSAGE "LAZY" CACHE STRING "")

if(CIE_ENABLE_GEO)
    set(CIE_ENABLE_LINALG ON CACHE BOOL "Configure linalg." FORCE)
endif()

if(CIE_ENABLE_LINALG)
    set(CIE_ENABLE_CIEUTILS ON CACHE BOOL "Configure cieutils." FORCE)
endif()

if(CIE_ENABLE_FEM)
    set(CIE_ENABLE_GEO ON CACHE BOOL "Configure geo." FORCE)
endif()

set(library_names "cieutils;linalg;geo;fem")
foreach(project_name ${library_names})
    string(TOUPPER ${project_name} project_name_upper)
    set(CIE_ENABLE_${project_name_upper} OFF CACHE BOOL "Configure ${project_name}.")

    if(CIE_ENABLE_${project_name_upper} MATCHES ON)
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/libraries/${project_name}")
    endif()
endforeach()

set(executable_names "benchmarks;bad_apple;fem_examples;geo_examples")
foreach(project_name ${executable_names})
    string(TOUPPER ${project_name} project_name_upper)
    set(CIE_ENABLE_${project_name_upper} OFF CACHE BOOL "Configure ${project_name}.")

    if(CIE_ENABLE_${project_name_upper} MATCHES ON)
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/executables/${project_name}")
    endif()
endforeach()
