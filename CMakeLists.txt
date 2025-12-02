cmake_minimum_required(VERSION 3.15.0)
project(cie LANGUAGES CXX VERSION 0.0.0)

set(library_names "cieutils;linalg;geo;fem")
foreach(project_name ${library_names})
    string(TOUPPER ${project_name} project_name_upper)
    set(CIE_ENABLE_${project_name_upper} OFF CACHE BOOL "Configure ${project_name}.")

    if(CIE_ENABLE_${project_name_upper} MATCHES ON)
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/libraries/${project_name}")
    endif()
endforeach()

set(executable_names "benchmarks;bad_apple")
foreach(project_name ${executable_names})
    string(TOUPPER ${project_name} project_name_upper)
    set(CIE_ENABLE_${project_name_upper} OFF CACHE BOOL "Configure ${project_name}.")

    if(CIE_ENABLE_${project_name_upper} MATCHES ON)
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/libraries/${project_name}")
    endif()
endforeach()
