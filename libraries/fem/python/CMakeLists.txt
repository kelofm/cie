if (${CIE_BUILD_PYTHON_MODULES})
    set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE 1)
    get_python_binding_name(binding_library)
    build_python_bindings()
    register_python_bindings(${binding_library})

    copy_python_modules()
    copy_scripts()
endif()
