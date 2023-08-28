// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/logging/inc/Logger.hpp"
#include "cmake_variables.hpp"

// --- Internal Includes ---
#include "packages/camera/inc/OrthographicCamera.hpp"
#include "packages/camera/inc/PerspectiveCamera.hpp"

// --- STL Includes ---
#include <iostream>


namespace cie::gl {


CIE_TEST_CASE("Camera", "[camera]")
{
    CIE_TEST_CASE_INIT("Camera")

    utils::Logger logger(getOutputPath()/ "Camera_test.txt");

    {
        CIE_TEST_CASE_INIT("OrthographicCamera")

        using CameraType = OrthographicCamera;
        AbsCamera::SharedPointer p_camera;

        CIE_TEST_REQUIRE_NOTHROW(p_camera = AbsCamera::SharedPointer(new CameraType(logger)));
    }

    {
        CIE_TEST_CASE_INIT("PerspectiveCamera")

        using CameraType = PerspectiveCamera;
        AbsCamera::SharedPointer p_camera;

        CIE_TEST_REQUIRE_NOTHROW(p_camera = AbsCamera::SharedPointer(new CameraType(logger)));
    }
}


} // namespace cie::gl