#ifndef CIE_GL_CALLBACK_TYPES_HPP
#define CIE_GL_CALLBACK_TYPES_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <functional>


namespace cie::gl {


/**
 * Enum type used by GLFW
 */
using KeyEnum = int;

/**
 * button (mouse button pressed or released)
 * action (GLFWPRESS or GLFWRELEASE)
 * modifiers (bit field for modifier keys)
 *      GLFW_MOD_SHIFT
 *      GLFW_MOD_CONTROL
 *      GLFW_MOD_ALT
 *      GLFW_MOD_SUPER
 */
using MouseButtonCallback = std::function<void(KeyEnum,KeyEnum,KeyEnum)>;

/**
 * window x position (range: [])
 * window y position (range: [])
 */
using CursorPositionCallback = std::function<void(Double,Double)>;

/**
 * Callback for cursor entering/leaving window area
 *  entered (GL_TRUE for entering and GL_FALSE for leaving)
 */
using CursorEnterCallback = std::function<void(KeyEnum)>;

/**
 * x offset
 * y offset
 */
using ScrollCallback = std::function<void(Double,Double)>;

/**
 * key (https://www.glfw.org/docs/3.0/group__keys.html)
 * scancode (platform-specific key identifiers)
 * action (GLFW_PRESS, GLFW_RELEASE, GLFW_REPEAT)
 * modifiers (bit field for modifier keys)
 *      GLFW_MOD_SHIFT
 *      GLFW_MOD_CONTROL
 *      GLFW_MOD_ALT
 *      GLFW_MOD_SUPER
 */
using KeyboardCallback = std::function<void(KeyEnum,KeyEnum,KeyEnum,KeyEnum)>;

/**
 * width
 * height
 */
using WindowResizeCallback = std::function<void(Size,Size)>;


} // namespace cie::gl


#endif