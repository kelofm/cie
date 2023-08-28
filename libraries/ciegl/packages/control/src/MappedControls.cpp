// --- External Includes ---
#include "glad/glad.h"
#include "GLFW/glfw3.h"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/io/inc/json.hpp"

// --- Internal Includes ---
#include "packages/control/inc/MappedControls.hpp"

// --- STL Includes ---
#include <functional>


namespace cie::gl {


std::unordered_map<std::string,KeyEnum> GLFW_KEY_ENUM_MAP = {
    std::make_pair("space",                GLFW_KEY_SPACE),
    std::make_pair("'",                    GLFW_KEY_APOSTROPHE),
    std::make_pair(",",                    GLFW_KEY_COMMA),
    std::make_pair("-",                    GLFW_KEY_MINUS),
    std::make_pair(".",                    GLFW_KEY_PERIOD),
    std::make_pair("/",                    GLFW_KEY_SLASH),
    std::make_pair("0",                    GLFW_KEY_0),
    std::make_pair("1",                    GLFW_KEY_1),
    std::make_pair("2",                    GLFW_KEY_2),
    std::make_pair("3",                    GLFW_KEY_3),
    std::make_pair("4",                    GLFW_KEY_4),
    std::make_pair("5",                    GLFW_KEY_5),
    std::make_pair("6",                    GLFW_KEY_6),
    std::make_pair("7",                    GLFW_KEY_7),
    std::make_pair("8",                    GLFW_KEY_8),
    std::make_pair("9",                    GLFW_KEY_9),
    std::make_pair(";",                    GLFW_KEY_SEMICOLON),
    std::make_pair("=",                    GLFW_KEY_EQUAL),
    std::make_pair("a",                    GLFW_KEY_A),
    std::make_pair("b",                    GLFW_KEY_B),
    std::make_pair("c",                    GLFW_KEY_C),
    std::make_pair("d",                    GLFW_KEY_D),
    std::make_pair("e",                    GLFW_KEY_E),
    std::make_pair("f",                    GLFW_KEY_F),
    std::make_pair("g",                    GLFW_KEY_G),
    std::make_pair("h",                    GLFW_KEY_H),
    std::make_pair("i",                    GLFW_KEY_I),
    std::make_pair("j",                    GLFW_KEY_J),
    std::make_pair("k",                    GLFW_KEY_K),
    std::make_pair("l",                    GLFW_KEY_L),
    std::make_pair("m",                    GLFW_KEY_M),
    std::make_pair("n",                    GLFW_KEY_N),
    std::make_pair("o",                    GLFW_KEY_O),
    std::make_pair("p",                    GLFW_KEY_P),
    std::make_pair("q",                    GLFW_KEY_Q),
    std::make_pair("r",                    GLFW_KEY_R),
    std::make_pair("s",                    GLFW_KEY_S),
    std::make_pair("t",                    GLFW_KEY_T),
    std::make_pair("u",                    GLFW_KEY_U),
    std::make_pair("v",                    GLFW_KEY_V),
    std::make_pair("w",                    GLFW_KEY_W),
    std::make_pair("x",                    GLFW_KEY_X),
    std::make_pair("y",                    GLFW_KEY_Y),
    std::make_pair("z",                    GLFW_KEY_Z),
    std::make_pair("[",                    GLFW_KEY_LEFT_BRACKET),
    std::make_pair("\\",                   GLFW_KEY_BACKSLASH),
    std::make_pair("]",                    GLFW_KEY_RIGHT_BRACKET),
    std::make_pair("`",                    GLFW_KEY_GRAVE_ACCENT),
    std::make_pair("esc",                  GLFW_KEY_ESCAPE),
    std::make_pair("enter",                GLFW_KEY_ENTER),
    std::make_pair("tab",                  GLFW_KEY_TAB),
    std::make_pair("backspace",            GLFW_KEY_BACKSPACE),
    std::make_pair("insert",               GLFW_KEY_INSERT),
    std::make_pair("delete",               GLFW_KEY_DELETE),
    std::make_pair("right",                GLFW_KEY_RIGHT),
    std::make_pair("left",                 GLFW_KEY_LEFT),
    std::make_pair("down",                 GLFW_KEY_DOWN),
    std::make_pair("up",                   GLFW_KEY_UP),
    std::make_pair("page_up",              GLFW_KEY_PAGE_UP),
    std::make_pair("page_down",            GLFW_KEY_PAGE_DOWN),
    std::make_pair("home",                 GLFW_KEY_HOME),
    std::make_pair("end",                  GLFW_KEY_END),
    std::make_pair("caps_lock",            GLFW_KEY_CAPS_LOCK),
    std::make_pair("scroll_lock",          GLFW_KEY_SCROLL_LOCK),
    std::make_pair("num_lock",             GLFW_KEY_NUM_LOCK),
    std::make_pair("print_screen",         GLFW_KEY_PRINT_SCREEN),
    std::make_pair("key_pause",            GLFW_KEY_PAUSE),
    std::make_pair("f1",                   GLFW_KEY_F1),
    std::make_pair("f2",                   GLFW_KEY_F2),
    std::make_pair("f3",                   GLFW_KEY_F3),
    std::make_pair("f4",                   GLFW_KEY_F4),
    std::make_pair("f5",                   GLFW_KEY_F5),
    std::make_pair("f6",                   GLFW_KEY_F6),
    std::make_pair("f7",                   GLFW_KEY_F7),
    std::make_pair("f8",                   GLFW_KEY_F8),
    std::make_pair("f9",                   GLFW_KEY_F9),
    std::make_pair("f10",                  GLFW_KEY_F10),
    std::make_pair("f11",                  GLFW_KEY_F11),
    std::make_pair("f12",                  GLFW_KEY_F12),
    std::make_pair("f13",                  GLFW_KEY_F13),
    std::make_pair("f14",                  GLFW_KEY_F14),
    std::make_pair("f15",                  GLFW_KEY_F15),
    std::make_pair("f16",                  GLFW_KEY_F16),
    std::make_pair("f17",                  GLFW_KEY_F17),
    std::make_pair("f18",                  GLFW_KEY_F18),
    std::make_pair("f19",                  GLFW_KEY_F19),
    std::make_pair("f20",                  GLFW_KEY_F20),
    std::make_pair("f21",                  GLFW_KEY_F21),
    std::make_pair("f22",                  GLFW_KEY_F22),
    std::make_pair("f23",                  GLFW_KEY_F23),
    std::make_pair("f24",                  GLFW_KEY_F24),
    std::make_pair("f25",                  GLFW_KEY_F25),
    std::make_pair("kp_0",                 GLFW_KEY_KP_0),
    std::make_pair("kp_1",                 GLFW_KEY_KP_1),
    std::make_pair("kp_2",                 GLFW_KEY_KP_2),
    std::make_pair("kp_3",                 GLFW_KEY_KP_3),
    std::make_pair("kp_4",                 GLFW_KEY_KP_4),
    std::make_pair("kp_5",                 GLFW_KEY_KP_5),
    std::make_pair("kp_6",                 GLFW_KEY_KP_6),
    std::make_pair("kp_7",                 GLFW_KEY_KP_7),
    std::make_pair("kp_8",                 GLFW_KEY_KP_8),
    std::make_pair("kp_9",                 GLFW_KEY_KP_9),
    std::make_pair("kp_decimal",           GLFW_KEY_KP_DECIMAL),
    std::make_pair("kp_divide",            GLFW_KEY_KP_DIVIDE),
    std::make_pair("kp_multiply",          GLFW_KEY_KP_MULTIPLY),
    std::make_pair("kp_subtract",          GLFW_KEY_KP_SUBTRACT),
    std::make_pair("kp_add",               GLFW_KEY_KP_ADD),
    std::make_pair("kp_enter",             GLFW_KEY_KP_ENTER),
    std::make_pair("kp_equal",             GLFW_KEY_KP_EQUAL),
    std::make_pair("left_shift",           GLFW_KEY_LEFT_SHIFT),
    std::make_pair("left_control",         GLFW_KEY_LEFT_CONTROL),
    std::make_pair("left_alt",             GLFW_KEY_LEFT_ALT),
    std::make_pair("left_super",           GLFW_KEY_LEFT_SUPER),
    std::make_pair("right_shift",          GLFW_KEY_RIGHT_SHIFT),
    std::make_pair("right_control",        GLFW_KEY_RIGHT_CONTROL),
    std::make_pair("right_alt",            GLFW_KEY_RIGHT_ALT),
    std::make_pair("right_super",          GLFW_KEY_RIGHT_SUPER),
    std::make_pair("menu",                 GLFW_KEY_MENU)
};


MappedControls::MappedControls() :
    _controlMap()
{
    _controlMap = MappedControls::map_type {
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_SPACE,           nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_APOSTROPHE,      nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_COMMA,           nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_MINUS,           nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_PERIOD,          nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_SLASH,           nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_0,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_1,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_2,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_3,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_4,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_5,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_6,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_7,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_8,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_9,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_SEMICOLON,       nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_EQUAL,           nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_A,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_B,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_C,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_D,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_E,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_G,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_H,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_I,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_J,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_K,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_L,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_M,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_N,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_O,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_P,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_Q,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_R,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_S,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_T,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_U,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_V,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_W,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_X,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_Y,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_Z,               nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_LEFT_BRACKET,    nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_BACKSLASH,       nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_RIGHT_BRACKET,   nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_GRAVE_ACCENT,    nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_ESCAPE,          nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_ENTER,           nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_TAB,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_BACKSPACE,       nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_INSERT,          nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_DELETE,          nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_RIGHT,           nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_LEFT,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_DOWN,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_UP,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_PAGE_UP,         nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_PAGE_DOWN,       nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_HOME,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_END,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_CAPS_LOCK,       nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_SCROLL_LOCK,     nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_NUM_LOCK,        nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_PRINT_SCREEN,    nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_PAUSE,           nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F1,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F2,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F3,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F4,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F5,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F6,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F7,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F8,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F9,              nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F10,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F11,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F12,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F13,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F14,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F15,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F16,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F17,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F18,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F19,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F20,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F21,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F22,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F23,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F24,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_F25,             nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_0,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_1,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_2,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_3,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_4,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_5,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_6,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_7,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_8,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_9,            nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_DECIMAL,      nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_DIVIDE,       nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_MULTIPLY,     nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_SUBTRACT,     nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_ADD,          nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_ENTER,        nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_KP_EQUAL,        nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_LEFT_SHIFT,      nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_LEFT_CONTROL,    nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_LEFT_ALT,        nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_LEFT_SUPER,      nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_RIGHT_SHIFT,     nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_RIGHT_CONTROL,   nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_RIGHT_ALT,       nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_RIGHT_SUPER,     nullptr),
        std::make_pair<KeyEnum,MappedControls::control_function>(GLFW_KEY_MENU,            nullptr)
    };
}


void MappedControls::configure(MappedControls::configuration_map_type& r_configMap)
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (auto& r_pair : r_configMap)
        for (auto key : r_pair.first.keys())
        {
            // Check whether a control is already assigned to this key
            auto it_map = this->_controlMap.find(key);

            if (it_map == this->_controlMap.end())
                CIE_THROW(Exception, "Cannot find '" + std::to_string(key) + "' is not found in the control map")

            if (it_map->second)
                CIE_THROW(Exception, "Control '" + std::to_string(key) + "' is overwritten")

            // Assign control
            it_map->second = r_pair.second;
        }

    CIE_END_EXCEPTION_TRACING
}


void MappedControls::configure(std::istream& r_configurationFile)
{
    CIE_BEGIN_EXCEPTION_TRACING

    io::JSONObject configuration(r_configurationFile);

    // Populate map from json
    MappedControls::configuration_contents contents;

    const auto it_end = configuration.end();
    for (auto it_node=configuration.begin(); it_node!=it_end; ++it_node)
    {
        auto& r_pair = contents.emplace_back(
            it_node.key(),
            ControlKey()
        );

        const auto node = *it_node;
        const Size controlSize = node.size();
        for (Size controlIndex=0; controlIndex<controlSize; ++controlIndex)
        {
            const auto key = node[controlIndex].as<std::string>();
            auto it_keyMap = GLFW_KEY_ENUM_MAP.find(key);
            if (it_keyMap == GLFW_KEY_ENUM_MAP.end())
                CIE_THROW(Exception, "Key '" + std::string(key) + "' not found in map")

            r_pair.second.add(it_keyMap->second);
        }
    }

    // Map function pointers to controls
    auto configurationMap = this->makeConfigurationMap(contents);
    this->configure(configurationMap);

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl