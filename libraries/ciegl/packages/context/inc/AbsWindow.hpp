#ifndef CIE_GL_ABS_WINDOW_HPP
#define CIE_GL_ABS_WINDOW_HPP

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/types/inc/NamedObject.hpp"
#include "packages/types/inc/IDObject.hpp"
#include "packages/logging/inc/Loggee.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/observer/inc/Subject.hpp"

// --- Internal Includes ---
#include "packages/scene/inc/Scene.hpp"
#include "packages/control/inc/callback_types.hpp"
#include "packages/texture/inc/Image.hpp"

// --- STL Includes ---
#include <memory>
#include <utility>
#include <list>
#include <functional>
#include <filesystem>


namespace cie::gl {


class AbsWindow :
    public utils::IDObject<Size>,
    public utils::Loggee,
    public utils::observer::Subject
{
private:
    using scene_container = std::list<ScenePtr>;

public:
    /// {key, action, modifier}
    using key_callback_function = std::function<void(KeyEnum,KeyEnum,KeyEnum)>;

    /// {x, y, button, action, modifier}
    using mouse_callback_function = std::function<void(Double,Double,KeyEnum,KeyEnum,KeyEnum)>;

    /// {width, height}
    using WindowSize = std::pair<Size,Size>;

    CIE_DEFINE_CLASS_POINTERS(AbsWindow)

public:
    CIE_DEFINE_CLASS_DEFAULT_MOVES(AbsWindow)

    virtual ~AbsWindow();

    /// Update all scenes and draw
    virtual bool update();

    /// Make window active
    virtual void activate() = 0;

    void setSize(WindowSize size);

    void setSize(Size width,
                 Size height);

    const std::pair<Size,Size>& getSize() const;

    double getAspectRatio() const;

    /// Construct a new scene and add it to this window
    template < class SceneType,
               class ...Args>
    std::shared_ptr<SceneType> makeScene(Args&&... args);

    /// Attempt to remove a scene from this window
    virtual void removeScene(ScenePtr p_scene);

    const scene_container& scenes() const;

    void beginLoop();

    void endLoop();

    /** Copy window data into an RGB image and write it to disk
     *  @param r_outputPath relative or full path of the new image file
     *  @param p_targetImage preallocated image. If none is given, allocation
     *  is performed internally
     */
    void screenshot(const std::filesystem::path& r_outputPath,
                    Image* p_targetImage = nullptr) const;

    void setMouseButtonCallback(MouseButtonCallback function);

    void setCursorPositionCallback(CursorPositionCallback function);

    void setCursorEnterCallback(CursorEnterCallback function);

    void setScrollCallback(ScrollCallback function);

    void setKeyboardCallback(KeyboardCallback function);

    void setWindowResizeCallback(WindowResizeCallback function);

    virtual void setName(std::string::const_iterator begin,
                         std::string::const_iterator end) override;

    using utils::Loggee::setName;

protected:
    AbsWindow(Size id,
              const std::string& r_name,
              Size width,
              Size height,
              utils::Logger& r_logger = utils::LoggerSingleton::get());

    /// Attempt to add a scene to this window
    virtual void addScene(ScenePtr p_scene, bool update = true);

    virtual void setSizeImpl(WindowSize size) = 0;

    virtual bool updateImpl() = 0;

    virtual void setWindowName(std::string::const_iterator begin,
                               std::string::const_iterator end) = 0;

    void triggerMouseButtonCallback(KeyEnum key, KeyEnum action, KeyEnum modifier);

    void triggerCursorPositionCallback(Double x, Double y);

    void triggerCursorEnterCallback(KeyEnum hasEntered);

    void triggerScrollCallback(Double xOffset, Double yOffset);

    void triggerKeyboardCallback(KeyEnum key, KeyEnum scancode, KeyEnum action, KeyEnum modifier);

    void triggerWindowResizeCallback(Size width, Size height);

protected:
    std::pair<Size,Size>    _size;

    scene_container         _scenes;

    bool                    _continueLooping;

private:
    MouseButtonCallback     _mouseButtonCallback;

    CursorPositionCallback  _cursorPositionCallback;

    CursorEnterCallback     _cursorEnterCallback;

    ScrollCallback          _scrollCallback;

    KeyboardCallback        _keyboardCallback;

    WindowResizeCallback    _windowResizeCallback;
};


} // namespace cie::gl

#include "packages/context/impl/AbsWindow_impl.hpp"

#endif