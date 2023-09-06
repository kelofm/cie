// --- External Includes ---
#include "glad/glad.h"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- Graphics Includes ---
#include "packages/control/inc/callback_types.hpp"
#include "packages/context/inc/AbsWindow.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::gl {


/* --- DEFAULT CALLBACKS --- */
namespace detail {

void defaultMouseButtonCallback(KeyEnum,
                                KeyEnum,
                                KeyEnum)
{
}

void defaultCursorPositionCallback(double,
                                   double)
{
}

void defaultCursorEnterCallback(KeyEnum)
{
}

void defaultScrollCallback(double,
                           double)
{
}

void defaultKeyboardCallback(KeyEnum,
                             KeyEnum,
                             KeyEnum,
                             KeyEnum)
{
}

} // namespace detail


AbsWindow::AbsWindow(Size id,
                     const std::string& r_name,
                     Size width,
                     Size height,
                     utils::Logger& r_logger) :
    utils::IDObject<Size>(id),
    utils::Loggee(r_logger, r_name),
    utils::observer::Subject(),
    _size(width, height),
    _scenes(),
    _continueLooping(false),
    _mouseButtonCallback(detail::defaultMouseButtonCallback),
    _cursorPositionCallback(detail::defaultCursorPositionCallback),
    _cursorEnterCallback(detail::defaultCursorEnterCallback),
    _scrollCallback(detail::defaultScrollCallback),
    _keyboardCallback(detail::defaultKeyboardCallback),
    _windowResizeCallback()
{
    // Default window resize callback
    _windowResizeCallback = [this](int width, int height){this->setSize(width, height);};
}


AbsWindow::~AbsWindow()
{
    auto localBlock = this->newBlock("close window");
    this->logID("", this->getID());
}


bool AbsWindow::update()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->activate();

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    this->notifyObservers();

    bool output = true;
    for (auto& rp_scene : _scenes)
        output = output && rp_scene->update();

    output = output && this->updateImpl();

    if (!output) [[unlikely]]
        _continueLooping = false;
    return output;

    CIE_END_EXCEPTION_TRACING
}


void AbsWindow::setSize(WindowSize size)
{
    CIE_BEGIN_EXCEPTION_TRACING
    this->setSizeImpl(size);
    _size = size;
    glViewport(0, 0, size.first, size.second);
    CIE_END_EXCEPTION_TRACING
}


void AbsWindow::setSize(Size width,
                        Size height)
{
    this->setSize(std::make_pair(width, height));

    // Update the aspect ratios of all cameras
    for (auto& rp_scene : _scenes)
        for (auto& rp_camera : rp_scene->cameras())
            rp_camera->setAspectRatio(this->getAspectRatio());

    // Trigger a redraw
    this->update();
}


const std::pair<Size,Size>& AbsWindow::getSize() const
{
    return _size;
}


double AbsWindow::getAspectRatio() const
{
    CIE_DIVISION_BY_ZERO_CHECK(_size.second != 0)
    return double(_size.first) / double(_size.second);
}


void AbsWindow::beginLoop()
{
    if (this->_continueLooping)
        this->log("Attempt to begin loop while looping",
                   LOG_TYPE_WARNING);

    this->_continueLooping = true;

    auto localBlock = this->newBlock("Event loop");

    while(this->_continueLooping)
    {this->update();}
}


void AbsWindow::endLoop()
{
    this->_continueLooping = false;
}


void AbsWindow::screenshot(const std::filesystem::path& r_outputPath,
                            Image* p_targetImage) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    bool hasTargetImage = false;

    if (!p_targetImage)
        p_targetImage = new Image(
            this->_size.first,
            this->_size.second,
            3
        );
    else
    {
        CIE_OUT_OF_RANGE_CHECK(p_targetImage->width() == this->_size.first)
        CIE_OUT_OF_RANGE_CHECK(p_targetImage->height() == this->_size.second)
        CIE_OUT_OF_RANGE_CHECK(p_targetImage->numberOfChannels() == 3)
        hasTargetImage = true;
    }

    glReadPixels(
        0, 0,
        this->_size.first, this->_size.second,
        GL_RGB,
        GL_UNSIGNED_BYTE,
        p_targetImage->data()
    );

    p_targetImage->verticalFlip();
    p_targetImage->write(r_outputPath);

    if (!hasTargetImage)
        delete p_targetImage;

    CIE_END_EXCEPTION_TRACING
}


void AbsWindow::addScene(ScenePtr p_scene, bool update)
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER(p_scene)

    auto it_scene = std::find(
        _scenes.begin(),
        _scenes.end(),
        p_scene
    );

    if (it_scene != _scenes.end())
        CIE_THROW(Exception, "Attempt to add existing Scene to Window")

    _scenes.push_back(p_scene);

    if (update)
        this->update();

    CIE_END_EXCEPTION_TRACING
}


void AbsWindow::removeScene(ScenePtr p_scene)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto it_scene = std::find(
        _scenes.begin(),
        _scenes.end(),
        p_scene
    );

    if (it_scene == _scenes.end())
        CIE_THROW(Exception, "Attempt to remove non-existent Scene from Window")

    CIE_CHECK(it_scene != _scenes.end(),
              "Attempt to remove scene '" << p_scene->name() << "' from window '" << this->name() << "' which does not own it")

    _scenes.erase(it_scene);

    CIE_END_EXCEPTION_TRACING
}


const typename AbsWindow::scene_container&
AbsWindow::scenes() const
{
    return _scenes;
}


void AbsWindow::setMouseButtonCallback(MouseButtonCallback function)
{this->_mouseButtonCallback = function;}


void AbsWindow::setCursorPositionCallback(CursorPositionCallback function)
{this->_cursorPositionCallback = function;}


void AbsWindow::setCursorEnterCallback(CursorEnterCallback function)
{this->_cursorEnterCallback = function;}


void AbsWindow::setScrollCallback(ScrollCallback function)
{this->_scrollCallback = function;}


void AbsWindow::setKeyboardCallback(KeyboardCallback function)
{this->_keyboardCallback = function;}


void AbsWindow::setWindowResizeCallback(WindowResizeCallback function)
{this->_windowResizeCallback = function;}


void AbsWindow::setName(std::string::const_iterator begin,
                        std::string::const_iterator end)
{
    CIE_BEGIN_EXCEPTION_TRACING

    utils::Loggee::setName(begin, end);
    this->setWindowName(begin, end);

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::gl