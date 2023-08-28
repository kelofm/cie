#ifndef CIE_GL_ABS_WINDOW_IMPL_HPP
#define CIE_GL_ABS_WINDOW_IMPL_HPP

namespace cie::gl {


template < class SceneType,
           class ...Args >
inline std::shared_ptr<SceneType>
AbsWindow::makeScene( Args&&... args )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto localBlock = this->newBlock( "construct new scene" );

    auto p_scene = std::make_shared<SceneType>( std::forward<Args>(args)... );

    this->log( p_scene->name() );

    this->addScene( p_scene );

    return p_scene;

    CIE_END_EXCEPTION_TRACING
}


inline void AbsWindow::triggerMouseButtonCallback(KeyEnum key, KeyEnum action, KeyEnum modifier)
{
    _mouseButtonCallback(key, action, modifier);
}


inline void AbsWindow::triggerCursorPositionCallback(Double x, Double y)
{
    _cursorPositionCallback(x, y);
}


inline void AbsWindow::triggerCursorEnterCallback(KeyEnum hasEntered)
{
    _cursorEnterCallback(hasEntered);
}


inline void AbsWindow::triggerScrollCallback(Double xOffset, Double yOffset)
{
    _scrollCallback(xOffset, yOffset);
}


inline void AbsWindow::triggerKeyboardCallback(KeyEnum key, KeyEnum scancode, KeyEnum action, KeyEnum modifier)
{
    _keyboardCallback(key, scancode, action, modifier);
}


inline void AbsWindow::triggerWindowResizeCallback(Size width, Size height)
{
    _windowResizeCallback(width, height);
}


} // namespace cie::gl

#endif