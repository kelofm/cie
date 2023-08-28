// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/context/inc/AbsContext.hpp"

// --- STL Includes ---
#include <algorithm>

namespace cie::gl {


AbsContext::AbsContext( Size versionMajor,
                        Size versionMinor,
                        Size MSAASamples,
                        utils::LoggerPtr p_logger ) :
    _version( versionMajor, versionMinor ),
    _MSAASamples( MSAASamples ),
    _windowCounter( 0 ),
    _p_logger( std::move(p_logger) )
{
}


AbsContext::~AbsContext()
{
    _p_logger->log( "Destroy context" );
}


const std::pair<Size,Size>&
AbsContext::version() const
{
    return _version;
}


const Size
AbsContext::MSAASamples() const
{
    return _MSAASamples;
}


AbsWindow::SharedPointer AbsContext::newWindow( Size width,
                                 Size height,
                                 const std::string& r_name )
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto scopedBlock = _p_logger->newBlock("new window '" + r_name);

    auto p_window = this->newWindowImpl(width,
                                        height,
                                        r_name );

    p_window->setID( this->_windowCounter++ );

    _p_logger->logs( "ID_", p_window->getID() );

    return p_window;

    CIE_END_EXCEPTION_TRACING
}


utils::Logger& AbsContext::logger()
{
    return *this->_p_logger;
}


} // namespace cie::gl