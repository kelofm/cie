// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/logging/inc/StreamInterceptor.hpp"

// --- STL Includes ---
#include <sstream>


namespace cie {


CIE_TEST_CASE("StreamInterceptor", "[logging]")
{
    CIE_TEST_CASE_INIT("StreamInterceptor")

    std::stringstream source, sink;
    auto p_source = OutputStream::SharedPointer(
        new OutputStream(&source, false)
    );
    auto p_sink = OutputStream::SharedPointer(
        new OutputStream(&sink, false)
    );

    {
        auto p_redirector = OutputStream::SharedPointer(
            new StreamInterceptor(p_source, p_sink)
        );
        CIE_TEST_CHECK_NOTHROW(source << "abc");
    }

    std::string sourceString = source.str();
    std::string sinkString   = sink.str();

    CIE_TEST_CHECK(sourceString.empty());
    CIE_TEST_CHECK(sinkString.size() == 3);
    CIE_TEST_CHECK(sinkString == "abc");
}


} // namespace cie