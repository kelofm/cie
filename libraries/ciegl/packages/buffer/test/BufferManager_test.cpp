// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/context/inc/GLFWContext.hpp"
#include "packages/context/inc/GLFWWindow.hpp"
#include "packages/buffer/inc/BufferManager.hpp"


namespace cie::gl {


CIE_TEST_CASE("BufferManager", "[buffer]")
{
    CIE_TEST_CASE_INIT("BufferManager")

    auto p_log = std::make_shared<utils::Logger>(getOutputPath()/ "BufferManager_test.txt");

    // Context
    auto p_context = GLFWContextSingleton::get(p_log);

    // Window
    AbsWindow::SharedPointer p_window;
    CIE_CHECK_POINTER(p_context.lock())
    CIE_TEST_REQUIRE_NOTHROW(p_window = p_context.lock()->newWindow());

    // Buffers
    BufferManager::SharedPointer p_bufferManager;
    VertexBuffer::SharedPointer p_vertexBuffer;
    ElementBuffer::SharedPointer p_elementBuffer;

    {
        CIE_TEST_CASE_INIT("create buffer manager")
        CIE_TEST_REQUIRE_NOTHROW(p_bufferManager.reset(new BufferManager(p_context.lock()->logger())));
    }

    {
        CIE_TEST_CASE_INIT("generate buffers")

        CIE_TEST_CHECK_NOTHROW(p_vertexBuffer = p_bufferManager->template makeBuffer<VertexBuffer>());
        CIE_TEST_CHECK_NOTHROW(p_elementBuffer = p_bufferManager->template makeBuffer<ElementBuffer>());
    }

    {
        CIE_TEST_CASE_INIT("bind buffers")

        CIE_TEST_CHECK_NOTHROW(p_bufferManager->bindBuffer(p_vertexBuffer));
        CIE_TEST_CHECK_NOTHROW(p_bufferManager->bindBuffer(p_elementBuffer));
    }
}


} // namespace cie::gl