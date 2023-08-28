// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/concurrency/inc/ThreadPoolSingleton.hpp"

// --- STL Includes ---
#include <vector>
#include <thread>
#include <chrono>
#include <algorithm>


namespace cie::mp {


CIE_TEST_CASE( "ThreadPool", "[concurrency]" )
{
    CIE_TEST_CASE_INIT( "ThreadPool" )

    {
        auto pool = ThreadPoolSingleton::get();
        std::vector<std::thread::id> threadIDs( 2 * pool.size() );

        for ( Size index=0; index<threadIDs.size(); ++index )
        {
            auto job = [index, &threadIDs]() -> void
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                threadIDs[index] = std::this_thread::get_id();
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            };

            pool.queueJob( job );
        }
        pool.barrier();

        for ( auto threadID : threadIDs )
            CIE_TEST_CHECK( std::count( threadIDs.begin(), threadIDs.end(), threadID ) == 2 );
    }

    {
        auto pool = ThreadPoolSingleton::firstPrivate<int>(10);
        const auto& r_poolStorage = pool.getStorage();
        CIE_TEST_REQUIRE(0 < r_poolStorage.size());
        for (const auto& r_storage : r_poolStorage)
            CIE_TEST_CHECK(r_storage.get<0>() == 10);
    }
}


} // namespace cie::mp