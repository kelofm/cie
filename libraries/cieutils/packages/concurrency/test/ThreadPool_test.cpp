// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/testing/inc/essentials.hpp"
#include "packages/concurrency/inc/ThreadPool.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <vector>
#include <thread>
#include <chrono>
#include <algorithm>


namespace cie::mp {


CIE_TEST_CASE("ThreadPool", "[concurrency]")
{
    CIE_TEST_CASE_INIT("ThreadPool")

    {
        ThreadPoolBase threads;
        ThreadPool<> pool(threads);
        std::vector<std::thread::id> threadIDs(2 * pool.size());

        for (Size index=0; index<threadIDs.size(); ++index) {
            auto job = [index, &threadIDs]() -> void {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                threadIDs[index] = std::this_thread::get_id();
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            };

            pool.queueJob(job);
        }
        pool.barrier();

        for (auto threadID : threadIDs) {
            CIE_TEST_CHECK(std::count(threadIDs.begin(), threadIDs.end(), threadID) == 2);
        }
    }

    {
        ThreadPoolBase threads;
        auto pool = ThreadPool<>(threads).firstPrivate<int>(10);
        const auto& r_poolStorage = pool.getStorage();
        CIE_TEST_REQUIRE(0 < r_poolStorage.size());
        for (const auto& r_storage : r_poolStorage) {
            CIE_TEST_CHECK(r_storage.get<0>() == 10);
        }
    }

    {
        const Size maxThreads = ThreadPoolBase::maxNumberOfThreads();
        constexpr std::size_t threadCount = 2;

        if (threadCount <= maxThreads) {
            CIE_TEST_CASE_INIT("parallel execution")
            using Clock = std::chrono::high_resolution_clock;
            using TimePoint = Clock::time_point;
            using Duration = std::chrono::milliseconds;
            constexpr Duration jobDuration = Duration(100);

            StaticArray<std::pair<TimePoint,TimePoint>,threadCount> ranges;
            {
                ThreadPoolBase threads(threadCount);
                CIE_TEST_REQUIRE(threads.size() == threadCount);
                for (unsigned i=0; i<threadCount; ++i) {
                    threads.queueJob([&ranges,i, jobDuration=jobDuration](){
                        ranges[i].first = Clock::now();
                        std::this_thread::sleep_for(jobDuration);
                        ranges[i].second = Clock::now();
                    });
                } // for i in range(2)
            }

            StaticArray<Duration,threadCount> durations;
            std::transform(ranges.begin(),
                        ranges.end(),
                        durations.begin(),
                        [](const auto pair) {
                                return std::chrono::duration_cast<Duration>(pair.second - pair.first);
                        });

            for (auto duration : durations) {
                CIE_TEST_CHECK(jobDuration <= duration);
            }
            const Duration total = std::accumulate(durations.begin(),
                                                   durations.end(),
                                                   Duration(0));

            const TimePoint begin = std::min_element(ranges.begin(),
                                                     ranges.end(),
                                                     [](auto left, auto right){
                                                         return left.first < right.first;
                                                     })->first;
            const TimePoint end = std::max_element(ranges.begin(),
                                                   ranges.end(),
                                                   [](auto left, auto right){
                                                       return left.second < right.second;
                                                   })->second;
            CIE_TEST_CHECK(std::chrono::duration_cast<Duration>(end - begin) < total);
        }
    }
}


} // namespace cie::mp