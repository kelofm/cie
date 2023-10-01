// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/concurrency/inc/MPI.hpp"
#include "packages/concurrency/inc/MPISingleton.hpp"
#include "cmake_variables.hpp"


namespace cie::mpi {


CIE_TEST_CASE("MPI", "[concurrency]")
{
    CIE_TEST_CASE_INIT("MPI")

    Ref<MPI> mpi = MPISingleton::get();
    Size numberOfRanks = mpi.size();

    for (int iteration=0; iteration<1e6; ++iteration) {
        if (mpi.isMaster()) {
            for (int id=1; id<static_cast<int>(numberOfRanks); ++id)
                mpi.send((id+1)*(id+1), RankID(id));
        } else {
            int message;
            mpi.receive(message, mpi.getMasterRankID());
            const int id = static_cast<int>(mpi.getRankID());
            CIE_TEST_CHECK(message == (id+1) * (id+1));
        }
    }
}


} // namespace cie::mpi
