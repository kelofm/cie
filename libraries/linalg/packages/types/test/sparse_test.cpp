// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- Internal Includes ---
#include "packages/types/inc/sparse.hpp"

// --- STL Includes ---
#include <vector>

namespace cie::linalg
{

/*
 *  ------------ Test quad connectivity -------------
 *
 *     33----34-----35----36-----37-----------39
 *     |            |            |            |
 *     |            26  28   30  32           |
 *     23    24     |            |            |
 *     |            25  27   29  31           |
 *     |            |            |            |
 *     17----18-----19--20---21--22-----------38
 *     |            |            |
 *     7   9   11   13    15     |
 *     |            |            16
 *     6   8   10   12    14     |
 *     |            |            |
 *     0---1----2---3------4-----5
 *
 *
 * location maps (with randomly chosen local numbering):
 *   [[ 0, 1, 2, 3, 6, 8, 10, 12, 7, 9, 11, 13, 17, 18, 19 ]
 *    [ 3, 12, 13, 19, 4, 14, 15, 21, 5, 16, 22, 20 ]
 *    [ 17, 23, 33, 34, 35, 26, 25, 19, 18, 24 ]
 *    [ 38, 22, 39, 37, 31, 32 ]
 *    [ 27, 29, 28, 30, 37, 36, 35, 26, 32, 31, 25, 22, 19, 21, 20 ]]
 */

CIE_TEST_CASE("allocateSparseMatrix_test", "[types]")
{
    std::vector<std::vector<size_t>> locationMaps
    {
        { 0, 1, 2, 3, 6, 8, 10, 12, 7, 9, 11, 13, 17, 18, 19 },
        { 3, 12, 13, 19, 4, 14, 15, 21, 5, 16, 22, 20 },
        { 17, 23, 33, 34, 35, 26, 25, 19, 18, 24 },
        { 38, 22, 39, 37, 31, 32 },
        { 27, 29, 28, 30, 37, 36, 35, 26, 32, 31, 25, 22, 19, 21, 20 }
    };

    CompressedSparseRowMatrix<size_t> sparseMatrix;

    CIE_TEST_REQUIRE_NOTHROW(sparseMatrix.allocate(locationMaps));

    CIE_TEST_REQUIRE(sparseMatrix.rowSize() == 40);
    CIE_TEST_REQUIRE(sparseMatrix.columnSize() == 40);

    auto sparseDataStructure = sparseMatrix.release();

    size_t* indices = std::get<0>(sparseDataStructure);
    size_t* indptr = std::get<1>(sparseDataStructure);
    double* data = std::get<2>(sparseDataStructure);

    for(auto& locationMap : locationMaps) {
        std::sort(locationMap.begin(), locationMap.end());
    }

    std::vector<std::vector<size_t>> dofElementMapping = {
        { 0 }   , { 0 }   , { 0 }      , { 0, 1 }, { 1 },          // dofs  0 -  4
        { 1 }   , { 0 }   , { 0 }      , { 0 }   , { 0 },          // dofs  5 -  9
        { 0 }   , { 0 }   , { 0, 1 }   , { 0, 1 }, { 1 },          // dofs 10 - 14
        { 1 }   , { 1 }   , { 0, 2 }   , { 0, 2 }, { 0, 1, 2, 4 }, // dofs 15 - 19
        { 1, 4 }, { 1, 4 }, { 1, 3, 4 }, { 2 }   , { 2 },          // dofs 20 - 24
        { 2, 4 }, { 2, 4 }, { 4 }      , { 4 }   , { 4 },          // dofs 25 - 29
        { 4 }   , { 3, 4 }, { 3, 4 }   , { 2 }   , { 2 },          // dofs 30 - 34
        { 2, 4 }, { 4 }   , { 3, 4 }   , { 3 }   , { 3 }           // dofs 35 - 39
    };

    std::vector<size_t> connectingDofSizes;
    std::vector<size_t> expectedIndices;

    for(const std::vector<size_t>& indices : dofElementMapping)
    {
        std::vector<size_t> concatenated;

        for(const auto& elementIndex : indices)
        {
            const auto& locationMap = locationMaps[elementIndex];

            concatenated.insert(concatenated.end(), locationMap.begin(), locationMap.end());
        }

        std::sort(concatenated.begin(), concatenated.end());
        concatenated.erase(std::unique(concatenated.begin(), concatenated.end()), concatenated.end());

        connectingDofSizes.push_back(concatenated.size());
        expectedIndices.insert(expectedIndices.end(), concatenated.begin(), concatenated.end());
    };

    size_t nnz = expectedIndices.size();

    CIE_TEST_REQUIRE(indptr[40] == nnz);

    for(size_t iEntry = 0; iEntry < nnz; ++iEntry)
    {
        CIE_TEST_CHECK(indices[iEntry] == expectedIndices[iEntry]);
        CIE_TEST_CHECK(data[iEntry] == 0.0);
    }

    size_t rowPtr = 0;

    for(size_t iDof = 0; iDof < 40; ++iDof)
    {
        CIE_TEST_CHECK(indptr[iDof] == rowPtr);

        rowPtr += connectingDofSizes[iDof];
    }

    CIE_TEST_CHECK(indptr[40] == rowPtr);

    delete[] indptr;
    delete[] indices;
    delete[] data;

} // allocateSparseMatrix_test

} // namespace cie::linalg