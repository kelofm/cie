// --- LinAlg Includes ---
#include "packages/utilities/inc/linalghelper.hpp"


namespace cie::linalg::linalghelper {


void write(const DynamicEigenMatrix<double>& matrix, std::ostream& out)
{
    Size rowSize = matrix.rowSize();
    Size columnSize = matrix.columnSize();

    for(Size i = 0; i < rowSize; ++i)
        linalghelper::writeRow([&](Size j){return matrix(i, j);},
                               columnSize,
                               out,
                               12);
}


} // namespace cie::linalg::linalghelper
