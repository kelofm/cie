// --- Internal Includes ---
#include "poisson2D/postprocessing.hpp"

// --- FEM Includes ---
#include "packages/io/inc/VTKHDF.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"
#include "packages/maths/inc/OuterProduct.hpp"

// --- STL Includes ---
#include <cstdint>


namespace cie::fem {


void postprocess(
    std::span<const Scalar> meshBase,
    std::span<const Scalar> meshLengths,
    CSRWrapper lhs,
    std::span<const Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Mesh> rMesh,
    std::span<const CellData> contiguousCellData,
    Ref<const BVH> rBVH,
    Ref<const Assembler> rAssembler,
    Ref<const utils::ArgParse::Results> rArguments,
    Ref<mp::ThreadPoolBase> rThreads) {
        CIE_BEGIN_EXCEPTION_TRACING
            const unsigned postprocessResolution = rArguments.get<std::size_t>("scatter-resolution");

            std::vector<Scalar> residual(rhs.begin(), rhs.end());
            {
                Eigen::Map<const Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>> lhsAdaptor(
                    lhs.rowCount,
                    lhs.columnCount,
                    lhs.entries.size(),
                    lhs.rowExtents.data(),
                    lhs.columnIndices.data(),
                    lhs.entries.data());
                Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(
                    solution.data(),
                    solution.size());
                Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> residualAdaptor(
                    residual.data(),
                    residual.size());
                residualAdaptor.noalias() -= lhsAdaptor * solutionAdaptor;
            }

            {
                cie::io::VTKHDF::Output io;

                // Write the mesh.
                io(
                    "mesh",
                    rMesh);

                // Write solution, load and residual.
                {
                    io.writeFieldVariables<Scalar>(
                        "mesh",
                        rMesh,
                        rAssembler,
                        {
                            {"state", solution},
                            {"load", rhs},
                            {"residual", residual}
                        });
                }

                // Write cell and DoF IDs.
                {
                    std::vector<std::size_t> cellIDs, dofIDs;
                    cellIDs.reserve(contiguousCellData.size());
                    dofIDs.reserve(contiguousCellData.size() * rMesh.data().ansatz(0).size());

                    for (Ref<const CellData> rCell : contiguousCellData) {
                        cellIDs.push_back(rCell.id());
                        const auto& rDoFIDs = rAssembler[rCell.id()];
                        std::copy(
                            rDoFIDs.begin(),
                            rDoFIDs.end(),
                            std::back_inserter(dofIDs));
                    } // for rCell in contiguousCellData

                    io.writeCellVariables<std::size_t>(
                        "mesh",
                        {
                            {"id", {1}, cellIDs},
                            {"dofs", {rMesh.data().ansatz(0).size()}, dofIDs}
                        });
                }

                // Pointwise sampling.
                {
                    auto logBlock = utils::LoggerSingleton::get().newBlock("scatter postprocess");
                    const Scalar postprocessDelta  = 1.0 / (postprocessResolution - 1);

                    const std::size_t sampleCount = intPow(postprocessResolution, 2);

                    std::vector<unsigned> cellIDs(sampleCount);
                    std::vector<Scalar> sampleCoordinates(sampleCount * Dimension);
                    std::vector<std::vector<Scalar>>
                        sampleState(polynomialOrder + 1),
                        sampleLoad(polynomialOrder + 1),
                        sampleResidual(polynomialOrder + 1);

                    std::vector<std::uint8_t> ansatzMask(intPow(polynomialOrder + 1, Dimension));
                    makeAnsatzMask<Dimension>(
                        polynomialOrder + 1,
                        ansatzMask);

                    for (std::size_t iOrder=0ul; iOrder<polynomialOrder + 1; ++iOrder) {
                        sampleState[iOrder].resize(sampleCount);
                        sampleLoad[iOrder].resize(sampleCount);
                        sampleResidual[iOrder].resize(sampleCount);
                    }

                    mp::ParallelFor<std::size_t>(rThreads).firstPrivate(std::vector<Scalar>(), std::vector<Scalar>()).execute(
                        sampleCount,
                        [&](const std::size_t iSample,
                            Ref<std::vector<Scalar>> rBuffer,
                            Ref<std::vector<Scalar>> rResults) -> void {
                                const std::size_t iSampleY = iSample / postprocessResolution;
                                const std::size_t iSampleX = iSample % postprocessResolution;
                                const auto physicalCoordinates = Kernel<Dimension,Scalar>::cast<PhysicalCoordinate<Scalar>>(std::span<Scalar,Dimension> {
                                    sampleCoordinates.data() + iSample * Dimension,
                                    Dimension});
                                physicalCoordinates[0] = meshBase[0] + meshLengths[0] * iSampleX * postprocessDelta;
                                physicalCoordinates[1] = meshBase[1] + meshLengths[1] * iSampleY * postprocessDelta;

                                // Find which cell the point lies in.
                                const auto iMaybeCellData = rBVH.makeView().find(
                                    Kernel<Dimension,Scalar>::decay(physicalCoordinates),
                                    std::span<const CellData>(contiguousCellData));

                                if (iMaybeCellData != contiguousCellData.size()) {
                                    Ref<const CellData> rCellData = contiguousCellData[iMaybeCellData];
                                    cellIDs[iSample] = rCellData.id();

                                    // Compute sample point in the cell's parametric space.
                                    StaticArray<ParametricCoordinate<Scalar>,Dimension> parametricCoordinates;
                                    rBuffer.resize(rCellData.makeSpatialTransform().bufferSize());
                                    rCellData.transform(
                                        physicalCoordinates,
                                        Kernel<Dimension,Scalar>::view(parametricCoordinates),
                                        rBuffer);

                                    // Evaluate the cell's ansatz functions at the parametric sample point.
                                    Ref<const Ansatz> rAnsatzSpace = rMesh.data().ansatz(rCellData.ansatzID());
                                    rResults.resize(rAnsatzSpace.size());
                                    rBuffer.resize(rAnsatzSpace.bufferSize());
                                    rAnsatzSpace.evaluate(
                                        Kernel<Dimension,Scalar>::decayView(
                                            parametricCoordinates),
                                        rResults,
                                        rBuffer);

                                    // Find the entries of the cell's DoFs in the global state vector.
                                    const auto& rGlobalIndices = rAssembler[rCellData.id()];

                                    for (std::size_t iOrder=0ul; iOrder<polynomialOrder; ++iOrder) {
                                        // Compute state as an indirect inner product of the solution vector
                                        // and the ansatz function values at the local corner coordinates.
                                        sampleState[iOrder][iSample] = 0.0;
                                        sampleLoad[iOrder][iSample] = 0.0;
                                        sampleResidual[iOrder][iSample] = 0.0;
                                        for (unsigned iFunction=0u; iFunction<rGlobalIndices.size(); ++iFunction) {
                                            if (ansatzMask[iFunction] == static_cast<std::uint8_t>(iOrder + 1)) {
                                                sampleState[iOrder][iSample] += solution[rGlobalIndices[iFunction]] * rResults[iFunction];
                                                sampleLoad[iOrder][iSample] += rhs[rGlobalIndices[iFunction]] * rResults[iFunction];
                                                sampleResidual[iOrder][iSample] += residual[rGlobalIndices[iFunction]] * rResults[iFunction];
                                            }
                                        } // for iFunction in range(rGlobalIndices.size())

                                        sampleState.back()[iSample] += sampleState[iOrder][iSample];
                                        sampleLoad.back()[iSample] += sampleLoad[iOrder][iSample];
                                        sampleResidual.back()[iSample] += sampleResidual[iOrder][iSample];
                                    }
                                } else {
                                    // Could not find the cell that contains the current sample point.
                                    for (std::size_t iOrder=0ul; iOrder<polynomialOrder; ++iOrder) {
                                        sampleState[iOrder][iSample] = NAN;
                                        sampleLoad[iOrder][iSample] = NAN;
                                        sampleResidual[iOrder][iSample] = NAN;
                                    }
                                }
                        });

                    io.writePointCloud<Scalar,Dimension>(
                        "samples",
                        sampleCoordinates,
                        /*gridSize=*/postprocessResolution);
                    io.writeFieldVariables<Scalar>(
                        "samples",
                        {
                            {"state", {1}, sampleState.back()},
                            {"load",  {1}, sampleLoad.back()},
                            {"residual", {1}, sampleResidual.back()}
                        });

                    for (std::size_t iOrder=0ul; iOrder<polynomialOrder; ++iOrder) {
                        const std::string stateName = std::format("state_p{}", iOrder + 1);
                        const std::string loadName = std::format("load_p{}", iOrder + 1);
                        const std::string residualName = std::format("residual_p{}", iOrder + 1);
                        io.writeFieldVariables<Scalar>(
                            "samples",
                            {
                                {stateName, {1}, sampleState[iOrder]},
                                {loadName,  {1}, sampleLoad[iOrder]},
                                {residualName, {1}, sampleResidual[iOrder]}
                            });
                    }
                }
            }
        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
