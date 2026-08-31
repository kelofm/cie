// --- Internal Includes ---
#include "poisson2D/postprocessing.hpp"

// --- FEM Includes ---
#include "packages/io/inc/VTKHDF.hpp"
#include "packages/postprocessing/inc/AnsatzSpacePostprocessor.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"

// --- STL Includes ---
#include <cstdint>


namespace cie::fem {


void sumComponents(
    std::span<const Scalar> components,
    std::size_t componentCount,
    std::span<Scalar> sums) {
            const std::size_t itemCount = components.size() / componentCount;
            for (std::size_t iItem=0ul; iItem<itemCount; ++iItem) {
                sums[iItem] = std::accumulate(
                    components.begin() + componentCount * iItem,
                    components.begin() + componentCount * (iItem + 1),
                    Scalar(0));
            }
}


void scatter(
    Ref<cie::io::VTKHDF::Output> rIO,
    std::string_view groupName,
    std::span<const Scalar> coordinates,
    std::vector<std::pair<
        std::string,
        std::span<const Scalar>
    >> fieldCoefficients,
    Ref<const Mesh> rMesh,
    Ref<const BVH> rBVH,
    Ref<const Assembler> rAssembler,
    bool writeComponents,
    Scalar nanReplacement) {
        CIE_BEGIN_EXCEPTION_TRACING
            // Select polynomial components to collect.
            std::vector<std::uint8_t> ansatzOrders(intPow(polynomialOrder + 1, Dimension));
            std::vector<std::uint8_t> requestedOrders;
            if (writeComponents) {
                makeAnsatzMask<Dimension,std::uint8_t>(
                    polynomialOrder + 1,
                    ansatzOrders);
                requestedOrders.resize(polynomialOrder);
                std::iota(
                    requestedOrders.begin(),
                    requestedOrders.end(),
                    std::uint8_t(1));
            } else {
                std::fill(
                    ansatzOrders.begin(),
                    ansatzOrders.end(),
                    std::uint8_t(0));
                requestedOrders.push_back(0);
            }

            // Allocate memory for the output values.
            const std::size_t pointCount = coordinates.size() / Dimension;
            std::vector<std::span<const Scalar>> input;
            input.reserve(fieldCoefficients.size());
            std::transform(
                fieldCoefficients.begin(),
                fieldCoefficients.end(),
                std::back_inserter(input),
                [] (const auto& rPair) {return rPair.second;});

            std::vector<Scalar> outputBuffer(fieldCoefficients.size() * requestedOrders.size() * pointCount);
            std::vector<std::span<Scalar>> output;
            output.reserve(fieldCoefficients.size());
            for (std::size_t iField=0ul; iField<fieldCoefficients.size(); ++iField)
                output.emplace_back(
                    outputBuffer.data() + iField * requestedOrders.size() * pointCount,
                    requestedOrders.size() * pointCount);

            // Perform the interpolation.
            AnsatzSpacePostprocessor<Ansatz> postprocessor(nanReplacement);
            postprocessor.interpolate(
                Kernel<Dimension,Scalar>::cast<PhysicalCoordinate<Scalar>>(coordinates),
                rMesh,
                rAssembler,
                rBVH,
                input,
                1,
                ansatzOrders,
                requestedOrders,
                output);

            for (std::size_t iField=0ul; iField<fieldCoefficients.size(); ++iField) {
                std::string fieldName = fieldCoefficients[iField].first;
                if (writeComponents) fieldName += "-components";
                rIO.writeFieldVariables<Scalar>(
                    groupName,
                    {{fieldName, {requestedOrders.size()}, output[iField]}});
            } // for iField

            // If components were requested, they need to be summed up
            // to output the complete field as well.
            if (writeComponents) {
                for (std::size_t iField=0ul; iField<fieldCoefficients.size(); ++iField) {
                    const std::size_t itemCount = output[iField].size() / requestedOrders.size();
                    for (std::size_t iItem=0ul; iItem<itemCount; ++iItem) {
                        output[iField][iItem] = std::accumulate(
                            output[iField].begin() + requestedOrders.size() * iItem,
                            output[iField].begin() + requestedOrders.size() * (iItem + 1),
                            Scalar(0));
                    } // for iItem

                    rIO.writeFieldVariables<Scalar>(
                        groupName,
                        {{
                            fieldCoefficients[iField].first,
                            {1},
                            {output[iField].data(), itemCount}
                        }});
                } // for iField
            } // if writeComponents
        CIE_END_EXCEPTION_TRACING
}


void postprocess(
    std::span<const Scalar> meshBase,
    std::span<const Scalar> meshLengths,
    linalg::CSRView<Scalar,int> lhs,
    std::span<const Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Mesh> rMesh,
    Ref<const BVH> rBVH,
    std::span<const BoundarySegment> boundarySegments,
    Ref<const Assembler> rAssembler,
    Ref<const cie::io::JSONObject> rConfiguration,
    Ref<mp::ThreadPoolBase> rThreads) {
        CIE_BEGIN_EXCEPTION_TRACING
            // Compute the residual.
            std::vector<Scalar> residual(rhs.begin(), rhs.end());
            {
                Eigen::Map<const Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>> lhsAdaptor(
                    lhs.rowCount(),
                    lhs.columnCount(),
                    lhs.entries().size(),
                    lhs.rowExtents().data(),
                    lhs.columnIndices().data(),
                    lhs.entries().data());
                Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(
                    solution.data(),
                    solution.size());
                Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> residualAdaptor(
                    residual.data(),
                    residual.size());
                residualAdaptor.noalias() -= lhsAdaptor * solutionAdaptor;
            }

            // Label ansatz functions by their polynomial orders.
            std::vector<std::uint8_t> ansatzOrders(intPow(polynomialOrder + 1, Dimension));
            makeAnsatzMask<Dimension,std::uint8_t>(
                polynomialOrder + 1,
                ansatzOrders);

            // Define which polynomial components of the solution we request.
            // In this case, all of them.
            std::vector<std::uint8_t> requestedOrders(polynomialOrder);
            std::iota(
                requestedOrders.begin(),
                requestedOrders.end(),
                std::uint8_t(1));

            {
                cie::io::VTKHDF::Output io;
                const auto discretizationConfig = rConfiguration["discretization"]["postprocessing"];

                // Write the mesh.
                io.writeCells(
                    "cells",
                    rMesh.data().cells());

                // Write solution, load and residual.
                {
                    std::vector<std::pair<std::string,std::span<const Scalar>>> fieldCoefficients;
                    if (discretizationConfig["state"].as<bool>())
                        fieldCoefficients.emplace_back("state", solution);
                    if (discretizationConfig["load"].as<bool>())
                        fieldCoefficients.emplace_back("load", rhs);
                    if (discretizationConfig["residual"].as<bool>())
                        fieldCoefficients.emplace_back("residual", residual);

                    io.writeFieldVariables<Scalar>(
                        "cells",
                        rMesh,
                        rMesh.data().cells(),
                        rAssembler,
                        fieldCoefficients);
                }

                // Write cell IDs.
                if (discretizationConfig["cell-id"].as<bool>()) {
                    std::vector<std::size_t> buffer;
                    buffer.reserve(rMesh.data().cells().size());
                    for (Ref<const Cell> rCell : rMesh.data().cells()) {
                        buffer.push_back(rCell.id());
                    }
                    io.writeCellVariables<std::size_t>(
                        "cells",
                        {{"id", {1}, buffer}});
                }

                // Write DoF IDs.
                if (discretizationConfig["dof-id"].as<bool>()) {
                    std::vector<std::size_t> buffer;
                    buffer.reserve(rMesh.data().cells().size() * rMesh.data().ansatz(0).size());
                    for (Ref<const Cell> rCell : rMesh.data().cells()) {
                        const auto& rDoFIDs = rAssembler[rCell.id()];
                        std::copy(
                            rDoFIDs.begin(),
                            rDoFIDs.end(),
                            std::back_inserter(buffer));
                    }
                    io.writeCellVariables<std::size_t>(
                        "cells",
                        {{"DoFs", {rMesh.data().ansatz(0).size()}, buffer}});
                }

                // Write triangulated domains.
                std::size_t iTesselatedDomain = 0ul;
                for (const auto domainConfig : rConfiguration["domains"]) {
                    if (domainConfig["type"].as<std::string>() == "default") continue;
                    const auto& [rDomainID, rTriangles] = rMesh.data().domainTriangles()[iTesselatedDomain++];

                    // Write triangulation.
                    const std::string domainName = std::format("domain-{}", rDomainID);
                    io.writeSTL<Scalar,Dimension>(
                        domainName,
                        rTriangles);

                    std::vector<std::pair<std::string,std::span<const Scalar>>> fieldCoefficients;
                    if (domainConfig["postprocessing"]["state"].as<bool>())
                        fieldCoefficients.emplace_back("state", solution);
                    if (domainConfig["postprocessing"]["load"].as<bool>())
                        fieldCoefficients.emplace_back("load", rhs);
                    if (domainConfig["postprocessing"]["residual"].as<bool>())
                        fieldCoefficients.emplace_back("residual", residual);

                    scatter(
                        io,
                        domainName,
                        rTriangles,
                        fieldCoefficients,
                        rMesh,
                        rBVH,
                        rAssembler,
                        domainConfig["postprocessing"]["write-components"].as<bool>(),
                        domainConfig["postprocessing"]["replace-nans"].is<std::nullptr_t>()
                            ? std::numeric_limits<Scalar>::quiet_NaN()
                            : domainConfig["postprocessing"]["replace-nans"].as<double>());
                } // for rDomain, rTriangles in rMesh.data().domainTriangles()

                // Write dirichlet constraints.
                {
                    std::vector<Scalar> segmentPoints(boundarySegments.size() * 2 * Dimension);
                    std::vector<Scalar> segmentStates(boundarySegments.size() * 2);
                    for (std::size_t iSegment=0ul; iSegment<boundarySegments.size(); ++iSegment) {
                        std::copy_n(
                            boundarySegments[iSegment].data(),
                            Dimension,
                            segmentPoints.data() + iSegment * 2 * Dimension);
                        std::copy_n(
                            boundarySegments[iSegment].data() + Dimension,
                            Dimension,
                            segmentPoints.data() + iSegment * 2 * Dimension + Dimension);
                        std::copy_n(
                            boundarySegments[iSegment].data() + Dimension + Dimension,
                            2,
                            segmentStates.data() + iSegment * 2);
                    }
                    io.writeSegments<Scalar,Dimension>(
                        "dirichlet-1d",
                        segmentPoints);
                    io.writeFieldVariables<Scalar>(
                        "dirichlet-1d",
                        {{"prescribed-state", {1}, segmentStates}});

                    const auto constraintConfig = rConfiguration["dirichlet-1d"]["postprocessing"];
                    std::vector<std::pair<std::string,std::span<const Scalar>>> fieldCoefficients;
                    if (constraintConfig["state"].as<bool>())
                        fieldCoefficients.emplace_back("state", solution);
                    if (constraintConfig["load"].as<bool>())
                        fieldCoefficients.emplace_back("load", rhs);
                    if (constraintConfig["residual"].as<bool>())
                        fieldCoefficients.emplace_back("residual", residual);
                    scatter(
                        io,
                        "dirichlet-1d",
                        segmentPoints,
                        fieldCoefficients,
                        rMesh,
                        rBVH,
                        rAssembler,
                        constraintConfig["write-components"].as<bool>(),
                        constraintConfig["replace-nans"].is<std::nullptr_t>()
                            ? std::numeric_limits<Scalar>::quiet_NaN()
                            : constraintConfig["replace-nans"].as<double>());
                }

                // Pointwise sampling.
                {
                    auto logBlock = utils::LoggerSingleton::get().newBlock("scatter");

                    // Parse configuration.
                    const auto localConfig = rConfiguration["discretization"]["postprocessing"];
                    std::array<std::size_t,Dimension> scatterResolution;

                    {
                        const auto resolutionConfig = discretizationConfig["resolution"];
                        if (resolutionConfig.is<std::size_t>() || resolutionConfig.is<double>()) {
                            std::fill(
                                scatterResolution.begin(),
                                scatterResolution.end(),
                                resolutionConfig.as<std::size_t>());
                        } else {
                            const auto resolution = resolutionConfig.as<std::vector<std::size_t>>();
                            CIE_CHECK(
                                resolution.size() == Dimension,
                                std::format(
                                    "expecting an array of size {}, but got one with {} items",
                                    Dimension,
                                    resolution.size()))
                            std::copy(
                                resolution.begin(),
                                resolution.end(),
                                scatterResolution.begin());
                        }
                    }

                    std::array<Scalar,Dimension> postprocessDelta;
                    for (std::size_t iDimension=0ul; iDimension<Dimension; ++iDimension)
                        postprocessDelta[iDimension] = meshLengths[iDimension] / (scatterResolution[iDimension] - 1);

                    const std::size_t sampleCount = std::accumulate(
                        scatterResolution.begin(),
                        scatterResolution.end(),
                        1ul,
                        std::multiplies<std::size_t>());
                    std::vector<unsigned> cellIDs(sampleCount);
                    std::vector<Scalar> sampleCoordinates(sampleCount * Dimension);

                    mp::ParallelFor<std::size_t>(rThreads).execute(
                        sampleCount,
                        [&](const std::size_t iSample) -> void {
                                const std::size_t iSampleY = iSample / scatterResolution[0];
                                const std::size_t iSampleX = iSample - iSampleY * scatterResolution[0];
                                const auto physicalCoordinates = Kernel<Dimension,Scalar>::cast<PhysicalCoordinate<Scalar>>(std::span<Scalar,Dimension> {
                                    sampleCoordinates.data() + iSample * Dimension,
                                    Dimension});
                                physicalCoordinates[0] = meshBase[0] + iSampleX * postprocessDelta[0];
                                physicalCoordinates[1] = meshBase[1] + iSampleY * postprocessDelta[1];

                                // Find which cell the point lies in.
                                const auto iMaybeCellData = rBVH.makeView().find(
                                    Kernel<Dimension,Scalar>::decay(physicalCoordinates),
                                    std::span<const Cell>(rMesh.data().cells()));

                                if (iMaybeCellData != rMesh.data().cells().size()) {
                                    Ref<const Cell> rCell = rMesh.data().cells()[iMaybeCellData];
                                    cellIDs[iSample] = rCell.id();
                                }
                        });

                    io.writePointCloud<Scalar,Dimension>(
                        "scatter",
                        sampleCoordinates,
                        scatterResolution);
                    std::vector<std::pair<std::string,std::span<const Scalar>>> fieldCoefficients;
                    if (discretizationConfig["state"].as<bool>())
                        fieldCoefficients.emplace_back("state", solution);
                    if (discretizationConfig["load"].as<bool>())
                        fieldCoefficients.emplace_back("load", rhs);
                    if (discretizationConfig["residual"].as<bool>())
                        fieldCoefficients.emplace_back("residual", residual);
                    scatter(
                        io,
                        "scatter",
                        sampleCoordinates,
                        fieldCoefficients,
                        rMesh,
                        rBVH,
                        rAssembler,
                        discretizationConfig["write-components"].as<bool>(),
                        discretizationConfig["replace-nans"].is<std::nullptr_t>()
                            ? std::numeric_limits<Scalar>::quiet_NaN()
                            : discretizationConfig["replace-nans"].as<double>());
                }
            }
        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
