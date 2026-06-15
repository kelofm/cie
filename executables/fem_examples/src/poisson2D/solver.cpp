// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/src/SparseCholesky/SimplicialCholesky.h"

// --- Internal Includes ---
#include "poisson2D/solver.hpp"

// --- Linalg Includes ---
#include "packages/constrained_solvers/inc/ConstrainedSolverFactory.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/solvers/inc/CompoundLinearSolverFactory.hpp"
#include "packages/utilities/inc/reorder.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/CSROperator.hpp"
#include "packages/solvers/inc/MaskedCSROperator.hpp"
#include "packages/solvers/inc/DiagonalOperator.hpp"
#include "packages/solvers/inc/JacobiOperator.hpp"
#include "packages/solvers/inc/ConjugateGradients.hpp"
#include "packages/solvers/inc/NestedProductOperator.hpp"
#include "packages/solvers/inc/MaskedIdentityOperator.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/solvers/inc/SYCLCSROperator.hpp"
#include "packages/solvers/inc/SYCLMaskedCSROperator.hpp"
#include "packages/utilities/inc/CSRUtility.hpp"

// --- Utility Includes ---
#include "packages/io/inc/MatrixMarket.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"


namespace cie::fem {


void solveMultigrid(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Assembler> rAssembler,
    Ref<mp::ThreadPoolBase> rThreads,
    Ref<const cie::io::JSONObject> rConfiguration) {
        using LinalgSpace = linalg::DefaultSpace<Scalar>;
        using Operator = linalg::LinearOperator<LinalgSpace>;
        auto pSpace = std::make_shared<LinalgSpace>(rThreads);

        using MaskScalar = std::uint16_t;
        using MaskSpace = linalg::DefaultSpace<MaskScalar>;
        auto pMaskSpace = std::make_shared<MaskSpace>();

        std::vector<MaskScalar> ansatzMask(rAssembler.dofCount());
        makeAnsatzMask<Dimension,MaskScalar>(
            rAssembler,
            polynomialOrder + 1,
            ansatzMask);

        auto pLhs = std::make_shared<linalg::CSROperator<int,Scalar>>(lhs, rThreads);
        auto residual = pSpace->makeVector(pSpace->size(rhs));

        // Compute the initial residual.
        pSpace->assign(residual, rhs);
        const Scalar initialResidualNorm = std::sqrt(pSpace->innerProduct(residual, residual));

        // Construct grids.
        struct Grid {
            std::shared_ptr<Operator> pOperator;
            std::shared_ptr<Operator> pRestriction;
            std::shared_ptr<Operator> pLhs;};
        std::vector<Grid> grids;

        auto pInverseDiagonal = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
            linalg::makeDiagonalOperator<Scalar,int,Scalar>(lhs, pSpace));

        // Lowest grid level is a proper linear solver.
        {
            const Scalar threshold = 2; // <== order + 1
            auto pGridLhs = std::make_shared<linalg::MaskedCSROperator<int,Scalar,Scalar,MaskScalar>>(
                lhs,
                ansatzMask,
                threshold,
                rThreads);
            linalg::ConjugateGradients<LinalgSpace>::Status settings {
                .iterationCount = rConfiguration["solver"]["max-iterations"].as<std::size_t>(),
                .absoluteResidual = rConfiguration["solver"]["absolute-residual"].as<double>(),
                .relativeResidual = rConfiguration["solver"]["relative-residual"].as<double>()};
            auto pOperator = std::make_shared<linalg::ConjugateGradients<LinalgSpace>>(
                pGridLhs,
                pSpace,
                pInverseDiagonal,
                settings,
                linalg::Verbosity::Warnings);
            auto pRestriction = std::make_shared<linalg::MaskedIdentityOperator<LinalgSpace,MaskSpace>>(
                pSpace,
                pMaskSpace,
                ansatzMask,
                threshold);
            grids.push_back(Grid {
                .pOperator = pOperator,
                .pRestriction = pRestriction,
                .pLhs = pGridLhs});
        }

        // The rest of the grids are jacobi smoothers.
        for (std::size_t iOrder=2ul; iOrder<polynomialOrder+1; ++iOrder) {
            const Scalar threshold = iOrder + 1;
            auto pGridLhs = std::make_shared<linalg::MaskedCSROperator<int,Scalar,Scalar,MaskScalar>>(
                lhs,
                ansatzMask,
                threshold,
                rThreads);
            auto pMask = std::make_shared<linalg::MaskedIdentityOperator<LinalgSpace,MaskSpace>>(
                pSpace,
                pMaskSpace,
                ansatzMask,
                threshold);
            auto pSmoother = std::make_shared<linalg::JacobiOperator<LinalgSpace>>(
                pSpace,
                pSpace->size(solution),
                pGridLhs,
                pInverseDiagonal,
                rConfiguration["smoother"]["iterations"].as<std::size_t>(),
                rConfiguration["smoother"]["relaxation"].as<double>());
            auto pOperator = std::make_shared<linalg::NestedProductOperator<LinalgSpace>>(
                pSpace,
                pSmoother,
                pMask,
                pSpace->size(solution));
            grids.push_back(Grid {
                .pOperator = pOperator,
                .pRestriction = pMask,
                .pLhs = pGridLhs});
        } // for iOrder in range(2, polynomialOrder + 1)

        Scalar residualNorm = initialResidualNorm;
        auto gridResidual = pSpace->makeVector(pSpace->size(rhs));
        auto solutionUpdate = pSpace->makeVector(pSpace->size(solution));

        const Scalar targetAbsoluteResidualNorm = rConfiguration["absolute-residual"].as<double>();
        const Scalar targetRelativeResidualNorm = rConfiguration["relative-residual"].as<double>();
        const std::size_t maxIterations = rConfiguration["max-iterations"].as<std::size_t>();

        for (
            std::size_t iIteration = 0ul;
            iIteration<maxIterations
                && targetAbsoluteResidualNorm <= residualNorm
                && targetRelativeResidualNorm <= (residualNorm / initialResidualNorm);
            ++iIteration) {
                for (auto itGrid=grids.rbegin(); itGrid!=grids.rend(); ++itGrid) {
                    itGrid->pRestriction->product(0, residual, 1, gridResidual);
                    pSpace->fill(solutionUpdate, 0);
                    itGrid->pOperator->product(0, gridResidual, 1, solutionUpdate);
                    pSpace->add(solution, solutionUpdate, 1);
                    itGrid->pLhs->product(1, solutionUpdate, -1, residual);
                } // for itOperator

                for (auto itGrid=grids.begin()+1; itGrid!=grids.end(); ++itGrid) {
                    itGrid->pRestriction->product(0, residual, 1, gridResidual);
                    pSpace->fill(solutionUpdate, 0);
                    itGrid->pOperator->product(0, gridResidual, 1, solutionUpdate);
                    pSpace->add(solution, solutionUpdate, 1);
                    itGrid->pLhs->product(1, solutionUpdate, -1, residual);
                } // for itOperator

                residualNorm = std::sqrt(pSpace->innerProduct(residual, residual));

                if (2 <= rConfiguration["verbosity"].as<int>())
                    std::cout << std::format("abs {:>10.4E} rel {:>10.4E}\n", residualNorm, residualNorm / initialResidualNorm);
        } // while not converged
}


#ifdef CIE_ENABLE_SYCL


void solveSYCLCG(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const cie::io::JSONObject> rConfiguration) {
        using LinalgSpace = linalg::SYCLSpace<Scalar>;
        auto pSpace = std::make_shared<LinalgSpace>(std::make_shared<sycl::queue>(sycl::default_selector_v));
        auto pIndexSpace = std::make_shared<linalg::SYCLSpace<int>>(pSpace->getQueue());

        // Copy the matrix to the device.
        auto rowExtents = pIndexSpace->makeVector(lhs.rowExtents().size());
        auto columnIndices = pIndexSpace->makeVector(lhs.columnIndices().size());
        auto entries = pSpace->makeVector(lhs.entries().size());

        pIndexSpace->assign(rowExtents, lhs.rowExtents());
        pIndexSpace->assign(columnIndices, lhs.columnIndices());
        pSpace->assign(entries, lhs.entries());

        linalg::CSRView<const Scalar,const int> deviceLHS(
            lhs.columnCount(),
            {rowExtents.get(), rowExtents.size()},
            {columnIndices.get(), columnIndices.size()},
            {entries.get(), entries.size()});

        // Copy the rhs and solution vectors to the device.
        auto deviceRHS = pSpace->makeVector(rhs.size());
        auto deviceSolution = pSpace->makeVector(solution.size());

        pSpace->assign(deviceRHS, rhs);
        pSpace->assign(deviceSolution, solution);

        // Build operators.
        auto pLinearOperator = std::make_shared<linalg::SYCLCSROperator<int,Scalar>>(deviceLHS, pSpace);
        auto pPreconditioner = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
            linalg::makeDiagonalOperator<Scalar,int,Scalar>(deviceLHS, pSpace));
        linalg::ConjugateGradients<LinalgSpace>::Statistics settings {
            .iterationCount = rConfiguration["max-iterations"].as<std::size_t>(),
            .absoluteResidual = rConfiguration["absolute-residual"].as<double>(),
            .relativeResidual = rConfiguration["relative-residual"].as<double>()};
        linalg::ConjugateGradients<LinalgSpace> solver(
            pLinearOperator,
            pSpace,
            pPreconditioner,
            settings,
            static_cast<linalg::Verbosity>(rConfiguration["verbosity"].as<int>()));

        // Solve the system.
        solver.product(0, deviceRHS, 1, deviceSolution);

        // Fetch the solution.
        pSpace->assign(solution, deviceSolution);
}


void solveSYCLMultigrid(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Assembler> rAssembler,
    Ref<const cie::io::JSONObject> rConfiguration) {
        using LinalgSpace = linalg::SYCLSpace<Scalar>;
        using Operator = linalg::LinearOperator<LinalgSpace>;

        auto pQueue = std::make_shared<sycl::queue>(sycl::default_selector_v);
        auto pSpace = std::make_shared<LinalgSpace>(pQueue);

        using MaskScalar = std::uint16_t;
        using MaskSpace = linalg::SYCLSpace<MaskScalar>;
        auto pMaskSpace = std::make_shared<MaskSpace>(pQueue);

        std::vector<MaskScalar> mask(rAssembler.dofCount());
        makeAnsatzMask<Dimension,MaskScalar>(
            rAssembler,
            polynomialOrder + 1,
            mask);

        // Copy the LHS matrix to the device.
        auto pIndexSpace = std::make_shared<linalg::SYCLSpace<int>>(pQueue);
        auto deviceRowExtents = pIndexSpace->makeVector(lhs.rowExtents().size());
        pIndexSpace->assign(deviceRowExtents, lhs.rowExtents());

        auto deviceColumnIndices = pIndexSpace->makeVector(lhs.columnIndices().size());
        pIndexSpace->assign(deviceColumnIndices, lhs.columnIndices());

        auto deviceEntries = pSpace->makeVector(lhs.entries().size());
        pSpace->assign(deviceEntries, lhs.entries());

        linalg::CSRView<const Scalar,const int> deviceLHS(
            lhs.columnCount(),
            {deviceRowExtents.get(), deviceRowExtents.size()},
            {deviceColumnIndices.get(), deviceColumnIndices.size()},
            {deviceEntries.get(), deviceEntries.size()});

        // Define arrays on the device.
        auto deviceResidual = pSpace->makeVector(rhs.size());
        pSpace->assign(deviceResidual, rhs);

        auto deviceSolution = pSpace->makeVector(solution.size());
        pSpace->assign(deviceSolution, solution);

        auto deviceMask = pMaskSpace->makeVector(mask.size());
        pMaskSpace->assign(deviceMask, mask);

        auto pLhs = std::make_shared<linalg::SYCLCSROperator<int,Scalar>>(lhs, pSpace);

        // Compute the initial residual.
        const Scalar initialResidualNorm = std::sqrt(pSpace->innerProduct(deviceResidual, deviceResidual));

        // Construct grids.
        struct Grid {
            std::shared_ptr<Operator> pOperator;
            std::shared_ptr<Operator> pRestriction;
            std::shared_ptr<Operator> pLhs;};
        std::vector<Grid> grids;

        auto pInverseDiagonal = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
            linalg::makeDiagonalOperator<Scalar,int,Scalar>(deviceLHS, pSpace));

        // Lowest grid level is a proper linear solver.
        {
            const MaskScalar threshold = 2; // <== order + 1
            auto pGridLhs = std::make_shared<linalg::SYCLMaskedCSROperator<int,Scalar,MaskScalar>>(
                deviceLHS,
                deviceMask,
                threshold,
                pSpace,
                pMaskSpace);
            linalg::ConjugateGradients<LinalgSpace>::Statistics settings {
                .iterationCount = rConfiguration["solver"]["max-iterations"].as<std::size_t>(),
                .absoluteResidual = rConfiguration["solver"]["absolute-residual"].as<double>(),
                .relativeResidual = rConfiguration["solver"]["relative-residual"].as<double>()};
            auto pOperator = std::make_shared<linalg::ConjugateGradients<LinalgSpace>>(
                pGridLhs,
                pSpace,
                pInverseDiagonal,
                settings,
                static_cast<linalg::Verbosity>(rConfiguration["solver"]["verbosity"].as<int>()));
            auto pRestriction = std::make_shared<linalg::MaskedIdentityOperator<LinalgSpace,MaskSpace>>(
                pSpace,
                pMaskSpace,
                deviceMask,
                threshold);
            grids.push_back(Grid {
                .pOperator = pOperator,
                .pRestriction = pRestriction,
                .pLhs = pGridLhs});
        }

        // The rest of the grids are jacobi smoothers.
        for (std::size_t iOrder=2ul; iOrder<polynomialOrder+1; ++iOrder) {
            const MaskScalar threshold = iOrder + 1;
            auto pGridLhs = std::make_shared<linalg::SYCLMaskedCSROperator<int,Scalar,MaskScalar>>(
                deviceLHS,
                deviceMask,
                threshold,
                pSpace,
                pMaskSpace);
            auto pMask = std::make_shared<linalg::MaskedIdentityOperator<LinalgSpace,MaskSpace>>(
                pSpace,
                pMaskSpace,
                deviceMask,
                threshold);
            auto pSmoother = std::make_shared<linalg::JacobiOperator<LinalgSpace>>(
                pSpace,
                pSpace->size(deviceSolution),
                pGridLhs,
                pInverseDiagonal,
                rConfiguration["smoother"]["iterations"].as<std::size_t>(),
                rConfiguration["smoother"]["relaxation"].as<double>());
            auto pOperator = std::make_shared<linalg::NestedProductOperator<LinalgSpace>>(
                pSpace,
                pSmoother,
                pMask,
                pSpace->size(deviceSolution));
            grids.push_back(Grid {
                .pOperator = pOperator,
                .pRestriction = pMask,
                .pLhs = pGridLhs});
        } // for iOrder in range(2, polynomialOrder + 1)

        Scalar residualNorm = initialResidualNorm;
        auto gridResidual = pSpace->makeVector(rhs.size());
        auto solutionUpdate = pSpace->makeVector(solution.size());

        const Scalar targetAbsoluteResidualNorm = rConfiguration["absolute-residual"].as<double>();
        const Scalar targetRelativeResidualNorm = rConfiguration["relative-residual"].as<double>();
        const std::size_t maxIterations = rConfiguration["max-iterations"].as<std::size_t>();

        for (
            std::size_t iIteration = 0ul;
            iIteration<maxIterations
                && targetAbsoluteResidualNorm <= residualNorm
                && targetRelativeResidualNorm <= (residualNorm / initialResidualNorm);
            ++iIteration) {
                for (auto itGrid=grids.rbegin(); itGrid!=grids.rend(); ++itGrid) {
                    itGrid->pRestriction->product(0, deviceResidual, 1, gridResidual);
                    pSpace->fill(solutionUpdate, 0);
                    itGrid->pOperator->product(0, gridResidual, 1, solutionUpdate);
                    pSpace->add(deviceSolution, solutionUpdate, 1);
                    itGrid->pLhs->product(1, solutionUpdate, -1, deviceResidual);
                } // for itOperator

                for (auto itGrid=grids.begin()+1; itGrid!=grids.end(); ++itGrid) {
                    itGrid->pRestriction->product(0, deviceResidual, 1, gridResidual);
                    pSpace->fill(solutionUpdate, 0);
                    itGrid->pOperator->product(0, gridResidual, 1, solutionUpdate);
                    pSpace->add(deviceSolution, solutionUpdate, 1);
                    itGrid->pLhs->product(1, solutionUpdate, -1, deviceResidual);
                } // for itOperator

                residualNorm = std::sqrt(pSpace->innerProduct(deviceResidual, deviceResidual));
                pSpace->assign(solution, deviceSolution);

                if (3 <= rConfiguration["verbosity"].as<int>())
                    std::cout << std::format("abs {:>10.4E} rel {:>10.4E}\n", residualNorm, residualNorm / initialResidualNorm);
        } // while not converged
}


#endif


void solveJacobi(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<mp::ThreadPoolBase> rThreads,
    Ref<const cie::io::JSONObject> rConfiguration) {
        using LinalgSpace = linalg::DefaultSpace<Scalar>;
        auto pSpace = std::make_shared<LinalgSpace>(rThreads);
        auto pLinearOperator = std::make_shared<linalg::CSROperator<int,Scalar>>(lhs, rThreads);
        auto pInverseDiagonal = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
            linalg::makeDiagonalOperator<Scalar,int,Scalar>(lhs, pSpace));

        const std::size_t iterations = 1e3;
        std::shared_ptr<linalg::LinearOperator<LinalgSpace>> pSmoother;
        pSmoother = std::make_shared<linalg::JacobiOperator<LinalgSpace>>(
            pSpace,
            pSpace->size(solution),
            pLinearOperator,
            pInverseDiagonal,
            rConfiguration["iterations"].as<std::size_t>(),
            rConfiguration["relaxation"].as<double>());
        pSmoother->product(0, rhs, 1, solution);

        // Compute residual.
        auto residual = pSpace->makeVector(pSpace->size(solution));
        pSpace->assign(residual, rhs);
        pLinearOperator->product(1, solution, -1, residual);
        const auto residualNorm = pSpace->innerProduct(residual, residual);
        const auto initialResidualNorm = pSpace->innerProduct(rhs, rhs);

        if (2 <= rConfiguration["verbosity"].as<int>())
            std::cout << std::format(
                "{} iterations {:.4E} residual\n",
                iterations, residualNorm / initialResidualNorm);
}


void solve(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<Scalar> rhs,
    linalg::CSRView<Scalar, int> constraintGradients,
    std::span<Scalar> constraintGaps,
    Ref<Assembler> rAssembler,
    Ref<mp::ThreadPoolBase> rThreads,
    Ref<const cie::io::JSONObject> rConfiguration) {
        CIE_BEGIN_EXCEPTION_TRACING
            // Parse and perform reordering if necessary.
            ReorderingStrategy reorderingStrategy = ReorderingStrategy::None;
            const auto reorderingConfiguration = rConfiguration["reordering"];
            if (!reorderingConfiguration.is<std::nullptr_t>()) {
                const std::string reorderingName = reorderingConfiguration.as<std::string>();
                if (reorderingName == "cuthill-mckee") reorderingStrategy = ReorderingStrategy::CuthillMcKee;
                else if (reorderingName == "reverse-cuthill-mckee") reorderingStrategy = ReorderingStrategy::ReverseCuthillMcKee;
                else if (reorderingName != "none") CIE_THROW(Exception, "unknown reordering strategy: " << reorderingName)
            }

            if (reorderingStrategy != ReorderingStrategy::None) {
                auto logBlock = utils::LoggerSingleton::get().newBlock("reorder");
                std::vector<int> reordering(solution.size());
                makeReordering<int,Scalar>(
                    reordering,
                    lhs.rowExtents(),
                    lhs.columnIndices(),
                    lhs.entries(),
                    reorderingStrategy,
                    rThreads);

                // Reorder unconstrained system.
                reorder<int,Scalar>(
                    reordering,
                    lhs.rowExtents(),
                    lhs.columnIndices(),
                    lhs.entries(),
                    rThreads);
                reorder<int,Scalar>(
                    reordering,
                    rhs,
                    rThreads);

                // Reorder constraints.
                reorder<int,Scalar>(
                    reordering,
                    constraintGradients.rowExtents(),
                    constraintGradients.columnIndices(),
                    constraintGradients.entries(),
                    rThreads);
                reorder<int,Scalar>(
                    reordering,
                    constraintGaps,
                    rThreads);

                // Reorder the assembler.
                std::vector<int> buffer(reordering.size());
                reverseReorder<int>(reordering, buffer, rThreads);
                rAssembler.reorder<int>(reordering);
            }

            // Output system components if requested.
            if (!rConfiguration["write-unconstrained-lhs"].is<std::nullptr_t>()) {
                std::ofstream file(rConfiguration["write-unconstrained-lhs"].as<std::string>());
                linalg::io::MatrixMarket::Output io(file);
                io(lhs);
            }

            if (!rConfiguration["write-unconstrained-rhs"].is<std::nullptr_t>()) {
                std::ofstream file(rConfiguration["write-unconstrained-rhs"].as<std::string>());
                linalg::io::MatrixMarket::Output io(file);
                io(rhs);
            }

            linalg::LinearSolverFactory<
                linalg::DefaultSpace<Scalar>,
                linalg::DefaultSpace<int>> linearSolverFactory;

            CIE_BEGIN_EXCEPTION_TRACING
                linearSolverFactory.load();
            CIE_END_EXCEPTION_TRACING

            linalg::ConstrainedSolverFactory<
                linalg::DefaultSpace<Scalar>,
                linalg::DefaultSpace<int>> constrainedSolverFactory;

            CIE_BEGIN_EXCEPTION_TRACING
                constrainedSolverFactory.load();
            CIE_END_EXCEPTION_TRACING

            CIE_BEGIN_EXCEPTION_TRACING
                const auto pMaybeConstrainedOperator = constrainedSolverFactory.make(
                    rConfiguration["solver"],
                    std::make_shared<linalg::DefaultSpace<Scalar>>(rThreads),
                    std::make_shared<linalg::DefaultSpace<int>>(rThreads),
                    constraintGradients,
                    constraintGaps,
                    lhs,
                    linearSolverFactory);
                CIE_CHECK(
                    pMaybeConstrainedOperator.has_value(),
                    std::format(
                        "\"{}\" does not name a registered constrained linear solver",
                        rConfiguration["solver"]["type"].as<std::string>()))
                pMaybeConstrainedOperator.value()->product(0, rhs, 1, solution);
            CIE_END_EXCEPTION_TRACING

            if (!rConfiguration["write-unconstrained-solution"].is<std::nullptr_t>()) {
                std::ofstream file(rConfiguration["write-unconstrained-solution"].as<std::string>());
                linalg::io::MatrixMarket::Output io(file);
                io(solution.data(), solution.size());
            }
        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
