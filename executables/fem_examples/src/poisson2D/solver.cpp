// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/src/SparseCholesky/SimplicialCholesky.h"

// --- Internal Includes ---
#include "poisson2D/solver.hpp"

// --- Linalg Includes ---
#include "hipSYCL/sycl/device_selector.hpp"
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

// --- Utility Includes ---
#include "packages/io/inc/MatrixMarket.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"


namespace cie::fem {


void solveEigenCG(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs) {
        using EigenSparseMatrix = Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>;
        Eigen::Map<EigenSparseMatrix> lhsAdaptor(
            lhs.rowCount(),
            lhs.columnCount(),
            lhs.entries().size(),
            lhs.rowExtents().data(),
            lhs.columnIndices().data(),
            lhs.entries().data());
        Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> rhsAdaptor(rhs.data(), rhs.size(), 1);
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(solution.data(), solution.size(), 1);
        Eigen::ConjugateGradient<
            EigenSparseMatrix,
            Eigen::Lower | Eigen::Upper,
            Eigen::DiagonalPreconditioner<Scalar>
        > solver;
        solver.setMaxIterations(int(1e3));
        solver.setTolerance(1e-6);

        solver.compute(lhsAdaptor);
        solutionAdaptor = solver.solve(rhsAdaptor);

        std::cout << solver.iterations() << " iterations "
                  << solver.error()      << " residual\n";
} // solveEigenCG


void solveEigenLLT(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs) {
        using EigenSparseMatrix = Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>;
        Eigen::Map<EigenSparseMatrix> lhsAdaptor(
            lhs.rowCount(),
            lhs.columnCount(),
            lhs.entries().size(),
            lhs.rowExtents().data(),
            lhs.columnIndices().data(),
            lhs.entries().data());
        Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> rhsAdaptor(rhs.data(), rhs.size(), 1);
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(solution.data(), solution.size(), 1);
        Eigen::SimplicialLLT<EigenSparseMatrix> solver;
        solver.compute(lhsAdaptor);
        solutionAdaptor = solver.solve(rhsAdaptor);
} // solveEigenLLT


void solveCG(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<mp::ThreadPoolBase> rThreads) {
        using LinalgSpace = linalg::DefaultSpace<Scalar,tags::SMP>;
        auto pSpace = std::make_shared<LinalgSpace>(rThreads);
        auto pLinearOperator = std::make_shared<linalg::CSROperator<int,Scalar>>(lhs, rThreads);
        std::shared_ptr<linalg::LinearOperator<LinalgSpace>> pPreconditioner;
        pPreconditioner = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
            linalg::makeDiagonalOperator<Scalar,int,Scalar>(lhs, pSpace));
        linalg::ConjugateGradients<LinalgSpace>::Statistics settings {
            .iterationCount = static_cast<std::size_t>(1e3),
            .absoluteResidual = 1e-6,
            .relativeResidual = 1e-6};
        linalg::ConjugateGradients<LinalgSpace> solver(
            pLinearOperator,
            pSpace,
            pPreconditioner,
            settings,
            3);
        solver.product(0, rhs, 1, solution);
        const auto stats = solver.getStats().value();
        std::cout
            << stats.iterationCount << " iterations "
            << stats.relativeResidual << " residual\n";
}


void solveSYCLCG(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs) {
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
            .iterationCount = static_cast<std::size_t>(1e3),
            .absoluteResidual = 1e-6,
            .relativeResidual = 1e-6};
        linalg::ConjugateGradients<LinalgSpace> solver(
            pLinearOperator,
            pSpace,
            pPreconditioner,
            settings,
            3);

        // Solve the system.
        solver.product(0, deviceRHS, 1, deviceSolution);

        // Fetch the solution.
        pSpace->assign(solution, deviceSolution);
}


void solveMultigrid(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Assembler> rAssembler,
    Ref<mp::ThreadPoolBase> rThreads) {
        using LinalgSpace = linalg::DefaultSpace<Scalar,tags::SMP>;
        using Operator = linalg::LinearOperator<LinalgSpace>;
        auto pSpace = std::make_shared<LinalgSpace>(rThreads);

        std::vector<Scalar> ansatzMask(rAssembler.dofCount());
        makeAnsatzMask<Dimension, Scalar>(
            rAssembler,
            polynomialOrder + 1,
            ansatzMask);

        auto pLhs = std::make_shared<linalg::CSROperator<int,Scalar>>(lhs, rThreads);
        auto residual = pSpace->makeVector(pSpace->size(rhs));

        // Compute the initial residual.
        pSpace->assign(residual, rhs);
        const Scalar initialResidualNorm = std::sqrt(pSpace->innerProduct(residual, residual));
        std::cout << std::format("initial abs {:.4E}\n", initialResidualNorm);

        // Construct grids.
        struct Grid {
            //std::vector<int> rowExtents;
            //std::vector<int> columnIndices;
            //std::vector<Scalar> entries;
            std::shared_ptr<Operator> pOperator;
            std::shared_ptr<Operator> pRestriction;
            std::shared_ptr<Operator> pLhs;};
        std::vector<Grid> grids;

        auto pInverseDiagonal = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
            linalg::makeDiagonalOperator<Scalar,int,Scalar>(lhs, pSpace));

        // Lowest grid level is a proper linear solver.
        {
            const Scalar threshold = 2; // <== order + 1
            auto pGridLhs = std::make_shared<linalg::MaskedCSROperator<int,Scalar,Scalar,Scalar>>(
                lhs,
                ansatzMask,
                threshold,
                rThreads);
            linalg::ConjugateGradients<LinalgSpace>::Statistics settings {
                .iterationCount = static_cast<std::size_t>(1e3),
                .absoluteResidual = 0,
                .relativeResidual = 5e-1};
            auto pOperator = std::make_shared<linalg::ConjugateGradients<LinalgSpace>>(
                pGridLhs,
                pSpace,
                pInverseDiagonal,
                settings,
                /*verbosity=*/3);
            auto pRestriction = std::make_shared<linalg::MaskedIdentityOperator<LinalgSpace>>(
                pSpace,
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
            auto pGridLhs = std::make_shared<linalg::MaskedCSROperator<int,Scalar,Scalar,Scalar>>(
                lhs,
                ansatzMask,
                threshold,
                rThreads);
            auto pMask = std::make_shared<linalg::MaskedIdentityOperator<LinalgSpace>>(
                pSpace,
                ansatzMask,
                threshold);
            auto pSmoother = std::make_shared<linalg::JacobiOperator<LinalgSpace>>(
                pSpace,
                pSpace->size(solution),
                pGridLhs,
                pInverseDiagonal,
                /*iterations=*/6,
                /*relaxation=*/2.0 / 3.0);
            auto pOperator = std::make_shared<linalg::NestedProductOperator<LinalgSpace>>(
                pSpace,
                pSmoother,
                pMask,
                pSpace->size(solution));
            auto pRestriction = std::make_shared<linalg::MaskedIdentityOperator<LinalgSpace>>(
                pSpace,
                ansatzMask,
                threshold);
            grids.push_back(Grid {
                .pOperator = pOperator,
                .pRestriction = pRestriction,
                .pLhs = pGridLhs});
        } // for iOrder in range(2, polynomialOrder + 1)

        Scalar residualNorm = initialResidualNorm;
        auto gridResidual = pSpace->makeVector(pSpace->size(rhs));
        auto solutionUpdate = pSpace->makeVector(pSpace->size(solution));

        //grids = std::vector<Grid>(1, grids.back());
        while (1e-6 * initialResidualNorm < residualNorm) {
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
            std::cout << std::format("abs {:>10.4E} rel {:>10.4E}\n", residualNorm, residualNorm / initialResidualNorm);
        } // while not converged
}


void solveJacobi(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<mp::ThreadPoolBase> rThreads) {
        using LinalgSpace = linalg::DefaultSpace<Scalar,tags::SMP>;
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
            /*iterations=*/iterations,
            /*relaxation=*/2.0 / 3.0);
        pSmoother->product(0, rhs, 1, solution);

        // Compute residual.
        auto residual = pSpace->makeVector(pSpace->size(solution));
        pSpace->assign(residual, rhs);
        pLinearOperator->product(1, solution, -1, residual);
        const auto residualNorm = pSpace->innerProduct(residual, residual);
        const auto initialResidualNorm = pSpace->innerProduct(rhs, rhs);

        std::cout << std::format(
            "{} iterations {:.4E} residual\n",
            iterations, residualNorm / initialResidualNorm);
}


void solve(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<Scalar> rhs,
    Ref<Assembler> rAssembler,
    Ref<mp::ThreadPoolBase> rThreads,
    Ref<const utils::ArgParse::Results> rArguments) {

        ReorderingStrategy reorderingStrategy = ReorderingStrategy::None;
        Ref<const std::string> reorderingName = rArguments.get<std::string>("reordering");
        if (reorderingName == "cuthill-mckee") reorderingStrategy = ReorderingStrategy::CuthillMcKee;
        else if (reorderingName == "reverse-cuthill-mckee") reorderingStrategy = ReorderingStrategy::ReverseCuthillMcKee;
        else if (reorderingName != "none") CIE_THROW(Exception, "unknown reordering strategy: " << reorderingName)

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
            std::vector<int> buffer(reordering.size());
            reverseReorder<int>(reordering, buffer, rThreads);
            rAssembler.reorder<int>(reordering);
        }

        if (rArguments.get<bool>("write-linear-system")) {
            {
                std::ofstream file("lhs.mm");
                cie::io::MatrixMarket::Output io(file);
                io(
                    lhs.rowCount(),
                    lhs.columnCount(),
                    lhs.entries().size(),
                    lhs.rowExtents().data(),
                    lhs.columnIndices().data(),
                    lhs.entries().data());
            }

            {
                std::ofstream file("rhs.mm");
                cie::io::MatrixMarket::Output io(file);
                io(rhs.data(), rhs.size());
            }
        }

        const std::string solver = rArguments.get<std::string>("solver");
        if (solver == "cg") {
            //solveCG(lhs, solution, rhs, rThreads);
            solveSYCLCG(lhs, solution, rhs);
        } else if (solver == "p-multigrid") {
            solveMultigrid(lhs, solution, rhs, rAssembler, rThreads);
        } else if (solver == "jacobi") {
            solveJacobi(lhs, solution, rhs, rThreads);
        } else if (solver == "eigen-cg") {
            solveEigenCG(lhs, solution, rhs);
        } else if (solver == "eigen-llt") {
            solveEigenLLT(lhs, solution, rhs);
        } else CIE_THROW(Exception, std::format("unknown solver \"{}\"", solver))

        if (rArguments.get<bool>("write-linear-system")) {
            std::ofstream file("solution.mm");
            cie::io::MatrixMarket::Output io(file);
            io(solution.data(), solution.size());
        }
}


} // namespace cie::fem
