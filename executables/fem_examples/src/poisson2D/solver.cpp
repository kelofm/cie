// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/src/SparseCholesky/SimplicialCholesky.h"

// --- Internal Includes ---
#include "poisson2D/solver.hpp"

// --- Linalg Includes ---
#include "packages/utilities/inc/reorder.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/CSROperator.hpp"
#include "packages/solvers/inc/MaskedCSROperator.hpp"
#include "packages/solvers/inc/DiagonalOperator.hpp"
#include "packages/solvers/inc/JacobiOperator.hpp"
#include "packages/solvers/inc/ConjugateGradients.hpp"
#include "packages/solvers/inc/MaskedJacobiOperator.hpp"
#include "packages/solvers/inc/MaskedIdentityOperator.hpp"

// --- Utility Includes ---
#include "packages/io/inc/MatrixMarket.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"


namespace cie::fem {


void solveEigenCG(
    CSRWrapper lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs) {
        using EigenSparseMatrix = Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>;
        Eigen::Map<EigenSparseMatrix> lhsAdaptor(
            lhs.rowCount,
            lhs.columnCount,
            lhs.entries.size(),
            lhs.rowExtents.data(),
            lhs.columnIndices.data(),
            lhs.entries.data());
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
    CSRWrapper lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs) {
        using EigenSparseMatrix = Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>;
        Eigen::Map<EigenSparseMatrix> lhsAdaptor(
            lhs.rowCount,
            lhs.columnCount,
            lhs.entries.size(),
            lhs.rowExtents.data(),
            lhs.columnIndices.data(),
            lhs.entries.data());
        Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> rhsAdaptor(rhs.data(), rhs.size(), 1);
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(solution.data(), solution.size(), 1);
        Eigen::SimplicialLLT<EigenSparseMatrix> solver;
        solver.compute(lhsAdaptor);
        solutionAdaptor = solver.solve(rhsAdaptor);
} // solveEigenLLT


void solveCG(
    CSRWrapper lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<mp::ThreadPoolBase> rThreads) {
        using LinalgSpace = linalg::DefaultSpace<Scalar,tags::SMP>;
        auto pSpace = std::make_shared<LinalgSpace>(rThreads);
        auto pLinearOperator = std::make_shared<linalg::CSROperator<int,Scalar>>(
            lhs.columnCount,
            lhs.rowExtents,
            lhs.columnIndices,
            lhs.entries,
            rThreads);
        std::shared_ptr<linalg::LinearOperator<LinalgSpace>> pPreconditioner;
        pPreconditioner = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
            linalg::makeDiagonalOperator<Scalar,int,Scalar>(
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
                pSpace));
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


void solveMultigridd(
    CSRWrapper lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Assembler> rAssembler,
    Ref<mp::ThreadPoolBase> rThreads) {
        using LinalgSpace = linalg::DefaultSpace<Scalar,tags::SMP>;
        auto pSpace = std::make_shared<LinalgSpace>(rThreads);

        std::vector<Scalar> ansatzMask(rAssembler.dofCount());
        makeAnsatzMask<Dimension, Scalar>(
            rAssembler,
            polynomialOrder + 1,
            ansatzMask);
        std::optional<Scalar> initialResidual;

        for (int iOrder=1; iOrder<polynomialOrder + 1; ++iOrder) {
            auto pLinearOperator = std::make_shared<linalg::MaskedCSROperator<int,Scalar,Scalar,Scalar>>(
                lhs.columnCount,
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
                ansatzMask,
                iOrder + 1,
                rThreads);
            std::shared_ptr<linalg::LinearOperator<LinalgSpace>> pPreconditioner;
            pPreconditioner = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
                linalg::makeDiagonalOperator<Scalar,int,Scalar>(
                    lhs.rowExtents,
                    lhs.columnIndices,
                    lhs.entries,
                    pSpace));
            linalg::ConjugateGradients<LinalgSpace>::Statistics settings {
                .iterationCount = static_cast<std::size_t>(1e3),
                .absoluteResidual = 1e-6,
                .relativeResidual = 1e-6};

            if (initialResidual.has_value()) {
                settings.absoluteResidual = 1e-6 * initialResidual.value();
            }

            linalg::ConjugateGradients<LinalgSpace> solver(
                pLinearOperator,
                pSpace,
                pPreconditioner,
                settings,
                /*verbosity=*/1);

            auto gridResidual = pSpace->makeVector(pSpace->size(rhs));
            for (std::size_t i=0; i<rhs.size(); ++i) gridResidual[i] = ansatzMask[i] < iOrder + 1 ? rhs[i] : 0.0;
            solver.product(0, gridResidual, 1, solution);
            const auto stats = solver.getStats().value();

            if (!initialResidual.has_value()) {
                initialResidual = stats.absoluteResidual / stats.relativeResidual;
            }

            std::cout << std::format(
                "{:>8} iterations {:>10.4E} abs res {:>10.4E} relative res\n",
                stats.iterationCount,
                stats.absoluteResidual,
                stats.relativeResidual);
        }
}


void solveMultigrid(
    CSRWrapper lhs,
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

        auto pLhs = std::make_shared<linalg::CSROperator<int,Scalar>>(
                lhs.columnCount,
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
                rThreads);
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

        // Lowest grid level is a proper linear solver.
        {
            const Scalar threshold = 2; // <== order + 1
            auto pGridLhs = std::make_shared<linalg::MaskedCSROperator<int,Scalar,Scalar,Scalar>>(
                lhs.columnCount,
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
                ansatzMask,
                threshold,
                rThreads);
            auto pPreconditioner = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>(
                linalg::makeDiagonalOperator<Scalar,int,Scalar>(
                    lhs.rowExtents,
                    lhs.columnIndices,
                    lhs.entries,
                    pSpace));
            linalg::ConjugateGradients<LinalgSpace>::Statistics settings {
                .iterationCount = static_cast<std::size_t>(1e3),
                .absoluteResidual = 1e-6,
                .relativeResidual = 5e-1};
            auto pOperator = std::make_shared<linalg::ConjugateGradients<LinalgSpace>>(
                pGridLhs,
                pSpace,
                pPreconditioner,
                settings,
                /*verbosity=*/1);
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
                lhs.columnCount,
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
                ansatzMask,
                threshold,
                rThreads);
            auto pOperator = std::make_shared<linalg::MaskedJacobiOperator<int,Scalar,Scalar,Scalar>>(
                lhs.columnCount,
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
                /*iterations=*/1,
                /*relaxation=*/7e-1,
                ansatzMask,
                threshold,
                pSpace);
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

            for (auto itGrid=grids.begin(); itGrid!=grids.end(); ++itGrid) {
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
    CSRWrapper lhs,
    std::span<Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<mp::ThreadPoolBase> rThreads) {
        using LinalgSpace = linalg::DefaultSpace<Scalar,tags::SMP>;
        auto pSpace = std::make_shared<LinalgSpace>(rThreads);
        auto pLinearOperator = std::make_shared<linalg::CSROperator<int,Scalar>>(
            lhs.columnCount,
            lhs.rowExtents,
            lhs.columnIndices,
            lhs.entries,
            rThreads);

        const std::size_t iterations = 1e3;
        std::shared_ptr<linalg::LinearOperator<LinalgSpace>> pSmoother;
        pSmoother = std::make_shared<linalg::JacobiOperator<int,Scalar>>(
            lhs.columnCount,
            lhs.rowExtents,
            lhs.columnIndices,
            lhs.entries,
            /*iterations=*/iterations,
            /*relaxation=*/9e-1,
            pSpace);
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
    CSRWrapper lhs,
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
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
                reorderingStrategy,
                rThreads);
            reorder<int,Scalar>(
                reordering,
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
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
                    lhs.rowCount,
                    lhs.columnCount,
                    lhs.entries.size(),
                    lhs.rowExtents.data(),
                    lhs.columnIndices.data(),
                    lhs.entries.data());
            }

            {
                std::ofstream file("rhs.mm");
                cie::io::MatrixMarket::Output io(file);
                io(rhs.data(), rhs.size());
            }
        }

        const std::string solver = rArguments.get<std::string>("solver");
        if (solver == "cg") {
            solveCG(lhs, solution, rhs, rThreads);
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
