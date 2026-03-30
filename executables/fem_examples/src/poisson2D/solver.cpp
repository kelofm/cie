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
#include "packages/solvers/inc/DiagonalOperator.hpp"
#include "packages/solvers/inc/ConjugateGradients.hpp"

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
        auto pLinalgSpace = std::make_shared<LinalgSpace>(rThreads);
        auto pLinearOperator = std::make_shared<linalg::CSROperator<int,Scalar>>(
            lhs.columnCount,
            lhs.rowExtents,
            lhs.columnIndices,
            lhs.entries,
            rThreads);
        auto pPreconditioner = std::make_shared<linalg::DiagonalOperator<LinalgSpace>>();
        *pPreconditioner = linalg::makeDiagonalOperator<Scalar,int,Scalar>(
            lhs.rowExtents,
            lhs.columnIndices,
            lhs.entries,
            pLinalgSpace);
        linalg::ConjugateGradients<LinalgSpace>::Statistics settings {
            .iterationCount = static_cast<std::size_t>(1e3),
            .absoluteResidual = 1e-6,
            .relativeResidual = 1e-6};
        linalg::ConjugateGradients<LinalgSpace> solver(
            pLinearOperator,
            pLinalgSpace,
            pPreconditioner,
            settings,
            3);
        solver.product(rhs, 1, solution);
        const auto stats = solver.getStats().value();
        std::cout
            << stats.iterationCount << " iterations "
            << stats.relativeResidual << " residual\n";
}


void solve(
    CSRWrapper lhs,
    std::span<Scalar> solution,
    std::span<Scalar> rhs,
    [[maybe_unused]] Ref<const Assembler> rAssembler,
    Ref<mp::ThreadPoolBase> rThreads,
    Ref<const utils::ArgParse::Results> rArguments) {
        std::vector<int> reordering(solution.size());
        ReorderingStrategy reorderingStrategy = ReorderingStrategy::None;
        Ref<const std::string> reorderingName = rArguments.get<std::string>("reordering");
        if (reorderingName == "cuthill-mckee") reorderingStrategy = ReorderingStrategy::CuthillMcKee;
        else if (reorderingName == "reverse-cuthill-mckee") reorderingStrategy = ReorderingStrategy::ReverseCuthillMcKee;
        else if (reorderingName != "none") CIE_THROW(Exception, "unknown reordering strategy: " << reorderingName)

        if (reorderingStrategy != ReorderingStrategy::None) {
            auto logBlock = utils::LoggerSingleton::get().newBlock("reorder");
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
        }

        if (rArguments.get<bool>("write-lhs")) {
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

        const std::string solver = rArguments.get<std::string>("solver");
        if (solver == "cg") {
            solveCG(lhs, solution, rhs, rThreads);
        } else if (solver == "eigen-cg") {
            solveEigenCG(lhs, solution, rhs);
        } else if (solver == "eigen-llt") {
            solveEigenLLT(lhs, solution, rhs);
        } else CIE_THROW(Exception, std::format("unknown solver \"{}\"", solver))

        if (reorderingStrategy != ReorderingStrategy::None) {
            auto logBlock = utils::LoggerSingleton::get().newBlock("reverse reorder");
            reverseReorder<int,Scalar>(
                reordering,
                lhs.rowExtents,
                lhs.columnIndices,
                lhs.entries,
                rThreads);
            reverseReorder<int,Scalar>(
                reordering,
                rhs,
                rThreads);
            reverseReorder<int,Scalar>(
                reordering,
                solution,
                rThreads);
        }
}


} // namespace cie::fem
