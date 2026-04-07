// --- Linalg Includes ---
#include "packages/solvers/inc/IterativeSolver.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"

// --- Utility Includes ---
#include "packages/types/inc/Color.hpp"

// --- STL Includes ---
#include <format>
#include <sstream>


namespace cie::linalg {


template <LinalgSpaceLike TS>
Ref<std::string> IterativeSolver<TS>::makeIterationReport(
    Ref<std::string> rOutput,
    int verbosity,
    ReportType reportType,
    Ref<const Statistics> rStatus,
    Ref<const Statistics> rSettings) const {
        rOutput.clear();
        std::stringstream message;
        const bool converged = rStatus.absoluteResidual < rSettings.absoluteResidual || rStatus.relativeResidual < rSettings.relativeResidual;

        if (2 <= verbosity) {
            if (rStatus.iterationCount == 0) {
                message << "+ --------- + ----------------- + ----------------- +\n"
                        << "| iteration | absolute residual | relative residual |\n"
                        << "+ --------- + ----------------- + ----------------- +";
            }

            if (reportType == ReportType::Termination || converged || (3 <= verbosity && reportType == ReportType::Iteration)) {
                if (rStatus.iterationCount == 0) {
                    message << "\n";
                }
                message << std::format(
                    "| {}{:>9}{} | {}{:>17.5E}{} | {}{:>17.5E}{} |",
                    converged ? RGBAColor::TUMGreen.ANSI() : RGBAColor::TUMOrange.ANSI(),
                    rStatus.iterationCount,
                    RGBAColor::White.ANSI(),
                    rStatus.absoluteResidual < rSettings.absoluteResidual ? RGBAColor::TUMGreen.ANSI() : RGBAColor::TUMOrange.ANSI(),
                    rStatus.absoluteResidual,
                    RGBAColor::White.ANSI(),
                    rStatus.relativeResidual < rSettings.relativeResidual ? RGBAColor::TUMGreen.ANSI() : RGBAColor::TUMOrange.ANSI(),
                    rStatus.relativeResidual,
                    RGBAColor::White.ANSI());
            }

            if (converged)
                message << "\n+ --------- + ----------------- + ----------------- +";
        }

        rOutput = message.str();
        return rOutput;
}


#define CIE_INSTANTIATE_ITERATIVE_SOLVER(T)                         \
    template class IterativeSolver<DefaultSpace<T,tags::Serial>>;   \
    template class IterativeSolver<DefaultSpace<T,tags::SMP>>;


CIE_INSTANTIATE_ITERATIVE_SOLVER(float)
CIE_INSTANTIATE_ITERATIVE_SOLVER(double)


#undef CIE_INSTANTIATE_ITERATIVE_SOLVER


#ifdef CIE_ENABLE_SYCL

template class IterativeSolver<SYCLSpace<float>>;
template class IterativeSolver<SYCLSpace<double>>;

#endif


} // namespace cie::linalg
