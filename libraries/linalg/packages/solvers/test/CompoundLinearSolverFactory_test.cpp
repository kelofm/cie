// --- Utility Includes ---
#include "packages/io/inc/json.hpp"
#include "packages/testing/inc/essentials.hpp"

// --- Linalg Includes ---
#include "packages/solvers/inc/CompoundLinearSolverFactory.hpp"


namespace cie::linalg {


CIE_TEST_CASE("CompoundLinearSolverFactory", "[solvers]") {
    CIE_TEST_CASE_INIT("CompoundLinearSolverFactory")
    using Scalar1 = double;
    using Index1  = int;
    using Scalar2 = float;
    using Index2  = int;

    using ScalarSpace1 = DefaultSpace<Scalar1>;
    using IndexSpace1  = DefaultSpace<Index1>;
    using Factory1 = LinearSolverFactory<ScalarSpace1,IndexSpace1>;

    using ScalarSpace2 = DefaultSpace<Scalar2>;
    using IndexSpace2  = DefaultSpace<Index2>;
    using Factory2 = LinearSolverFactory<ScalarSpace2,IndexSpace2>;

    using Factory = CompoundLinearSolverFactory<Factory1,Factory2>;

    const auto pScalarSpace1 = std::make_shared<ScalarSpace1>();
    const auto pIndexSpace1 = std::make_shared<IndexSpace1>();

    const auto pScalarSpace2 = std::make_shared<ScalarSpace2>();
    const auto pIndexSpace2 = std::make_shared<IndexSpace2>();

    Ref<Factory1> rFactory1 = Factory1::Singleton::get();
    Ref<Factory2> rFactory2 = Factory2::Singleton::get();
    Factory factory(rFactory1, rFactory2);

    const std::vector<Index1> rowExtents {0, 1}, columnIndices {0};
    const std::vector<Scalar1> entries {1.0};
    CSRView<const Scalar1,const Index1> matrix(
        1, rowExtents, columnIndices, entries);
    cie::io::JSONObject configuration (std::string {R"({"type" : "cg"})"});
    CIE_TEST_CHECK_THROWS(factory.make(
        configuration,
        {pScalarSpace1, pScalarSpace2},
        {pIndexSpace1, pIndexSpace2},
        matrix));
    CIE_TEST_CHECK_NOTHROW(factory.make(
        configuration,
        {pScalarSpace1, pScalarSpace2},
        {pIndexSpace1, pIndexSpace2},
        matrix,
        /*allowAmbiguous=*/true));
    const auto pMaybeSolver = factory.make(
        configuration,
        {pScalarSpace1, pScalarSpace2},
        {pIndexSpace1, pIndexSpace2},
        matrix,
        /*allowAmbiguous=*/true);

    CIE_TEST_REQUIRE(pMaybeSolver.has_value());
} // CIE_TEST_CASE("CompoundLinearSolverFactory")


} // namespace cie::linalg
