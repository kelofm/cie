#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinalgSpace.hpp"
#include "packages/solvers/inc/LoggedOperator.hpp"
#include "packages/io/inc/StatusStream.hpp"
#include "packages/utilities/inc/CSRView.hpp"


namespace cie::linalg {


template <class TScalar>
class ConstrainedStatusStream : public StatusStream<TScalar> {
public:
    using typename StatusStream<TScalar>::Status;

    ConstrainedStatusStream();

    ConstrainedStatusStream(
        Status constraintConfiguration,
        Status linearSystemConfiguration,
        OptionalRef<std::ostream> rMaybeStream = {},
        utils::Comparison<TScalar> scalarComparison = {});

    ConstrainedStatusStream(ConstrainedStatusStream&&) noexcept;

    ~ConstrainedStatusStream() override;

    ConstrainedStatusStream& operator=(ConstrainedStatusStream&&) noexcept;

    void report(
        StatusReportType reportType,
        Verbosity verbosity) override;

    [[nodiscard]] Ref<const Status> constraintConfiguration() const;

    [[nodiscard]] Ref<const Status> constraintStatus() const;

    void configureConstraints(Status configuration);

    void submitConstraintStatus(Status status);

private:
    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class ConstrainedStatusStream


template <LinalgSpaceLike TScalarSpace, LinalgSpaceLike TIndexSpace>
class ConstrainedOperator : public LoggedOperator<TScalarSpace> {
public:
    using ScalarSpace = TScalarSpace;

    using IndexSpace = TIndexSpace;

    using Scalar = typename ScalarSpace::Value;

    using Index = typename IndexSpace::Value;

    ConstrainedOperator() noexcept = default;

    ConstrainedOperator(
        std::shared_ptr<ScalarSpace> pScalarSpace,
        std::shared_ptr<IndexSpace> pIndexSpace,
        CSRView<const Scalar,const Index> constraintGradients,
        typename ScalarSpace::ConstVectorView constraintGaps);

protected:
    std::shared_ptr<const ScalarSpace> scalarSpace() const {
        return _pScalarSpace;
    }

    std::shared_ptr<ScalarSpace> scalarSpace() {
        return _pScalarSpace;
    }

    std::shared_ptr<const IndexSpace> indexSpace() const {
        return _pIndexSpace;
    }

    std::shared_ptr<IndexSpace> indexSpace() {
        return _pIndexSpace;
    }

    CSRView<const Scalar,const Index> constraintGradients() const noexcept {
        return _constraintGradients;
    }

    typename TScalarSpace::ConstVectorView constraintGaps() const noexcept {
        return _constraintGaps;
    }

private:
    std::shared_ptr<ScalarSpace> _pScalarSpace;

    std::shared_ptr<IndexSpace> _pIndexSpace;

    CSRView<const Scalar,const Index> _constraintGradients;

    typename TScalarSpace::ConstVectorView _constraintGaps;
}; // class ConstrainedOperator


} // namespace cie::linalg
