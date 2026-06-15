#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/io/inc/StatusStream.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class LoggedOperator : public LinearOperator<TSpace> {
public:
    constexpr LoggedOperator()
        : LoggedOperator(std::make_shared<StatusStream<typename TSpace::Value>>())
    {}

    explicit LoggedOperator(Ref<const std::shared_ptr<StatusStream<typename TSpace::Value>>> rpLogger)
        : _pLogger(rpLogger)
    {}

    [[nodiscard]] std::shared_ptr<StatusStream<typename TSpace::Value>> streamLogger() const {
        CIE_CHECK_POINTER(_pLogger)
        return _pLogger;
    }

    virtual void setLogger(Ref<const std::shared_ptr<StatusStream<typename TSpace::Value>>> rpLogger) {
        _pLogger = rpLogger;
    }

private:
    std::shared_ptr<StatusStream<typename TSpace::Value>> _pLogger;
}; // class LoggedOperator


template <LinearOperatorLike TOperator>
class LoggedOperatorWrapper
    :   public TOperator,
        public LoggedOperator<typename TOperator::LinalgSpace> {
public:
    template <class ...TArgs>
    LoggedOperatorWrapper(TArgs&& ... rArgs)
        : TOperator(std::forward<TArgs>(rArgs)...)
    {}

    using typename TOperator::LinalgSpace;

    /// @copydoc LinearOperator::product
    void product(
        typename LinalgSpace::Value inScale,
        typename LinalgSpace::ConstVectorView in,
        typename LinalgSpace::Value outScale,
        typename LinalgSpace::VectorView out) override {
            TOperator::product(inScale, in, outScale, out);
    }
}; // class LoggedOperatorWrapper


template <LinearOperatorLike TOperator, class ...TArgs>
std::shared_ptr<LoggedOperatorWrapper<TOperator>>
makeLoggedOperator(TArgs&&... rArgs) {
    return std::make_shared<LoggedOperatorWrapper<TOperator>>(
        std::forward<TArgs>(rArgs)...);
}


} // namespace cie::linalg
