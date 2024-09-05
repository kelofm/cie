// --- Utility Includes ---
#include "packages/io/inc/MatrixMarket.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <iostream> // std::cin, std::cout


namespace cie::utils::io {


struct MatrixMarket::Output::Impl {
    Ptr<std::ostream> _pStream;

    MatrixMarket::Settings _settings;
}; // struct MatrixMarket::Output::Impl


MatrixMarket::Output::Output()
    : Output(std::cout, Settings {})
{
}


MatrixMarket::Output::Output(Ref<std::ostream> rStream,
                             Settings settings)
    : _pImpl(new Impl {._pStream = &rStream,
                       ._settings = settings})
{
}


MatrixMarket::Output::~Output() = default;


#define CIE_MATRIX_MARKET_OUTPUT_DEFINITION(TIndex, TValue)                                                         \
    Ref<MatrixMarket::Output> MatrixMarket::Output::operator()(const TIndex rowCount,                               \
                                                               const TIndex columnCount,                            \
                                                               const TIndex nonzeroCount,                           \
                                                               Ptr<const TIndex> pRowExtents,                       \
                                                               Ptr<const TIndex> pColumnIndices,                    \
                                                               Ptr<const TValue> pNonzeros)                         \
    {                                                                                                               \
        CIE_BEGIN_EXCEPTION_TRACING                                                                                 \
                                                                                                                    \
        Ref<std::ostream> rStream = *_pImpl->_pStream;                                                              \
                                                                                                                    \
        rStream << "%%MatrixMarket matrix coordinate real general\n";                                               \
        rStream << rowCount << ' ' << columnCount << ' ' << nonzeroCount << '\n';                                   \
        for (TIndex iRow=0; iRow<rowCount; ++iRow) {                                                                \
            const TIndex iBegin = pRowExtents[iRow];                                                                \
            const TIndex iEnd = pRowExtents[iRow + 1];                                                              \
            for (TIndex iNonzero=iBegin; iNonzero<iEnd; ++iNonzero) {                                               \
                rStream << iRow + 1 << ' ' << pColumnIndices[iNonzero] + 1 << ' ' << pNonzeros[iNonzero] << '\n';   \
            } /* for iNonzero in range(iBegin, iEnd) */                                                             \
        } /* for iRow in range(rowCount) */                                                                         \
                                                                                                                    \
        return *this;                                                                                               \
        CIE_END_EXCEPTION_TRACING                                                                                   \
    }

#define CIE_MATRIX_MARKET_OUTPUT_DEFINITION_FOR_VALUE(TValue)                                                       \
    Ref<MatrixMarket::Output> MatrixMarket::Output::operator()(Ptr<const TValue> itBegin,                           \
                                                               const std::size_t size)                              \
    {                                                                                                               \
        CIE_BEGIN_EXCEPTION_TRACING                                                                                 \
        Ref<std::ostream> rStream = *_pImpl->_pStream;                                                              \
        rStream << "%%MatrixMarket matrix array real general\n"                                                     \
                << size << ' ' << 1 << '\n';                                                                        \
        Ptr<const TValue> itEnd = itBegin + size;                                                                   \
        for (; itBegin!=itEnd; ++itBegin) rStream << *itBegin << '\n';                                              \
        return *this;                                                                                               \
        CIE_END_EXCEPTION_TRACING                                                                                   \
    }                                                                                                               \
                                                                                                                    \
    CIE_MATRIX_MARKET_OUTPUT_DEFINITION(int, TValue)                                                                \
    CIE_MATRIX_MARKET_OUTPUT_DEFINITION(unsigned, TValue)                                                           \
    CIE_MATRIX_MARKET_OUTPUT_DEFINITION(std::size_t, TValue)


CIE_MATRIX_MARKET_OUTPUT_DEFINITION_FOR_VALUE(double)
#undef CIE_MATRIX_MARKET_OUTPUT_DEFINITION
#undef CIE_MATRIX_MARKET_OUTPUT_DEFINITION_FOR_VALUE


} // namespace cie::utils::io

