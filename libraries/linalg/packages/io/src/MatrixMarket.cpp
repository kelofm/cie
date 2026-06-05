// --- LinAlg Includes ---
#include "packages/io/inc/MatrixMarket.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <iostream> // std::cin, std::cout


namespace cie::linalg::io {


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


#define CIE_MATRIX_MARKET_OUTPUT_DEFINITION(I, T)                                                                   \
    Ref<MatrixMarket::Output> MatrixMarket::Output::operator()(                                                     \
        const I rowCount,                                                                                           \
        const I columnCount,                                                                                        \
        const I nonzeroCount,                                                                                       \
        Ptr<const I> pRowExtents,                                                                                   \
        Ptr<const I> pColumnIndices,                                                                                \
        Ptr<const T> pNonzeros) {                                                                                   \
        CIE_BEGIN_EXCEPTION_TRACING                                                                                 \
                                                                                                                    \
        Ref<std::ostream> rStream = *_pImpl->_pStream;                                                              \
                                                                                                                    \
        rStream << "%%MatrixMarket matrix coordinate real general\n";                                               \
        rStream << rowCount << ' ' << columnCount << ' ' << nonzeroCount << '\n';                                   \
        for (I iRow=0; iRow<rowCount; ++iRow) {                                                                     \
            const I iBegin = pRowExtents[iRow];                                                                     \
            const I iEnd = pRowExtents[iRow + 1];                                                                   \
            for (I iNonzero=iBegin; iNonzero<iEnd; ++iNonzero) {                                                    \
                rStream << iRow + 1 << ' ' << pColumnIndices[iNonzero] + 1 << ' ' << pNonzeros[iNonzero] << '\n';   \
            } /* for iNonzero in range(iBegin, iEnd) */                                                             \
        } /* for iRow in range(rowCount) */                                                                         \
                                                                                                                    \
        return *this;                                                                                               \
        CIE_END_EXCEPTION_TRACING                                                                                   \
    }                                                                                                               \
    Ref<MatrixMarket::Output> MatrixMarket::Output::operator()(CSRView<const T,const I> matrix) {                   \
        return this->operator()(                                                                                    \
            matrix.rowCount(), matrix.columnCount(), matrix.entries().size(),                                       \
            matrix.rowExtents().data(), matrix.columnIndices().data(), matrix.entries().data());                    \
    }

#define CIE_MATRIX_MARKET_OUTPUT_DEFINITION_FOR_VALUE(T)                    \
    Ref<MatrixMarket::Output> MatrixMarket::Output::operator()(             \
        Ptr<const T> itBegin,                                               \
        const std::size_t size) {                                           \
            CIE_BEGIN_EXCEPTION_TRACING                                     \
            Ref<std::ostream> rStream = *_pImpl->_pStream;                  \
            rStream << "%%MatrixMarket matrix array real general\n"         \
                    << size << ' ' << 1 << '\n';                            \
            Ptr<const T> itEnd = itBegin + size;                            \
            for (; itBegin!=itEnd; ++itBegin) rStream << *itBegin << '\n';  \
            return *this;                                                   \
            CIE_END_EXCEPTION_TRACING                                       \
    }                                                                       \
    Ref<MatrixMarket::Output> MatrixMarket::Output::operator()(             \
        std::span<const T> array) {                                         \
            return this->operator()(array.data(), array.size());            \
    }                                                                       \
    CIE_MATRIX_MARKET_OUTPUT_DEFINITION(int, T)                             \
    CIE_MATRIX_MARKET_OUTPUT_DEFINITION(unsigned, T)                        \
    CIE_MATRIX_MARKET_OUTPUT_DEFINITION(std::size_t, T)


CIE_MATRIX_MARKET_OUTPUT_DEFINITION_FOR_VALUE(float)
CIE_MATRIX_MARKET_OUTPUT_DEFINITION_FOR_VALUE(double)
#undef CIE_MATRIX_MARKET_OUTPUT_DEFINITION
#undef CIE_MATRIX_MARKET_OUTPUT_DEFINITION_FOR_VALUE


} // namespace cie::io

