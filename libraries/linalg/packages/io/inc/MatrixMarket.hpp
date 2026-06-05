#pragma once

// --- LinAlg Includes ---
#include "packages/utilities/inc/CSRView.hpp"

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <memory> // unique_ptr
#include <iosfwd> // istream, ostream
#include <span> // span


namespace cie::linalg::io {


/// @brief IO class for reading and writing matrices in @a MatrixMarket format.
/// @ingroup cieutils
struct MatrixMarket {
    struct Settings
    {}; // struct Settings

    class Input
    {}; // class Input

    class Output
    {
    public:
        /// @brief Construct a matrix market output object writing to @p stdout.
        Output();

        /// @brief Construct a matrix market output object writing to the provided stream.
        Output(Ref<std::ostream> rStream, Settings settings = {});

        /// @brief Default destructor required by PIMPL.
        ~Output();

        #define CIE_MATRIX_MARKET_OUTPUT_INTERFACE(I, T)    \
            Ref<Output> operator()(                         \
                const I rowCount,                           \
                const I columnCount,                        \
                const I nonzeroCount,                       \
                Ptr<const I> pRowExtents,                   \
                Ptr<const I> pColumnIndices,                \
                Ptr<const T> pNonzeros);                    \
            Ref<Output> operator()(CSRView<const T,const I> matrix);

        #define CIE_MATRIX_MARKET_OUTPUT_INTERFACE_FOR_VALUE(T) \
            Ref<Output> operator()(                             \
                Ptr<const T> itBegin,                           \
                const std::size_t size);                        \
            Ref<Output> operator()(std::span<const T> array);   \
            CIE_MATRIX_MARKET_OUTPUT_INTERFACE(int, T)          \
            CIE_MATRIX_MARKET_OUTPUT_INTERFACE(unsigned, T)     \
            CIE_MATRIX_MARKET_OUTPUT_INTERFACE(std::size_t, T)

        CIE_MATRIX_MARKET_OUTPUT_INTERFACE_FOR_VALUE(float)
        CIE_MATRIX_MARKET_OUTPUT_INTERFACE_FOR_VALUE(double)
        #undef CIE_MATRIX_MARKET_OUTPUT_INTERFACE_FOR_VALUE
        #undef CIE_MATRIX_MARKET_OUTPUT_INTERFACE

    private:
        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Output
}; // struct MatrixMarket


} // namespace cie::io
