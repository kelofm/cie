#ifndef CIE_UTILS_IO_MATRIX_MARKET_HPP
#define CIE_UTILS_IO_MATRIX_MARKET_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <memory> // unique_ptr
#include <iosfwd> // istream, ostream


namespace cie::io {


/// @brief IO class for reading and writing matrices in @a MatrixMarket format.
/// @ingroup cieutils
struct MatrixMarket
{
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

        #define CIE_MATRIX_MARKET_OUTPUT_INTERFACE(TIndex, TValue)      \
            Ref<Output> operator()(const TIndex rowCount,               \
                                   const TIndex columnCount,            \
                                   const TIndex nonzeroCount,           \
                                   Ptr<const TIndex> pRowExtents,       \
                                   Ptr<const TIndex> pColumnIndices,    \
                                   Ptr<const TValue> pNonzeros)

        #define CIE_MATRIX_MARKET_OUTPUT_INTERFACE_FOR_VALUE(TValue)    \
            Ref<Output> operator()(Ptr<const TValue> itBegin,           \
                                   const std::size_t size);             \
                                                                        \
            CIE_MATRIX_MARKET_OUTPUT_INTERFACE(int, TValue);            \
            CIE_MATRIX_MARKET_OUTPUT_INTERFACE(unsigned, TValue);       \
            CIE_MATRIX_MARKET_OUTPUT_INTERFACE(std::size_t, TValue)

        CIE_MATRIX_MARKET_OUTPUT_INTERFACE_FOR_VALUE(float);
        CIE_MATRIX_MARKET_OUTPUT_INTERFACE_FOR_VALUE(double);
        #undef CIE_MATRIX_MARKET_OUTPUT_INTERFACE_FOR_VALUE
        #undef CIE_MATRIX_MARKET_OUTPUT_INTERFACE

    private:
        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Output
}; // struct MatrixMarket


} // namespace cie::io


#endif
