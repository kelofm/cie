#ifndef CIE_UTILS_IO_MATRIX_MARKET_HPP
#define CIE_UTILS_IO_MATRIX_MARKET_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <memory> // unique_ptr
#include <iosfwd> // istream, ostream


namespace cie::utils::io {


/// @brief IO class for reading and writing matrices in @a MatrixMarket format.
/// @ingroup cieutils
struct MatrixMarket
{
public:
    struct Settings
    {}; // struct Settings

    class Input
    {}; // class Input

    class Output
    {
    public:
        /// @brief Construct a matrix market input object reading from @a stdin.
        Output();

        /// @brief Construct a matrix market input object reading from the provided stream.
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

        CIE_MATRIX_MARKET_OUTPUT_INTERFACE(unsigned, double);

        CIE_MATRIX_MARKET_OUTPUT_INTERFACE(std::size_t, double);

        #undef CIE_MATRIX_MARKET_OUTPUT_INTERFACE

    private:
        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Output
}; // struct MatrixMarket


} // namespace cie::utils::io


#endif
