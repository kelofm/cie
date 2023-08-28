#ifndef CIE_LINALG_SPARSE_HPP
#define CIE_LINALG_SPARSE_HPP

// --- Internal Includes ---
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <vector>
#include <tuple>
#include <cstddef>


namespace cie::linalg {

template<concepts::Numeric IndexType, concepts::Numeric ValueType = Double>
class CompressedSparseRowMatrix
{
public:
    typedef IndexType index_type;
    typedef ValueType value_type;

    explicit CompressedSparseRowMatrix();

    void allocate(const std::vector<std::vector<Size>>& locationMaps);

    Size rowSize() const;

    Size columnSize() const;

    IndexType nnz() const;

    ValueType operator()(Size i, Size j) const;
    std::vector<ValueType> operator*(const std::vector<ValueType>& vector);

    void scatter(const DynamicEigenMatrix<ValueType>& elementMatrix,
                 const std::vector<Size>& locationMap);

    //! Obtain sparse data structure and set matrix to zero
    std::tuple<IndexType*, IndexType*, ValueType*, Size> release();

private:
    IndexType* indices_;

    IndexType* indptr_;

    ValueType* data_;

    Size rowSize_, columnSize_;

    void cleanup();

public:
    ~CompressedSparseRowMatrix();

    // disable copying
    CompressedSparseRowMatrix(const CompressedSparseRowMatrix&) = delete;

    CompressedSparseRowMatrix& operator=(const CompressedSparseRowMatrix&) = delete;
};

}

#include "packages/types/impl/sparse_impl.hpp"

#endif // CIE_LINALG_SPARSE_HPP
