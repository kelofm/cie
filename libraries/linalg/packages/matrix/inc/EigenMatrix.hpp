#ifndef CIE_LINALG_EIGEN_MATRIX_HPP
#define CIE_LINALG_EIGEN_MATRIX_HPP

// --- External Inlcudes ---
#include "Eigen/Dense"

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- STL Includes ---
#include <ostream>


namespace cie::linalg {


template <class TMatrix>
class EigenMatrix
{
public:
    using Wrapped = TMatrix;

    using value_type = typename TMatrix::value_type;

    using size_type  = Size;

    static const int RowTag    = TMatrix::RowsAtCompileTime;

    static const int ColumnTag = TMatrix::ColsAtCompileTime;

public:
    EigenMatrix() noexcept = default;

    EigenMatrix(const EigenMatrix<TMatrix>& r_rhs);

    EigenMatrix(EigenMatrix<TMatrix>&& r_rhs) noexcept;

    EigenMatrix(TMatrix&& r_wrapped);

    EigenMatrix(const TMatrix& r_wrapped);

    EigenMatrix(size_type rowSize, size_type columnSize)
    requires concepts::ConstructibleFrom<TMatrix, size_type, size_type>;

    EigenMatrix(size_type size)
    requires (RowTag == 1 || ColumnTag == 1);

    EigenMatrix(Ref<const std::initializer_list<value_type>> r_components);

    EigenMatrix<TMatrix>& operator=(EigenMatrix<TMatrix>&& r_rhs) = default;

    EigenMatrix<TMatrix>& operator=(const EigenMatrix<TMatrix>& r_rhs) = default;

    EigenMatrix<TMatrix>& operator=(TMatrix&& r_wrapped);

    EigenMatrix<TMatrix>& operator=(const TMatrix& r_wrapped);

    EigenMatrix<TMatrix>& operator=(const TMatrix&& r_wrapped);

    Ref<value_type> at(size_type index) noexcept;

    Ref<const value_type> at(size_type) const noexcept;

    value_type& at(size_type rowIndex,
                   size_type columnIndex);

    value_type at(size_type rowIndex,
                  size_type columnIndex) const;

    Ref<value_type> operator[](size_type index) noexcept;

    Ref<const value_type> operator[](size_type index) const noexcept;

    value_type& operator()(size_type rowIndex,
                           size_type columnIndex);

    value_type operator()(size_type rowIndex,
                          size_type columnIndex) const;

    void resize(size_type rowSize, size_type columnSize)
    requires concepts::detail::HasResize<TMatrix,size_type,size_type>;

    void resize(size_type rowSize, size_type columnSize)
    requires (!concepts::detail::HasResize<TMatrix,size_type,size_type>);

    void resize(size_type size)
    requires (RowTag == 1 || ColumnTag == 1);

    Size rowSize() const noexcept;

    Size columnSize() const noexcept;

    Size size() const noexcept;

    Ptr<const value_type> begin() const noexcept;

    Ptr<value_type> begin() noexcept;

    Ptr<const value_type> end() const noexcept;

    Ptr<value_type> end() noexcept;

    Ptr<const value_type> data() const noexcept;

    Ptr<value_type> data() noexcept;

    TMatrix& wrapped();

    const TMatrix& wrapped() const;

    explicit operator TMatrix& ();

    explicit operator const TMatrix& () const;

    static EigenMatrix makeZeroMatrix();

    static EigenMatrix makeIdentityMatrix();

protected:
    TMatrix _wrapped;
};


template <class TMatrix>
std::ostream& operator<<(std::ostream& r_stream, const EigenMatrix<TMatrix>& r_matrix);


} // namespace cie::linalg


namespace cie::concepts {

template <class T>
concept EigenMatrix
= requires (T instance)
{
    typename T::value_type;
    typename T::size_type;
    {T::RowTag};
    {T::ColumnTag};
    {instance.at(typename T::size_type(), typename T::size_type())};
    {instance(typename T::size_type(), typename T::size_type())};
    {instance.rowSize()};
    {instance.columnSize()};
};

} // namespace cie::concepts


#include "packages/matrix/impl/EigenMatrix_impl.hpp"

#endif
