#ifndef CIE_LINALG_EIGEN_MATRIX_IMPL_HPP
#define CIE_LINALG_EIGEN_MATRIX_IMPL_HPP

// --- Linalg Includes ---
#include "packages/matrix/inc/EigenMatrix.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include <packages/macros/inc/exceptions.hpp>


namespace cie::linalg {


template <class TMatrix>
EigenMatrix<TMatrix>::EigenMatrix(const EigenMatrix<TMatrix>& r_rhs)
    : _wrapped(r_rhs._wrapped)
{
}


template <class TMatrix>
EigenMatrix<TMatrix>::EigenMatrix(EigenMatrix<TMatrix>&& r_rhs) noexcept
    : _wrapped(std::move(r_rhs._wrapped))
{
}


template <class TMatrix>
EigenMatrix<TMatrix>::EigenMatrix(TMatrix&& r_wrapped)
    : _wrapped(std::move(r_wrapped))
{
}


template <class TMatrix>
EigenMatrix<TMatrix>::EigenMatrix(const TMatrix& r_wrapped)
    : _wrapped(r_wrapped)
{
}


template <class TMatrix>
EigenMatrix<TMatrix>::EigenMatrix(typename EigenMatrix<TMatrix>::size_type rowSize,
                                     typename EigenMatrix<TMatrix>::size_type columnSize)
requires concepts::ConstructibleFrom<TMatrix, typename EigenMatrix<TMatrix>::size_type, typename EigenMatrix<TMatrix>::size_type>
    : _wrapped(rowSize, columnSize)
{
}


template <class TMatrix>
EigenMatrix<TMatrix>::EigenMatrix(size_type size)
requires (EigenMatrix<TMatrix>::RowTag == 1 || EigenMatrix<TMatrix>::ColumnTag == 1)
    : _wrapped()
{
    this->resize(size);
}


template <class TMatrix>
EigenMatrix<TMatrix>::EigenMatrix(Ref<const std::initializer_list<value_type>> r_components)
    : _wrapped()
{
    if constexpr (RowTag != Eigen::Dynamic && ColumnTag != Eigen::Dynamic) {
        this->resize(RowTag, ColumnTag);
    } else if constexpr (RowTag != Eigen::Dynamic) {
        const size_type columns = r_components.size() / RowTag;
        CIE_OUT_OF_RANGE_CHECK(RowTag * columns == r_components.size())
        this->resize(RowTag, columns);
    }  else if constexpr (ColumnTag != Eigen::Dynamic) {
        const size_type rows = r_components.size() / ColumnTag;
        CIE_OUT_OF_RANGE_CHECK(rows * ColumnTag == r_components.size())
        this->resize(rows, ColumnTag);
    } else {
        this->resize(r_components.size(),1);
    }
    std::copy(r_components.begin(),
              r_components.end(),
              this->begin());
}



template <class TMatrix>
inline EigenMatrix<TMatrix>&
EigenMatrix<TMatrix>::operator=(TMatrix&& r_wrapped)
{
    CIE_BEGIN_EXCEPTION_TRACING
    _wrapped = std::move(r_wrapped);
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <class TMatrix>
inline EigenMatrix<TMatrix>&
EigenMatrix<TMatrix>::operator=(const TMatrix& r_wrapped)
{
    CIE_BEGIN_EXCEPTION_TRACING
    _wrapped = r_wrapped;
    return *this;
    CIE_END_EXCEPTION_TRACING
}


template <class TMatrix>
inline EigenMatrix<TMatrix>&
EigenMatrix<TMatrix>::operator=(const TMatrix&& r_wrapped)
{
    CIE_BEGIN_EXCEPTION_TRACING

    _wrapped = r_wrapped;
    return *this;

    CIE_END_EXCEPTION_TRACING
}


template <class TMatrix>
inline Ref<typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::at(size_type index) noexcept
{
    return *(this->begin() + index);
}


template <class TMatrix>
inline Ref<const typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::at(size_type index) const noexcept
{
    return *(this->begin() + index);
}


template <class TMatrix>
inline Ref<typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::operator[](size_type index) noexcept
{
    return *(this->begin() + index);
}


template <class TMatrix>
inline Ref<const typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::operator[](size_type index) const noexcept
{
    return *(this->begin() + index);
}


template <class TMatrix>
inline typename EigenMatrix<TMatrix>::value_type&
EigenMatrix<TMatrix>::at(typename EigenMatrix<TMatrix>::size_type rowIndex,
                             typename EigenMatrix<TMatrix>::size_type columnIndex)
{
    CIE_OUT_OF_RANGE_CHECK(rowIndex < this->rowSize())
    CIE_OUT_OF_RANGE_CHECK(columnIndex < this->columnSize())

    return _wrapped(rowIndex, columnIndex);
}


template <class TMatrix>
inline typename EigenMatrix<TMatrix>::value_type
EigenMatrix<TMatrix>::at(typename EigenMatrix<TMatrix>::size_type rowIndex,
                             typename EigenMatrix<TMatrix>::size_type columnIndex) const
{
    CIE_OUT_OF_RANGE_CHECK(rowIndex < this->rowSize())
    CIE_OUT_OF_RANGE_CHECK(columnIndex < this->columnSize())

    return _wrapped(rowIndex, columnIndex);
}


template <class TMatrix>
inline typename EigenMatrix<TMatrix>::value_type&
EigenMatrix<TMatrix>::operator()(typename EigenMatrix<TMatrix>::size_type rowIndex,
                                     typename EigenMatrix<TMatrix>::size_type columnIndex)
{
    return this->at(rowIndex, columnIndex);
}


template <class TMatrix>
inline typename EigenMatrix<TMatrix>::value_type
EigenMatrix<TMatrix>::operator()(typename EigenMatrix<TMatrix>::size_type rowIndex,
                                     typename EigenMatrix<TMatrix>::size_type columnIndex) const
{
    return this->at(rowIndex, columnIndex);
}


template <class TMatrix>
inline void
EigenMatrix<TMatrix>::resize(typename EigenMatrix<TMatrix>::size_type rowSize,
                                 typename EigenMatrix<TMatrix>::size_type columnSize)
requires concepts::detail::HasResize<TMatrix,typename EigenMatrix<TMatrix>::size_type,typename EigenMatrix<TMatrix>::size_type>
{
    CIE_BEGIN_EXCEPTION_TRACING
    _wrapped.resize(rowSize, columnSize);
    CIE_END_EXCEPTION_TRACING
}


template <class TMatrix>
inline void
EigenMatrix<TMatrix>::resize(typename EigenMatrix<TMatrix>::size_type rowSize,
                                 typename EigenMatrix<TMatrix>::size_type columnSize)
requires (!concepts::detail::HasResize<TMatrix,typename EigenMatrix<TMatrix>::size_type,typename EigenMatrix<TMatrix>::size_type>)
{
    CIE_CHECK(rowSize == this->rowSize() && columnSize == this->columnSize(),
              "Attempt to resize a static matrix from ("
              << this->rowSize() << 'x' << this->columnSize()
              << ") to (" << rowSize << 'x' << columnSize << ')')
}


template <class TMatrix>
inline void
EigenMatrix<TMatrix>::resize(size_type size)
requires (EigenMatrix<TMatrix>::RowTag == 1 || EigenMatrix<TMatrix>::ColumnTag == 1)
{
    CIE_BEGIN_EXCEPTION_TRACING
    if constexpr (ColumnTag == 1) {
        this->resize(size, 1);
    } else {
        this->resize(1, size);
    }
    CIE_END_EXCEPTION_TRACING
}


template <class TMatrix>
inline typename EigenMatrix<TMatrix>::size_type
EigenMatrix<TMatrix>::rowSize() const noexcept
{
    return _wrapped.rows();
}


template <class TMatrix>
inline typename EigenMatrix<TMatrix>::size_type
EigenMatrix<TMatrix>::columnSize() const noexcept
{
    return _wrapped.cols();
}


template <class TMatrix>
inline typename EigenMatrix<TMatrix>::size_type
EigenMatrix<TMatrix>::size() const noexcept
{
    if constexpr (RowTag != Eigen::Dynamic && ColumnTag != Eigen::Dynamic) {
        return RowTag * ColumnTag;
    } else {
        return this->rowSize() * this->columnSize();
    }
}


template <class TMatrix>
inline Ptr<const typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::begin() const noexcept
{
    return this->data();
}


template <class TMatrix>
inline Ptr<typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::begin() noexcept
{
    return this->data();
}


template <class TMatrix>
inline Ptr<const typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::end() const noexcept
{
    return this->data() + this->size();
}


template <class TMatrix>
inline Ptr<typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::end() noexcept
{
    return this->data() + this->size();
}


template <class TMatrix>
inline Ptr<const typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::data() const noexcept
{
    return _wrapped.data();
}


template <class TMatrix>
inline Ptr<typename EigenMatrix<TMatrix>::value_type>
EigenMatrix<TMatrix>::data() noexcept
{
    return _wrapped.data();
}


template <class TMatrix>
inline TMatrix&
EigenMatrix<TMatrix>::wrapped()
{
    return _wrapped;
}


template <class TMatrix>
inline const TMatrix&
EigenMatrix<TMatrix>::wrapped() const
{
    return _wrapped;
}


template <class TMatrix>
inline
EigenMatrix<TMatrix>::operator TMatrix&()
{
    return this->wrapped();
}


template <class TMatrix>
inline
EigenMatrix<TMatrix>::operator const TMatrix&() const
{
    return this->wrapped();
}


template <class TMatrix>
inline std::ostream&
operator<<(std::ostream& r_stream, const EigenMatrix<TMatrix>& r_matrix)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return r_stream << r_matrix.wrapped();

    CIE_END_EXCEPTION_TRACING
}


template <class TMatrix>
EigenMatrix<TMatrix>
EigenMatrix<TMatrix>::makeZeroMatrix()
{
    return EigenMatrix<TMatrix>(TMatrix::Zero());
}


template <class TMatrix>
EigenMatrix<TMatrix>
EigenMatrix<TMatrix>::makeIdentityMatrix()
{
    return EigenMatrix<TMatrix>(TMatrix::Identity());
}


} // namespace cie::linalg


#endif