#ifndef CIE_LINALG_TYPES_MATRIX_WRAPPER_IMPL_HPP
#define CIE_LINALG_TYPES_MATRIX_WRAPPER_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"


namespace cie::linalg {


template <class MT, class VT, class IT, class CIT>
MatrixWrapper<MT,VT,IT,CIT>::MatrixWrapper( const MatrixWrapper<MT,VT,IT,CIT>& r_rhs ) :
    _wrapped( r_rhs._wrapped )
{
}


template <class MT, class VT, class IT, class CIT>
MatrixWrapper<MT,VT,IT,CIT>::MatrixWrapper( MatrixWrapper<MT,VT,IT,CIT>&& r_rhs ) :
    _wrapped( std::move(r_rhs._wrapped) )
{
}


template <class MT, class VT, class IT, class CIT>
MatrixWrapper<MT,VT,IT,CIT>::MatrixWrapper( const MatrixWrapper<MT,VT,IT,CIT>&& r_rhs ) :
    _wrapped( std::move(r_rhs._wrapped) )
{
}


template <class MT, class VT, class IT, class CIT>
template <class ...Args>
MatrixWrapper<MT,VT,IT,CIT>::MatrixWrapper( Args&&... r_args ) :
    _wrapped( std::forward<Args>(r_args)... )
{
}


template <class MT, class VT, class IT, class CIT>
inline IT
MatrixWrapper<MT,VT,IT,CIT>::begin()
{
    return IT( &this->_wrapped.operator()( 0, 0 ) );
}


template <class MT, class VT, class IT, class CIT>
inline CIT
MatrixWrapper<MT,VT,IT,CIT>::begin() const
{
    return CIT( &this->_wrapped.operator()( 0, 0 ) );
}


template <class MT, class VT, class IT, class CIT>
inline IT
MatrixWrapper<MT,VT,IT,CIT>::end()
{
    return IT( this->begin() + this->size() );
}


template <class MT, class VT, class IT, class CIT>
inline CIT
MatrixWrapper<MT,VT,IT,CIT>::end() const
{
    return CIT( this->begin() + this->size() );
}


template <class MT, class VT, class IT, class CIT>
inline typename MatrixWrapper<MT,VT,IT,CIT>::size_type
MatrixWrapper<MT,VT,IT,CIT>::size() const
{
    return this->_wrapped.size();
}


template <class MT, class VT, class IT, class CIT>
inline typename MatrixWrapper<MT,VT,IT,CIT>::size_type
MatrixWrapper<MT,VT,IT,CIT>::numberOfRows() const
{
    return this->_wrapped.rows();
}


template <class MT, class VT, class IT, class CIT>
inline typename MatrixWrapper<MT,VT,IT,CIT>::size_type
MatrixWrapper<MT,VT,IT,CIT>::numberOfColumns() const
{
    return this->_wrapped.cols();
}


template <class MT, class VT, class IT, class CIT>
inline VT&
MatrixWrapper<MT,VT,IT,CIT>::operator()( typename MatrixWrapper<MT,VT,IT,CIT>::size_type rowIndex,
                                         typename MatrixWrapper<MT,VT,IT,CIT>::size_type columnIndex )
{
    CIE_OUT_OF_RANGE_CHECK( rowIndex < this->numberOfRows() )
    CIE_OUT_OF_RANGE_CHECK( columnIndex < this->numberOfColumns() )

    return this->_wrapped( rowIndex, columnIndex );
}


template <class MT, class VT, class IT, class CIT>
inline VT
MatrixWrapper<MT,VT,IT,CIT>::operator()( typename MatrixWrapper<MT,VT,IT,CIT>::size_type rowIndex,
                                         typename MatrixWrapper<MT,VT,IT,CIT>::size_type columnIndex ) const
{
    CIE_OUT_OF_RANGE_CHECK( rowIndex < this->numberOfRows() )
    CIE_OUT_OF_RANGE_CHECK( columnIndex < this->numberOfColumns() )

    return this->_wrapped( rowIndex, columnIndex );
}


template <class MT, class VT, class IT, class CIT>
MatrixWrapper<MT,VT,IT,CIT>::operator const MT&() const
{
    return this->_wrapped;
}


template <class MT, class VT, class IT, class CIT>
MatrixWrapper<MT,VT,IT,CIT>::operator MT&()
{
    return this->_wrapped;
}


} // namespace cie::linalg

#endif