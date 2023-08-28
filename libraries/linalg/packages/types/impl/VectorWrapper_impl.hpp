#ifndef CIE_LINALG_TYPES_VECTOR_WRAPPER_IMPL_HPP
#define CIE_LINALG_TYPES_VECTOR_WRAPPER_IMPL_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/resize.hpp"


namespace cie::linalg {


template <class WT, class VT, class IT, class CIT>
template <class ...Args>
VectorWrapper<WT,VT,IT,CIT>::VectorWrapper( Args&&... r_args ) :
    _wrapped( std::forward<Args>(r_args)... )
{
}


template <class WT, class VT, class IT, class CIT>
inline Size
VectorWrapper<WT,VT,IT,CIT>::size() const
{
    return this->_wrapped.size();
} 


template <class WT, class VT, class IT, class CIT>
inline void
VectorWrapper<WT,VT,IT,CIT>::resize( Size newSize )
{
    CIE_BEGIN_EXCEPTION_TRACING

    utils::resize( this->_wrapped, newSize );

    CIE_END_EXCEPTION_TRACING
}



template <class WT, class VT, class IT, class CIT>
typename VectorWrapper<WT,VT,IT,CIT>::const_iterator
VectorWrapper<WT,VT,IT,CIT>::begin() const
{
    return &this->_wrapped( 0 );
}


template <class WT, class VT, class IT, class CIT>
typename VectorWrapper<WT,VT,IT,CIT>::iterator
VectorWrapper<WT,VT,IT,CIT>::begin()
{
    return &this->_wrapped( 0 );
}


template <class WT, class VT, class IT, class CIT>
typename VectorWrapper<WT,VT,IT,CIT>::const_iterator
VectorWrapper<WT,VT,IT,CIT>::end() const
{
    return &this->_wrapped( this->_wrapped.size() - 1 );
}


template <class WT, class VT, class IT, class CIT>
typename VectorWrapper<WT,VT,IT,CIT>::iterator
VectorWrapper<WT,VT,IT,CIT>::end()
{
    return &this->_wrapped( this->_wrapped.size() - 1 );
}


template <class WT, class VT, class IT, class CIT>
inline VT
VectorWrapper<WT,VT,IT,CIT>::operator()( Size index ) const
{
    return this->_wrapped[index];
}


template <class WT, class VT, class IT, class CIT>
inline VT&
VectorWrapper<WT,VT,IT,CIT>::operator()( Size index )
{
    return this->_wrapped[index];
}


template <class WT, class VT, class IT, class CIT>
inline VT
VectorWrapper<WT,VT,IT,CIT>::operator[]( Size index ) const
{
    return this->_wrapped[index];
}


template <class WT, class VT, class IT, class CIT>
inline VT&
VectorWrapper<WT,VT,IT,CIT>::operator[]( Size index )
{
    return this->_wrapped[index];
}



} // namespace cie::linalg


#endif