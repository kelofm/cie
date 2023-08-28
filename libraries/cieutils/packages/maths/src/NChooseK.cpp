// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/maths/inc/NChooseK.hpp"

// --- STL Includes ---
#include <numeric>
#include <algorithm>


namespace cie::utils {


NChooseK::NChooseK( Size n, Size k ) :
    _state( k ),
    _mask( n, 0 )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK(
        k <= n,
        "Cannot choose " + std::to_string(k) + " from " + std::to_string(n) + " elements!"
    )

    this->reset();

    CIE_END_EXCEPTION_TRACING
}


bool NChooseK::next()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( !std::prev_permutation(this->_mask.begin(), this->_mask.end()) )
        return false;

    auto it_state = this->_state.begin();

    for ( Size i=0; i<this->n(); ++i )
        if ( this->_mask[i] )
            *it_state++ = i;

    return true;

    CIE_END_EXCEPTION_TRACING
}


void NChooseK::reset()
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::iota(
        this->_state.begin(),
        this->_state.end(),
        0
    );

    std::fill_n(
        this->_mask.begin(),
        this->k(),
        1
    );

    std::fill(
        this->_mask.begin() + this->k(),
        this->_mask.end(),
        0
    );

    CIE_END_EXCEPTION_TRACING
}


NChooseK::state_container::const_iterator NChooseK::begin() const
{
    return this->_state.begin();
}


NChooseK::state_container::const_iterator NChooseK::end() const
{
    return this->_state.end();
}


Size NChooseK::n() const
{
    return this->_mask.size();
}


Size NChooseK::k() const
{
    return this->_state.size();
}


Size NChooseK::numberOfPermutations() const
{
    return this->nChooseK( this->n(), this->k() );
}


Size NChooseK::nChooseK( Size n, Size k )
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( k == 0 )
        return 1;

    return n * NChooseK::nChooseK( n-1, k-1 ) / k;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils