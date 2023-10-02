#ifndef CIE_MPI_RANKID_IMPL_HPP
#define CIE_MPI_RANKID_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/RankID.hpp"


namespace cie::mpi {


inline RankID::RankID(int id)
    : _id(id)
{
}


inline bool RankID::isMain() const noexcept
{
    return _id == 0;
}


inline bool operator==(RankID left, RankID right) noexcept
{
    return left._id == right._id;
}


inline RankID::operator int() const
{
    return _id;
}


inline int RankID::getInt() const noexcept
{
    return _id;
}


} // namespace cie::mpi


#endif
