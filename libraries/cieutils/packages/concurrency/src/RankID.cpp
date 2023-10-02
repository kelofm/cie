// --- Utility Includes ---
#include "packages/concurrency/inc/RankID.hpp"

// --- STL Includes ---
#include <iostream>


namespace cie::mpi {


RankID::operator std::string() const
{
    return std::to_string(_id);
}


Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, RankID id)
{
    return r_stream << int(id);
}


} // namespace cie::mpi
