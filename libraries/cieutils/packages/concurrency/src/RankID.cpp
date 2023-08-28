// --- Utility Includes ---
#include "packages/concurrency/inc/RankID.hpp"


namespace cie::mpi {


RankID::operator std::string() const
{
    return std::to_string(_id);
}


Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, RankID id)
{
    return r_stream << (std::string)id;
}


} // namespace cie::mpi
