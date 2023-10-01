#ifndef CIE_MPI_RANK_ID_HPP
#define CIE_MPI_RANK_ID_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <iosfwd>


namespace cie::mpi {


class RankID
{
private:
    friend class MPI;

    friend class MPISerialImpl;

    friend class MPIDistributedImpl;

public:
    explicit RankID(int id);

    RankID(RankID&& r_rhs) noexcept = default;

    RankID(const RankID& r_rhs) noexcept = default;

    RankID& operator=(RankID&& r_rhs) noexcept = default;

    RankID& operator=(const RankID& r_rhs) noexcept = default;

    bool isMain() const noexcept;

    friend bool operator==(RankID left, RankID right) noexcept;

    explicit operator std::string() const;

    explicit operator int() const;

private:
    int getInt() const noexcept;

private:
    int _id;
};


Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, RankID id);


using Color = int;


} // namespace cie::mpi

#include "packages/concurrency/impl/RankID_impl.hpp"

#endif
