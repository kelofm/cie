#ifndef CIE_MPI_RANK_ID_HPP
#define CIE_MPI_RANK_ID_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"


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

    friend bool operator==(RankID left, RankID right) noexcept;

    explicit operator std::string() const;

    operator int() const;

private:
    int getInt() const noexcept;

    Ref<int> getInt() noexcept;

private:
    int _id;
};


using Color = int;


} // namespace cie::mpi

#include "packages/concurrency/impl/RankID_impl.hpp"

#endif
