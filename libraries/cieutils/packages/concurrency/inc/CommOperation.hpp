#ifndef CIE_UTILS_MPI_COMM_OPERATION_HPP
#define CIE_UTILS_MPI_COMM_OPERATION_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/MPI.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::mpi {


/** @brief Base class for wrapped RAII MPI communication.
 *  @details Derived classes should implement their represented communication
 *           such that it begins during construction and finished at destruction.
 */
class CommOperation
{
public:
    CommOperation(RightRef<std::shared_ptr<MPI>> r_mpi);

    CommOperation(CommOperation&& r_rhs) noexcept = default;

    CommOperation& operator=(CommOperation&& r_rhs) noexcept = default;

    virtual bool isFinished() const = 0;

    virtual ~CommOperation() = 0;

protected:
    Ref<const MPI> getMPI() const noexcept;

    Ref<MPI> getMPI() noexcept;

    std::shared_ptr<const MPI> shareMPI() const;

    std::shared_ptr<MPI> shareMPI();

private:
    CommOperation(const CommOperation&) = delete;

    CommOperation& operator=(const CommOperation&) = delete;

private:
    std::shared_ptr<MPI> _p_mpi;
}; // class CommOperation


} // namespace cie::mpi


#endif
