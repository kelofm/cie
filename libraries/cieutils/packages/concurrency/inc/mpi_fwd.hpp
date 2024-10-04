#ifndef CIE_UTILS_MPI_FWD_HPP
#define CIE_UTILS_MPI_FWD_HPP


/** @brief This is terrible.
 *
 *  @details Pretty much every meaningful MPI call requires a communicator
 *           of type @p MPI_Comm, except @p MPI_Comm is not a type but a typedef
 *           (at least in OpenMPI), which means that the underlying type
 *           must be forward declared, and then MPI_Comm needs to be typedef'd
 *           to that. Problem is, the underlying type is not part of the
 *           public API. As a result, there are two options here:
 *           0) to hell with forward declarations, @code #include "mpi.h" @endcode and
 *              bleed the whole fucking MPI interface into the project.
 *           1) Restrict the set of supported MPI implementations and find
 *              all underlying types, then forward declare them based on
 *              which version CiE is compiled with.
 *
 *           Great, @todo time. The following is only valid for the OpenMPI
 *           implementation.
 */
#ifdef CIEUTILS_COMPILE_MPI
    typedef struct ompi_communicator_t *MPI_Comm; // <== OpenMPI
#else
    struct MPI_Comm {};
#endif


#endif
