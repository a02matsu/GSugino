module parallel
use mpi
implicit none

!include 'mpif.h'
integer, save :: NPROCS,MYRANK,IERR
integer, save :: ISTATUS(MPI_STATUS_SIZE) ! for MPI_RECV
integer :: MPI_GENRAND_SREPR, MPI_LOCAL_LABEL


end module parallel
