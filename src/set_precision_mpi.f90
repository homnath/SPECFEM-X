! precision parameters
! REVISION:
! HNG, July 07,2011
module set_precision_mpi
use set_precision
use mpi
implicit none
integer,parameter :: MPI_KINT=MPI_INTEGER
integer,parameter :: MPI_KREAL=MPI_DOUBLE_PRECISION
end module set_precision_mpi
!============================================
