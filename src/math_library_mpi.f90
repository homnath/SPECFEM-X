! MPI math library
module math_library_mpi
use global,only:nproc,logunit
use math_constants
use set_precision_mpi

private :: iallgather_scal,fallgather_scal
private :: iminscal,fminscal
private :: iminvec,fminvec
private :: isumscal,fsumscal
private :: imaxscal,fmaxscal
private :: imaxvec,fmaxvec

! all gather of a scalar in all processors
interface allgather_scal
  module procedure iallgather_scal
  module procedure fallgather_scal
end interface

! global sum of a scalar in all processors
interface sumscal
  module procedure isumscal
  module procedure fsumscal
end interface

! global maximum of a scalar in all processors
interface minscal
  module procedure iminscal
  module procedure fminscal
end interface

! global maximum of a scalar in all processors
interface maxscal
  module procedure imaxscal
  module procedure fmaxscal
end interface

! global maximum of a vector in all processors
interface maxvec
  module procedure imaxvec
  module procedure fmaxvec
end interface

! global minimum of a scalar in all processors
interface minvec
  module procedure iminvec
  module procedure fminvec
end interface
contains
!=======================================================
!=======================================================

function iallgather_scal(scal) result(allscal)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer,dimension(1,0:nproc-1) :: allscal
integer :: scal_mpi(1) ! input for MPI_ALLGATHER are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLGATHER(scal_mpi,1,MPI_INTEGER,allscal,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
return
end function iallgather_scal
!=======================================================

function fallgather_scal(scal) result(allscal)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal),dimension(1,0:nproc-1) :: allscal
real(kind=kreal) :: scal_mpi(1) ! input for MPI_ALLGATHER are arrays
real(kind=kreal) :: ierr

scal_mpi(1)=scal
call MPI_ALLGATHER(scal_mpi,1,MPI_KREAL,allscal,1,MPI_KREAL,MPI_COMM_WORLD,ierr)
return
end function fallgather_scal
!=======================================================

function iminscal(scal) result(gmin)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gmin
integer :: scal_mpi(1),gmin_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gmin_mpi,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
gmin=gmin_mpi(1)
return
end function iminscal
!=======================================================

function fminscal(scal) result(gmin)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gmin
real(kind=kreal) :: scal_mpi(1),gmin_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gmin_mpi,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
gmin=gmin_mpi(1)
return
end function fminscal
!=======================================================

function imaxscal(scal) result(gmax)
!
! this finds a global maximum of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gmax
integer :: scal_mpi(1),gmax_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr
scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gmax_mpi,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
gmax=gmax_mpi(1)
return
end function imaxscal
!=======================================================

function fmaxscal(scal) result(gmax)
!
! this finds a global maximum of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gmax
real(kind=kreal) :: scal_mpi(1),gmax_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gmax_mpi,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
gmax=gmax_mpi(1)
return
end function fmaxscal
!=======================================================

function imaxvec(vec,n) result(gmax)
implicit none
integer,intent(in) :: n
integer,intent(in)::vec(:)
integer :: gmax(n) ! Global
integer :: ncomp
integer :: lvec(n) ! Local
integer :: ierr

ncomp=size(vec)
if(ncomp.ne.n)then
 write(logunit,'(a,1x,i0,1x,i0)')'ERROR: size of the "vec" doesn''t match!',ncomp,n
 stop
endif
lvec=vec
call MPI_ALLREDUCE(lvec,gmax,n,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
return
end function imaxvec
!=======================================================

function fmaxvec(vec,n) result(gmax)
implicit none
integer,intent(in) :: n
real(kind=kreal),intent(in)::vec(n)
real(kind=kreal) :: gmax(n) ! Global
integer :: ncomp
real(kind=kreal) :: lvec(n) ! Local
integer :: ierr

ncomp=size(vec)
if(ncomp.ne.n)then
 write(logunit,'(a,1x,i0,1x,i0)')'ERROR: size of the "vec" doesn''t match!',ncomp,n
 stop
endif
lvec=vec
call MPI_ALLREDUCE(lvec,gmax,n,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
return
end function fmaxvec
!=======================================================

function iminvec(vec,n) result(gmin)
implicit none
integer,intent(in) :: n
integer,intent(in)::vec(n)
integer :: gmin(n) ! Global
integer :: ierr

call MPI_ALLREDUCE(vec,gmin,n,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
return
end function iminvec
!=======================================================

function fminvec(vec,n) result(gmin)
implicit none
integer,intent(in) :: n
real(kind=kreal),intent(in)::vec(n)
real(kind=kreal) :: gmin(n) ! Global
integer :: ierr

call MPI_ALLREDUCE(vec,gmin,n,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
return
end function fminvec
!=======================================================

function isumscal(scal) result(gsum)
!
! this finds a global summation of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gsum ! global
integer :: scal_mpi(1),gsum_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gsum_mpi,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
gsum=gsum_mpi(1)
return
end function isumscal
!=======================================================

function fsumscal(scal) result(gsum)
!
! this finds a global summation of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gsum ! global
real(kind=kreal) :: scal_mpi(1),gsum_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

scal_mpi(1)=scal
call MPI_ALLREDUCE(scal_mpi,gsum_mpi,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
gsum=gsum_mpi(1)
return
end function fsumscal
!=======================================================

function dot_product_par(vec1,vec2) result(gdot)
!
! this finds global dot product of TWO vectors across the processors
!
implicit none
real(kind=kreal),intent(in)::vec1(:),vec2(:)
real(kind=kreal) :: gdot ! global
real(kind=kreal) :: ldot_mpi(1),gdot_mpi(1) ! input for MPI_ALLREDUCE are arrays
integer :: ierr

! find local dot
ldot_mpi(1)=dot_product(vec1,vec2)
call MPI_ALLREDUCE(ldot_mpi,gdot_mpi,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
gdot=gdot_mpi(1)
return
end function dot_product_par
!=======================================================

end module math_library_mpi
