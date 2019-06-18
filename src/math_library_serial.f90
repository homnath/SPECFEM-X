! Serial math library this equivalent to math_library_mpi
module math_library_serial
use math_constants
use global,only:nproc
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
!-------------------------------------------------------------------------------

function iallgather_scal(scal) result(allscal)
!
! this finds a global minimum of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer,dimension(1,0:nproc-1) :: allscal
integer :: scal_mpi(1) ! input for MPI_ALLGATHER are arrays
integer :: ierr

allscal=scal
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

allscal=scal
return
end function fallgather_scal
!=======================================================

function iminscal(scal) result(gmin)
!
! this finds a summation of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gmin

gmin=scal

return
end function iminscal
!=======================================================

function fminscal(scal) result(gmin)
!
! this finds a summation of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gmin

gmin=scal

return
end function fminscal
!=======================================================

function imaxscal(scal) result(gmax)
!
! this finds a summation of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gmax

gmax=scal

return
end function imaxscal
!=======================================================

function fmaxscal(scal) result(gmax)
!
! this finds a summation of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gmax

gmax=scal

return
end function fmaxscal
!=======================================================

function imaxvec(vec,n) result(gmax)
implicit none
integer,intent(in) :: n
integer,intent(in)::vec(:)
integer :: lmax(n),gmax(n) ! local and global

lmax=vec

gmax=lmax

return
end function imaxvec
!=======================================================

function fmaxvec(vec,n) result(gmax)
implicit none
integer,intent(in) :: n
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: lmax(n),gmax(n) ! local and global

lmax=vec

gmax=lmax

return
end function fmaxvec
!=======================================================

function iminvec(vec,n) result(gmin)
implicit none
integer,intent(in) :: n
integer,intent(in)::vec(:)
integer :: lmin(n),gmin(n) ! local and global

lmin=vec

gmin=lmin

return
end function iminvec
!=======================================================

function fminvec(vec,n) result(gmin)
implicit none
integer,intent(in) :: n
real(kind=kreal),intent(in)::vec(:)
real(kind=kreal) :: lmin(n),gmin(n) ! local and global

lmin=vec

gmin=lmin

return
end function fminvec
!=======================================================

function isumscal(scal) result(gsum)
!
! this finds a summation of a scalar across the processors
!
implicit none
integer,intent(in)::scal
integer :: gsum

gsum=scal

return
end function isumscal
!=======================================================

function fsumscal(scal) result(gsum)
!
! this finds a summation of a scalar across the processors
!
implicit none
real(kind=kreal),intent(in)::scal
real(kind=kreal) :: gsum

gsum=scal

return
end function fsumscal
!=======================================================

function dot_product_par(vec1,vec2) result(gdot)
!
! this finds dot product of two vectors across the processors
!
implicit none
real(kind=kreal),intent(in)::vec1(:),vec2(:)
real(kind=kreal) :: ldot,gdot

! find local dot
ldot=dot_product(vec1,vec2)

gdot=ldot

return
end function dot_product_par
!=======================================================

end module math_library_serial
