module map_location
contains
!-------------------------------------------------------------------------------

!DESCRIPTION
! This routine maps the original location to the natural coordinates [-1 1].
!REVISION
! HNG, Feb 19,2016
!TODO
! - better to put in math_library.f90?
subroutine map_point2naturalhex8(xg,xp,xip,x,iter,errx)
use global
use math_constants
use math_library,only:norm,determinant,invert
use shape_library,only : shape_function_hex8p,dshape_function_hex8p
implicit none
real(kind=kreal),intent(in) :: xg(ndim,8) ! coordinates of the geometrical nodes
real(kind=kreal),intent(in) :: xp(ndim)
double precision,intent(out) :: xip(ndim)
double precision,intent(out) :: x(ndim)
integer,intent(out) :: iter
double precision,intent(out) :: errx

integer,parameter :: maxiter=10
double precision :: tol=1d-12

double precision :: jac(ndim,ndim),detjac,xi(ndim),dx(ndim,1),dxi(ndim,1)

double precision :: shape_hex8(8),dshape_hex8(8,ndim)

! scale the tolerance otherwise it may not reach the tolerance for large models
!if(.not.devel_nondim)then
!  tol=tol*maxsize_elmt
!else
!  !tol=tol
!endif
! NOTE: if we use scaled tolerance, sometimes it converges quickly but the
!  absoulte error is still very large. On the other hand,
!  if we use nonscaled tolerance,
!  convergence may not occur for very large models. Instead, we choose to use
!  the nonscaled tolerance, and if it does not converge it will run for the
!  fixed number of iterations maxiter, error is very low even for 5 iterations.

! initial guess
xi=(/ ZERO, ZERO, ZERO /)
do iter=1,maxiter
  call shape_function_hex8p(8,xi(1),xi(2),xi(3),shape_hex8)
  x(1)=dot_product(xg(1,:),shape_hex8)
  x(2)=dot_product(xg(2,:),shape_hex8)
  x(3)=dot_product(xg(3,:),shape_hex8)
  dx(:,1)=xp-x
  ! scale the error otherwise it may not reach the tolerance for large models
  ! SEE the note above
  errx=norm(dx(:,1))
  if(errx.le.tol)then
    ! correction for the range
    where(xi.gt. ONE)xi= ONE
    where(xi.lt.-ONE)xi=-ONE
    xip=xi
    return
  endif
  call dshape_function_hex8p(8,xi(1),xi(2),xi(3),dshape_hex8)
  jac=matmul(xg,dshape_hex8)
  detjac=determinant(jac)
  call invert(jac)
  dxi=matmul(jac,dx)
  xi=xi+dxi(:,1)
  ! correction for the range
  if(xi(1).gt.one)xi(1)=one
  if(xi(1).lt.-one)xi(1)=-one
  if(xi(2).gt.one)xi(2)=one
  if(xi(2).lt.-one)xi(2)=-one
  if(xi(3).gt.one)xi(3)=one
  if(xi(3).lt.-one)xi(3)=-one
enddo
!print*,norm(dx(:,1))
!print*,'xp:',xp
!print*,'xg:',xg
!write(errtag,*)'ERROR: point mapping doesn''t converge! Point may be outside   &
!& the domain'
return
end subroutine map_point2naturalhex8
!===============================================================================

!DESCRIPTION
! This routine maps the original location to the natural coordinates [-1 1].
!REVISION
! HNG, Oct 08,2018
!TODO
! - better to put in math_library.f90?
subroutine map_point2naturalquad4(xg,xp,xip,x,iter,errx)
use global
use math_constants
use math_library,only:norm,determinant,invert
use shape_library,only : shape_function_quad4p,dshape_function_quad4p
implicit none
real(kind=kreal),intent(in) :: xg(2,4) ! coordinates of the geometrical nodes
real(kind=kreal),intent(in) :: xp(2)
double precision,intent(out) :: xip(2)
double precision,intent(out) :: x(2)
integer,intent(out) :: iter
double precision,intent(out) :: errx

integer,parameter :: maxiter=10
double precision :: tol=1d-12

double precision :: jac(2,2),detjac,xi(2),dx(2,1),dxi(2,1)

double precision :: shape_quad4(4),dshape_quad4(4,2)

! scale the tolerance otherwise it may not reach the tolerance for large models
!if(.not.devel_nondim)then
!  tol=tol*maxsize_elmt
!else
!  !tol=tol
!endif
! NOTE: if we use scaled tolerance, sometimes it converges quickly but the
!  absoulte error is still very large. On the other hand,
!  if we use nonscaled tolerance,
!  convergence may not occur for very large models. Instead, we choose to use
!  the nonscaled tolerance, and if it does not converge it will run for the
!  fixed number of iterations maxiter, error is very low even for 5 iterations.

! initial guess
xi=(/ ZERO, ZERO /)
do iter=1,maxiter
  call shape_function_quad4p(4,xi(1),xi(2),shape_quad4)
  x(1)=dot_product(xg(1,:),shape_quad4)
  x(2)=dot_product(xg(2,:),shape_quad4)
  dx(:,1)=xp-x
  ! scale the error otherwise it may not reach the tolerance for large models
  ! SEE the note above
  errx=norm(dx(:,1))
  if(errx.le.tol)then
    ! correction for the range
    where(xi.gt. ONE)xi= ONE
    where(xi.lt.-ONE)xi=-ONE
    xip=xi
    return
  endif
  call dshape_function_quad4p(4,xi(1),xi(2),dshape_quad4)
  jac=matmul(xg,dshape_quad4)
  detjac=determinant(jac)
  call invert(jac)
  dxi=matmul(jac,dx)
  xi=xi+dxi(:,1)
  ! correction for the range
  if(xi(1).gt.one)xi(1)=one
  if(xi(1).lt.-one)xi(1)=-one
  if(xi(2).gt.one)xi(2)=one
  if(xi(2).lt.-one)xi(2)=-one
enddo
!print*,norm(dx(:,1))
!print*,'xp:',xp
!print*,'xg:',xg
!write(errtag,*)'ERROR: point mapping doesn''t converge! Point may be outside   &
!& the domain'
return
end subroutine map_point2naturalquad4
!===============================================================================

end module map_location
!===============================================================================
