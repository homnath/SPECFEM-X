! DESCRIPTION
!  This module contains routines to compute shape functions and 
!  their derivatives for hexhedral and quadrilateral elements.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Jul 12,2011; HNG, Apr 09,2010
! TODO
module shape_library
use set_precision
contains
!-------------------------------------------------------------------------------

! This subroutines computes the shape fucntions at gll
! points. the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention
subroutine shape_function_hex8(ngnod,ngllx,nglly,ngllz,xigll,etagll,      &
zetagll,shape_hex8)
use set_precision
use math_constants
implicit none

integer,intent(in) :: ngnod,ngllx,nglly,ngllz

! gauss-lobatto-legendre points of integration
double precision :: xigll(ngllx)
double precision :: etagll(nglly)
double precision :: zetagll(ngllz)

! 3d shape functions
double precision :: shape_hex8(ngnod,ngllx,nglly,ngllz)

integer :: i,j,k,i_gnod

! location of the nodes of the 3d quadrilateral elements
!double precision xi,eta,gamma
double precision :: xip,xim,etap,etam,zetap,zetam

! for checking the 3d shape functions
double precision :: sum_shape

double precision, parameter :: one_eighth = 0.125d0

! check that the parameter file is correct
if(ngnod /= 8)then
  write(*,*)'ERROR: elements must have 8 geometrical nodes!'
  stop
endif

! compute shape functions
! case of a 3d 8-node element (dhatt-touzot p. 115)
do k=1,ngllz
  do j=1,nglly
    do i=1,ngllx

      !xi = xigll(i)
      !eta = etagll(j)
      !gamma = zetagll(k)

      xip = one + xigll(i)
      xim = one - xigll(i)

      etap = one + etagll(j)
      etam = one - etagll(j)

      zetap = one + zetagll(k)
      zetam = one - zetagll(k)

      shape_hex8(1,i,j,k) = one_eighth*xim*etam*zetam
      shape_hex8(2,i,j,k) = one_eighth*xip*etam*zetam
      shape_hex8(3,i,j,k) = one_eighth*xip*etap*zetam
      shape_hex8(4,i,j,k) = one_eighth*xim*etap*zetam
      shape_hex8(5,i,j,k) = one_eighth*xim*etam*zetap
      shape_hex8(6,i,j,k) = one_eighth*xip*etam*zetap
      shape_hex8(7,i,j,k) = one_eighth*xip*etap*zetap
      shape_hex8(8,i,j,k) = one_eighth*xim*etap*zetap
    enddo
  enddo
enddo

! check the shape functions and their derivatives

do k=1,ngllz
  do j=1,nglly
    do i=1,ngllx

      sum_shape = zero

      do i_gnod=1,ngnod
        sum_shape = sum_shape + shape_hex8(i_gnod,i,j,k)
      enddo

      ! sum of shape functions should be one
      if(abs(sum_shape-one) >  zerotol)then
        write(*,*)'ERROR: error shape functions!'
        stop
      endif
    enddo
  enddo
enddo

end subroutine shape_function_hex8
!===============================================================================

! This subroutines computes derivatives of the shape fucntions at gll
! points. the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention
subroutine dshape_function_hex8(ngnod,ngllx,nglly,ngllz,xigll,etagll,     &
zetagll,dshape_hex8)
use set_precision
use math_constants
implicit none
integer,intent(in) :: ngnod,ngllx,nglly,ngllz

! gauss-lobatto-legendre points of integration
double precision :: xigll(ngllx)
double precision :: etagll(nglly)
double precision :: zetagll(ngllz)

! derivatives of the 3d shape functions
double precision :: dshape_hex8(3,ngnod,ngllx*nglly*ngllz)

integer :: i,j,k,i_gnod
integer :: igll,ngll

double precision :: xip,xim,etap,etam,zetap,zetam

! for checking the 3d shape functions
double precision :: sum_dshapexi,sum_dshapeeta,sum_dshapezeta

double precision, parameter :: one_eighth = 0.125_kreal

! check that the parameter file is correct
if(ngnod /= 8)then
  write(*,*)'ERROR: elements must have 8 geometrical nodes!'
  stop
endif

ngll=ngllx*nglly*ngllz

! compute the derivatives of 3d shape functions
igll=0
do k=1,ngllz
  zetap = one + zetagll(k)
  zetam = one - zetagll(k)
  do j=1,nglly
    etap = one + etagll(j)
    etam = one - etagll(j)
    do i=1,ngllx
      igll=igll+1

      !xi = xigll(i)
      !eta = etagll(j)
      !gamma = zetagll(k)

      xip = one + xigll(i)
      xim = one - xigll(i)

      dshape_hex8(1,1,igll) = - one_eighth*etam*zetam
      dshape_hex8(1,2,igll) = one_eighth*etam*zetam
      dshape_hex8(1,3,igll) = one_eighth*etap*zetam
      dshape_hex8(1,4,igll) = - one_eighth*etap*zetam
      dshape_hex8(1,5,igll) = - one_eighth*etam*zetap
      dshape_hex8(1,6,igll) = one_eighth*etam*zetap
      dshape_hex8(1,7,igll) = one_eighth*etap*zetap
      dshape_hex8(1,8,igll) = - one_eighth*etap*zetap

      dshape_hex8(2,1,igll) = - one_eighth*xim*zetam
      dshape_hex8(2,2,igll) = - one_eighth*xip*zetam
      dshape_hex8(2,3,igll) = one_eighth*xip*zetam
      dshape_hex8(2,4,igll) = one_eighth*xim*zetam
      dshape_hex8(2,5,igll) = - one_eighth*xim*zetap
      dshape_hex8(2,6,igll) = - one_eighth*xip*zetap
      dshape_hex8(2,7,igll) = one_eighth*xip*zetap
      dshape_hex8(2,8,igll) = one_eighth*xim*zetap

      dshape_hex8(3,1,igll) = - one_eighth*xim*etam
      dshape_hex8(3,2,igll) = - one_eighth*xip*etam
      dshape_hex8(3,3,igll) = - one_eighth*xip*etap
      dshape_hex8(3,4,igll) = - one_eighth*xim*etap
      dshape_hex8(3,5,igll) = one_eighth*xim*etam
      dshape_hex8(3,6,igll) = one_eighth*xip*etam
      dshape_hex8(3,7,igll) = one_eighth*xip*etap
      dshape_hex8(3,8,igll) = one_eighth*xim*etap

    enddo
  enddo
enddo

! check the shape functions and their derivatives

do i=1,ngll
      sum_dshapexi = zero
      sum_dshapeeta = zero
      sum_dshapezeta = zero

      do i_gnod=1,ngnod
        sum_dshapexi = sum_dshapexi + dshape_hex8(1,i_gnod,i)
        sum_dshapeeta = sum_dshapeeta + dshape_hex8(2,i_gnod,i)
        sum_dshapezeta = sum_dshapezeta + dshape_hex8(3,i_gnod,i)
      enddo

      ! sum of derivative of shape functions should be zero
      if(abs(sum_dshapexi) >  zerotol)then
        write(*,*)'ERROR: derivative xi shape functions!'
        stop
      endif
      if(abs(sum_dshapeeta) >  zerotol)then
        write(*,*)'ERROR: derivative eta shape functions!'
        stop
      endif
      if(abs(sum_dshapezeta) >  zerotol)then
        write(*,*)'ERROR: derivative gamma shape functions!'
        stop
      endif
enddo

end subroutine dshape_function_hex8
!===============================================================================

! This subroutines computes the shape fucntions at a given point (xi,eta,zeta)
! points. the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention
subroutine shape_function_hex8p(ngnod,xi,eta,zeta,shape_hex8)
use set_precision
use math_constants
implicit none

integer,intent(in) :: ngnod

! given point
double precision :: xi
double precision :: eta
double precision :: zeta

! 3d shape functions
double precision :: shape_hex8(ngnod)

integer :: i_gnod

double precision :: xip,xim,etap,etam,zetap,zetam

! for checking the 3d shape functions
double precision :: sum_shape

double precision, parameter :: one_eighth = 0.125d0

! check that the parameter file is correct
if(ngnod /= 8)then
  write(*,*)'ERROR: elements must have 8 geometrical nodes!'
  stop
endif

! compute shape functions
! case of a 3d 8-node element (dhatt-touzot p. 115)
! check point
if(abs(xi).gt.one .or. abs(eta).gt.one .or. abs(zeta).gt.one)then
  write(*,*)'WARNING: point is outside the domain (-1,1)!',xi,eta,zeta
  !stop
endif

xip = one + xi
xim = one - xi

etap = one + eta
etam = one - eta

zetap = one + zeta
zetam = one - zeta
shape_hex8=zero
shape_hex8(1) = one_eighth*xim*etam*zetam
shape_hex8(2) = one_eighth*xip*etam*zetam
shape_hex8(3) = one_eighth*xip*etap*zetam
shape_hex8(4) = one_eighth*xim*etap*zetam
shape_hex8(5) = one_eighth*xim*etam*zetap
shape_hex8(6) = one_eighth*xip*etam*zetap
shape_hex8(7) = one_eighth*xip*etap*zetap
shape_hex8(8) = one_eighth*xim*etap*zetap

! check the shape functions and their derivatives
sum_shape = zero

do i_gnod=1,ngnod
  sum_shape = sum_shape + shape_hex8(i_gnod)
enddo

! sum of shape functions should be one
if(abs(sum_shape-one) >  zerotol)then
  write(*,*)'ERROR: error shape functions!'
  stop
endif

end subroutine shape_function_hex8p
!===============================================================================

! This subroutines computes derivatives of the shape fucntions at given point.
! the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention
!NOTE: dimension of dshape_hex8 is (ngnod,3) NOT (3,ngnode)
subroutine dshape_function_hex8p(ngnod,xi,eta,zeta,dshape_hex8)
use set_precision
use math_constants
implicit none
integer,intent(in) :: ngnod

! given point
double precision :: xi
double precision :: eta
double precision :: zeta

! derivatives of the 3d shape functions
double precision :: dshape_hex8(ngnod,3)

integer :: i_gnod

double precision :: xip,xim,etap,etam,zetap,zetam

! for checking the 3d shape functions
double precision :: sum_dshapexi,sum_dshapeeta,sum_dshapezeta

double precision, parameter :: one_eighth = 0.125_kreal

! check that the parameter file is correct
if(ngnod /= 8)then
  write(*,*)'ERROR: elements must have 8 geometrical nodes!'
  stop
endif

! compute the derivatives of 3d shape functions
zetap = one + zeta
zetam = one - zeta
etap =  one + eta
etam =  one - eta
xip =   one + xi
xim =   one - xi

dshape_hex8=zero

dshape_hex8(1,1) = - one_eighth*etam*zetam
dshape_hex8(2,1) = one_eighth*etam*zetam
dshape_hex8(3,1) = one_eighth*etap*zetam
dshape_hex8(4,1) = - one_eighth*etap*zetam
dshape_hex8(5,1) = - one_eighth*etam*zetap
dshape_hex8(6,1) = one_eighth*etam*zetap
dshape_hex8(7,1) = one_eighth*etap*zetap
dshape_hex8(8,1) = - one_eighth*etap*zetap

dshape_hex8(1,2) = - one_eighth*xim*zetam
dshape_hex8(2,2) = - one_eighth*xip*zetam
dshape_hex8(3,2) = one_eighth*xip*zetam
dshape_hex8(4,2) = one_eighth*xim*zetam
dshape_hex8(5,2) = - one_eighth*xim*zetap
dshape_hex8(6,2) = - one_eighth*xip*zetap
dshape_hex8(7,2) = one_eighth*xip*zetap
dshape_hex8(8,2) = one_eighth*xim*zetap

dshape_hex8(1,3) = - one_eighth*xim*etam
dshape_hex8(2,3) = - one_eighth*xip*etam
dshape_hex8(3,3) = - one_eighth*xip*etap
dshape_hex8(4,3) = - one_eighth*xim*etap
dshape_hex8(5,3) = one_eighth*xim*etam
dshape_hex8(6,3) = one_eighth*xip*etam
dshape_hex8(7,3) = one_eighth*xip*etap
dshape_hex8(8,3) = one_eighth*xim*etap

! check the shape functions and their derivatives
sum_dshapexi = zero
sum_dshapeeta = zero
sum_dshapezeta = zero

do i_gnod=1,ngnod
  sum_dshapexi = sum_dshapexi + dshape_hex8(i_gnod,1)
  sum_dshapeeta = sum_dshapeeta + dshape_hex8(i_gnod,2)
  sum_dshapezeta = sum_dshapezeta + dshape_hex8(i_gnod,3)
enddo

! sum of derivative of shape functions should be zero
if(abs(sum_dshapexi) >  zerotol)then
  write(*,*)'ERROR: derivative xi shape functions!'
  stop
endif
if(abs(sum_dshapeeta) >  zerotol)then
  write(*,*)'ERROR: derivative eta shape functions!'
  stop
endif
if(abs(sum_dshapezeta) >  zerotol)then
  write(*,*)'ERROR: derivative gamma shape functions!'
  stop
endif

end subroutine dshape_function_hex8p
!===============================================================================

! This subroutines computes derivatives of the shape fucntions at gll
! points. the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention
subroutine dshape_function_quad4(ngnod2d,ngllx,nglly,xigll,etagll,dshape_quad4)
use set_precision
use math_constants
implicit none
integer,intent(in) :: ngnod2d,ngllx,nglly

! gauss-lobatto-legendre points of integration
double precision :: xigll(ngllx)
double precision :: etagll(nglly)

! derivatives of the 3d shape functions
double precision :: dshape_quad4(2,ngnod2d,ngllx*nglly)

integer :: i,j,i_gnod
integer :: igll,ngll

! location of the nodes of the 2d quadrilateral elements
!double precision :: xi,eta
double precision :: xip,xim,etap,etam

! for checking the 3d shape functions
double precision :: sum_dshapexi,sum_dshapeeta

double precision, parameter :: one_fourth = 0.25_kreal

! check that the parameter file is correct
if(ngnod2d /= 4)then
  write(*,*)'ERROR: element faces must have 4 geometrical nodes!'
  stop
endif

ngll=ngllx*nglly

! compute the derivatives of 2d shape functions
igll=0
do j=1,nglly
  etap = one + etagll(j)
  etam = one - etagll(j)
  do i=1,ngllx
    igll=igll+1

    xip = one + xigll(i)
    xim = one - xigll(i)

    ! corner nodes
    !shape_quad4(1,igll) = one_fourth*xim*etam
    !shape_quad4(2,igll) = one_fourth*xip*etam
    !shape_quad4(3,igll) = one_fourth*xip*etap
    !shape_quad4(4,igll) = one_fourth*xim*etap

    dshape_quad4(1,1,igll) = -one_fourth*etam
    dshape_quad4(1,2,igll) = one_fourth*etam
    dshape_quad4(1,3,igll) = one_fourth*etap
    dshape_quad4(1,4,igll) = -one_fourth*etap

    dshape_quad4(2,1,igll) = -one_fourth*xim
    dshape_quad4(2,2,igll) = -one_fourth*xip
    dshape_quad4(2,3,igll) = one_fourth*xip
    dshape_quad4(2,4,igll) = one_fourth*xim

    enddo
  enddo

  ! check the shape functions and their derivatives
  do i=1,ngll
    sum_dshapexi = zero
    sum_dshapeeta = zero

    do i_gnod=1,ngnod2d
      sum_dshapexi = sum_dshapexi + dshape_quad4(1,i_gnod,i)
      sum_dshapeeta = sum_dshapeeta + dshape_quad4(2,i_gnod,i)
    enddo

    ! sum of derivative of shape functions should be zero
    if(abs(sum_dshapexi) >  zerotol)then
      write(*,*)'ERROR: derivative xi shape functions!'
      stop
    endif
    if(abs(sum_dshapeeta) >  zerotol)then
      write(*,*)'ERROR: derivative eta shape functions!'
      stop
    endif
  enddo

end subroutine dshape_function_quad4
!===============================================================================

! This subroutines computes the shape fucntions at a given &
! point. The 4-noded quadrilateral is conformed to the exodus/cubit numbering
! convention
subroutine shape_function_quad4p(ngnod2d,xi,eta,shape_quad4)
use set_precision
use math_constants
implicit none
integer,intent(in) :: ngnod2d

! Given point
double precision :: xi
double precision :: eta

! The 3d shape functions
double precision :: shape_quad4(ngnod2d)

double precision :: xip,xim,etap,etam

double precision, parameter :: one_fourth = 0.25_kreal

! check that the parameter file is correct
if(ngnod2d /= 4)then
  write(*,*)'ERROR: element faces must have 4 geometrical nodes!'
  stop
endif

! compute the derivatives of 2d shape functions
etap = one + eta
etam = one - eta

xip = one + xi
xim = one - xi

shape_quad4(1) = one_fourth*xim*etam
shape_quad4(2) = one_fourth*xip*etam
shape_quad4(3) = one_fourth*xip*etap
shape_quad4(4) = one_fourth*xim*etap

! sum of shape functions should be one
if(abs(sum(shape_quad4)-one)  >  zerotol)then
  write(*,*)'ERROR: sum of shape functions!'
  stop
endif

end subroutine shape_function_quad4p
!===============================================================================

! This subroutines computes derivatives of the shape fucntions at a given &
! point. The 4-noded quadrilateral is conformed to the exodus/cubit numbering
! convention
subroutine dshape_function_quad4p(ngnod2d,xi,eta,dshape_quad4)
use set_precision
use math_constants
implicit none
integer,intent(in) :: ngnod2d

! Given point
double precision :: xi
double precision :: eta

! derivatives of the 3d shape functions
double precision :: dshape_quad4(ngnod2d,2)

integer :: i_gnod
! location of the nodes of the 2d quadrilateral elements
!double precision :: xi,eta
double precision :: xip,xim,etap,etam

! for checking the 3d shape functions
double precision :: sum_dshapexi,sum_dshapeeta

double precision, parameter :: one_fourth = 0.25_kreal

! check that the parameter file is correct
if(ngnod2d /= 4)then
  write(*,*)'ERROR: element faces must have 4 geometrical nodes!'
  stop
endif

! compute the derivatives of 2d shape functions
etap = one + eta
etam = one - eta

xip = one + xi
xim = one - xi

! corner nodes
!shape_quad4(1,igll) = one_fourth*xim*etam
!shape_quad4(2,igll) = one_fourth*xip*etam
!shape_quad4(3,igll) = one_fourth*xip*etap
!shape_quad4(4,igll) = one_fourth*xim*etap

dshape_quad4(1,1) = -one_fourth*etam
dshape_quad4(2,1) = one_fourth*etam
dshape_quad4(3,1) = one_fourth*etap
dshape_quad4(4,1) = -one_fourth*etap
                
dshape_quad4(1,2) = -one_fourth*xim
dshape_quad4(2,2) = -one_fourth*xip
dshape_quad4(3,2) = one_fourth*xip
dshape_quad4(4,2) = one_fourth*xim

! check the shape functions and their derivatives
sum_dshapexi = zero
sum_dshapeeta = zero

do i_gnod=1,ngnod2d
  sum_dshapexi = sum_dshapexi + dshape_quad4(i_gnod,1)
  sum_dshapeeta = sum_dshapeeta + dshape_quad4(i_gnod,2)
enddo

! sum of derivative of shape functions should be zero
if(abs(sum_dshapexi) >  zerotol)then
  write(*,*)'ERROR: derivative xi shape functions!'
  stop
endif
if(abs(sum_dshapeeta) >  zerotol)then
  write(*,*)'ERROR: derivative eta shape functions!'
  stop
endif

end subroutine dshape_function_quad4p
!===============================================================================

end module shape_library
!===============================================================================
