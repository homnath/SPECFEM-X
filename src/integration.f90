! AUTHOR
!   Hom Nath Gharti
module integration
use set_precision
implicit none
character(len=250),private :: myfname=' => integration.f90'
character(len=500),private :: errsrc

! 3D
real(kind=kreal),allocatable :: gll_weights(:)
real(kind=kreal),allocatable :: lagrange_gll(:,:),dlagrange_gll(:,:,:)
real(kind=kreal),allocatable :: dshape_hex8(:,:,:)

! 2D
real(kind=kreal),dimension(:,:,:),allocatable :: dshape_quad4_xy,              &
dshape_quad4_yz,dshape_quad4_zx
real(kind=kreal),dimension(:),allocatable :: gll_weights_xy,                   &
gll_weights_yz,gll_weights_zx
real(kind=kreal),dimension(:,:),allocatable :: lagrange_gll_xy,                &
lagrange_gll_yz,lagrange_gll_zx
real(kind=kreal),dimension(:,:,:),allocatable :: dlagrange_gll_xy,             &
dlagrange_gll_yz,dlagrange_gll_zx

contains

!-------------------------------------------------------------------------------
subroutine prepare_integration(errcode,errtag)
use global,only:ndim,ngllx,nglly,ngllz,ngll,ngnode
use gll_library,only:gllpx,gllpy,gllpz,gll_quadrature
use shape_library,only:dshape_function_hex8
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

errtag="ERROR: unknown!"
errcode=-1

! get derivatives of shape functions for 8-noded hex
allocate(dshape_hex8(ndim,ngnode,ngll))
allocate(gll_weights(ngll),lagrange_gll(ngll,ngll),                            &
         dlagrange_gll(ndim,ngll,ngll))

call dshape_function_hex8(ngnode,ngllx,nglly,ngllz,gllpx,gllpy,gllpz,  &
dshape_hex8)

! compute gauss-lobatto-legendre quadrature information
call gll_quadrature(ndim,ngllx,nglly,ngllz,ngll,gll_weights,        &
lagrange_gll,dlagrange_gll)

errcode=0
return

end subroutine prepare_integration
!===============================================================================

subroutine cleanup_integration(errcode,errtag)
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

errtag="ERROR: unknown!"
errcode=-1

deallocate(dshape_hex8)
deallocate(gll_weights)
deallocate(lagrange_gll)
deallocate(dlagrange_gll)

errcode=0
return

end subroutine cleanup_integration
!===============================================================================

subroutine prepare_integration2d(errcode,errtag)
use global,only:ngllx,nglly,ngllz,ngllxy,ngllyz,ngllzx
use gll_library,only:gllpx,gllpy,gllpz,gll_quadrature2d,zwgljd
use shape_library,only:dshape_function_quad4
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal
!real(kind=kreal),dimension(ngllx) :: xigll,wxgll !double precision
!real(kind=kreal),dimension(nglly) :: etagll,wygll !double precision
!real(kind=kreal),dimension(ngllz) :: zetagll,wzgll !double precision
real(kind=kreal),dimension(:,:),allocatable :: gll_points_xy,gll_points_yz,    &
gll_points_zx

errtag="ERROR: unknown!"
errcode=-1

!! compute GLL points and weights
!call zwgljd(xigll,wxgll,ngllx,jacobi_alpha,jacobi_beta)
!call zwgljd(etagll,wygll,nglly,jacobi_alpha,jacobi_beta)
!call zwgljd(zetagll,wzgll,ngllz,jacobi_alpha,jacobi_beta)

allocate(dshape_quad4_xy(2,4,ngllxy),dshape_quad4_yz(2,4,ngllyz),              &
dshape_quad4_zx(2,4,ngllzx))
allocate(gll_weights_xy(ngllxy),gll_points_xy(2,ngllxy),                       &
lagrange_gll_xy(ngllxy,ngllxy),dlagrange_gll_xy(2,ngllxy,ngllxy))
allocate(gll_weights_yz(ngllyz),gll_points_yz(2,ngllyz),                       &
lagrange_gll_yz(ngllyz,ngllyz),dlagrange_gll_yz(2,ngllyz,ngllyz))
allocate(gll_weights_zx(ngllzx),gll_points_zx(2,ngllzx),                       &
lagrange_gll_zx(ngllzx,ngllzx),dlagrange_gll_zx(2,ngllzx,ngllzx))

call dshape_function_quad4(4,ngllx,nglly,gllpx,gllpy,dshape_quad4_xy)
call dshape_function_quad4(4,nglly,ngllz,gllpy,gllpz,dshape_quad4_yz)
call dshape_function_quad4(4,ngllz,ngllx,gllpz,gllpx,dshape_quad4_zx)

call gll_quadrature2d(2,ngllx,nglly,ngllxy,gll_points_xy,gll_weights_xy,       &
lagrange_gll_xy,dlagrange_gll_xy)
call gll_quadrature2d(2,nglly,ngllz,ngllyz,gll_points_yz,gll_weights_yz,       &
lagrange_gll_yz,dlagrange_gll_yz)
call gll_quadrature2d(2,ngllz,ngllx,ngllzx,gll_points_zx,gll_weights_zx,       &
lagrange_gll_zx,dlagrange_gll_zx)

deallocate(gll_points_xy,gll_points_yz,gll_points_zx)

errcode=0
return

end subroutine prepare_integration2d
!===============================================================================

subroutine cleanup_integration2d(errcode,errtag)
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

errtag="ERROR: unknown!"
errcode=-1

if(allocated(dshape_quad4_xy))deallocate(dshape_quad4_xy)
if(allocated(dshape_quad4_yz))deallocate(dshape_quad4_yz)
if(allocated(dshape_quad4_zx))deallocate(dshape_quad4_zx)
if(allocated(dshape_quad4_xy))deallocate(dshape_quad4_xy)
if(allocated(dshape_quad4_yz))deallocate(dshape_quad4_yz)
if(allocated(dshape_quad4_zx))deallocate(dshape_quad4_zx)
if(allocated(gll_weights_xy))deallocate(gll_weights_xy)
if(allocated(lagrange_gll_xy))deallocate(lagrange_gll_xy)
if(allocated(dlagrange_gll_xy))deallocate(dlagrange_gll_xy)
if(allocated(gll_weights_yz))deallocate(gll_weights_yz)
if(allocated(lagrange_gll_yz))deallocate(lagrange_gll_yz)
if(allocated(dlagrange_gll_yz))deallocate(dlagrange_gll_yz)
if(allocated(gll_weights_zx))deallocate(gll_weights_zx)
if(allocated(lagrange_gll_zx))deallocate(lagrange_gll_zx)
if(allocated(dlagrange_gll_zx))deallocate(dlagrange_gll_zx)

errcode=0
return

end subroutine cleanup_integration2d
!===============================================================================

end module integration
!===============================================================================
