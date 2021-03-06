! This module contains preprocessing library routines.
! REVISION:
!   HNG, Jul 07,2011
module preprocess
implicit none
character(len=250),private :: myfname=" => preprocess.f90"
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

! This subrotine computes the stiffness matrix, and
! body loads contributed by mass density or magnetiztion.
! TODO: optional precoditioner,optional assembly of stiffness
subroutine compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
use set_precision
use global,only:myrank,NDIM,nst,nelmt,ngll,nedof,nedofu,nedofphi,nenode,ngnode,&
ngllx,nglly,ngllz,ngll,g_coord,gdof_elmt,g_num,mat_domain,mat_id,massdens_elmt,&
bulkmod_elmt,shearmod_elmt,isempty_blk,rho_blk,ym_blk,magnetization_elmt,&
infinite_iface,infinite_face_idir,pole_coord0,pole_coord1,       &
pole_type,pole_axis,axis_range,ISDISP_DOF,ISPOT_DOF,storederiv,    &
POT_TYPE,PGRAVITY,PMAGNETIC,    &
storejw,devel_nondim, &
isdxval,isdyval,isdzval,devel_gaminf,infquad, &
edofu,edofphi,grav0_nodal,dgrav0_elmt,ISGRAV0,&
imat_to_imatmag,magnetization_blk,ismagnet_blk
use element,only:hex8_gnode,map2exodus_hex8
use elastic,only:compute_cmat_elastic
use math_constants,only:HALF,ONE,ZERO,FOUR,GRAV_CONS,PI
use math_library,only:determinant,invert,issymmetric
use weakform
use shape_library
use gll_library
use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights,       &
prepare_integration
use infinite_element
!use ieee_arithmetic
implicit none
real(kind=kreal),intent(out) :: storekmat(:,:,:)
real(kind=kreal),intent(out) :: rhoload(0:)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
real(kind=kreal),parameter :: FOUR_PI_G=FOUR*PI*GRAV_CONS
real(kind=kreal),parameter :: FOUR_PI_G_INV=ONE/FOUR_PI_G
integer :: i,i_gll
integer :: i_elmt,ielmt,imat
integer :: imatve,mdomain
integer :: num(nenode),egdof(nedof)
real(kind=kreal) :: cmat(nst,nst)
real(kind=kreal) :: detjac !determinant of Jacobian
real(kind=kreal) :: coord(ngnode,NDIM),jac(NDIM,NDIM),eload(nedof)
real(kind=kreal) :: dgmat(NDIM,NDIM),eg0(ngll,NDIM),g0(NDIM),dg0(6)

real(kind=kreal) :: kmat(nedof,nedof),xval
real(kind=kreal) :: interpf(NGLL),deriv(NDIM,nenode)

! magnetization
integer :: imatmag
real(kind=kreal) :: M(NDIM),Mgll(NDIM,ngll)
real(kind=kreal) :: divM,dxMx,dyMy,dzMz

real(kind=kreal) :: bmatu(NST,NEDOFU),wmat_gradphi(NEDOFU,NDIM),      &
rmat_gradphi(NDIM,NEDOFPHI),wmat_sPE(NEDOFPHI,NDIM),rmat_sPE(NDIM,NEDOFU), &
rmat_term1(NDIM,NEDOFU),wmat_term1(NEDOFU,NDIM),                           &
rmat_term2(NDIM,NEDOFU),wmat_term2(NEDOFU,NDIM),                           &
rmat_term3(1,NEDOFU),wmat_term3(NEDOFU,1),                                 &
rmat_term4(1,NEDOFU),wmat_term4(NEDOFU,1)

real(kind=kreal) :: tratio

! for infinite elements
integer,parameter :: nginf=8
integer :: i_face,idir
real(kind=kreal) :: tcoord
real(kind=kreal) :: gaminf
real(kind=kreal) :: polex(4,NDIM)
real(kind=kreal) :: coordinf(nginf,NDIM)
! jacw=jacobian*weight
real(kind=kreal) :: jacw
real(kind=kreal),allocatable :: shape_infinite(:,:),dshape_infinite(:,:,:)
real(kind=kreal),allocatable :: lagrange_gl(:,:),dlagrange_gl(:,:,:)
real(kind=kreal),allocatable :: GLw(:)
integer :: nip,nipinf
logical :: isinf ! flag to check if the element if on the infinite domain
logical :: isfaces(6)

errtag="ERROR: unknown!"
errcode=-1
errsrc=trim(myfname)//' => compute_stiffness_elastic'

! use ng for Gauss quadrature and ngll fro GLL-Radau quadrature
gaminf=2.00_kreal !1.99_kreal ! 2.0: X1 in the mid-position. gaminf must be > 1.0
if(devel_gaminf.gt.ONE .and. devel_gaminf.lt.FOUR)gaminf=devel_gaminf
nipinf=ngll

allocate(shape_infinite(nipinf,nginf),dshape_infinite(NDIM,nipinf,nginf))
allocate(lagrange_gl(nipinf,ngll),dlagrange_gl(NDIM,nipinf,ngll))

allocate(GLw(nipinf))

storekmat=zero
rhoload=zero
storederiv=zero
storejw=zero
! Purely elastic elements
! Viscoelastic elements are elastic at time = 0
! Following loops through nelmt_elas+nelmt_viscoelas
do i_elmt=1,nelmt
  ielmt=i_elmt
  num=g_num(:,ielmt)
  imat=mat_id(ielmt)
  mdomain=mat_domain(imat) 
  coord=transpose(g_coord(:,num(hex8_gnode)))
  nip=ngll

  ! set magnetization at GLL points
  Mgll=ZERO
  if(POT_TYPE==PMAGNETIC)then
    if(ismagnet_blk(imat))then
      ! ONLY FOR BLOCK PROPERTIES
      !imatmag=imat_to_imatmag(imat)
      !M=magnetization_blk(:,imatmag)
      !do i_gll=1,ngll 
      !  Mgll(:,i_gll)=M
      !enddo
      do i_gll=1,ngll 
        Mgll(:,i_gll)=magnetization_elmt(:,i_gll,ielmt)
      enddo
    endif
  endif
  !H=dgrav0_elmt(:,:,ielmt)
  eg0=transpose(grav0_nodal(:,num))

  isfaces=infinite_iface(:,ielmt)
  isinf=any(isfaces)
  if(count(isfaces).gt.1)isinf=.false.
  ! set coordinates
  if(isinf)then
    !coordinf=transpose(g_coord(:,num(gnodinf)))
    coordinf=transpose(g_coord(:,num(hex8_gnode)))

    ! Loop through infinite faces
    do i_face=1,6
      if(.not.isfaces(i_face))cycle
      idir=infinite_face_idir(i_face,ielmt)
      if(i_face==1)then
        ! ymin face
        ! Set X coordinate of the pole
        if(trim(pole_type)=='plane')then
          polex(1,:)=coordinf(3,:)!coordinf(1,:)
          polex(2,:)=coordinf(4,:)!coordinf(4,:)
          polex(3,:)=coordinf(7,:)!coordinf(5,:)
          polex(4,:)=coordinf(8,:)!coordinf(8,:)
          polex(:,idir)=pole_coord0(idir)
        elseif(trim(pole_type)=='axis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          polex(1,pole_axis)=coordinf(3,pole_axis)
          polex(2,pole_axis)=coordinf(4,pole_axis)
          polex(3,pole_axis)=coordinf(7,pole_axis)
          polex(4,pole_axis)=coordinf(8,pole_axis)
        elseif(trim(pole_type)=='pointaxis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          tcoord=coordinf(3,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(1,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(1,:)=pole_coord1
          endif
          tcoord=coordinf(4,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(2,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(2,:)=pole_coord1
          endif
          tcoord=coordinf(7,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(3,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(3,:)=pole_coord1
          endif
          tcoord=coordinf(8,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(4,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(4,:)=pole_coord1
          endif
        else
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
        endif
        ! X2
        coordinf(2,:)=polex(1,:)+gaminf*(coordinf(3,:)-polex(1,:))
        coordinf(1,:)=polex(2,:)+gaminf*(coordinf(4,:)-polex(2,:))
        coordinf(6,:)=polex(3,:)+gaminf*(coordinf(7,:)-polex(3,:))
        coordinf(5,:)=polex(4,:)+gaminf*(coordinf(8,:)-polex(4,:))
        ! Set near face coordinates to the pole coordinates
        ! X0
        coordinf(3,:)=polex(1,:)
        coordinf(4,:)=polex(2,:)
        coordinf(7,:)=polex(3,:)
        coordinf(8,:)=polex(4,:)
      
      elseif(i_face==2)then
        ! xmax face
        ! Set X coordinate of the pole
        if(trim(pole_type)=='plane')then
          polex(1,:)=coordinf(1,:)!coordinf(1,:)
          polex(2,:)=coordinf(4,:)!coordinf(4,:)
          polex(3,:)=coordinf(5,:)!coordinf(5,:)
          polex(4,:)=coordinf(8,:)!coordinf(8,:)
          polex(:,idir)=pole_coord0(idir)
        elseif(trim(pole_type)=='axis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          polex(1,pole_axis)=coordinf(1,pole_axis)
          polex(2,pole_axis)=coordinf(4,pole_axis)
          polex(3,pole_axis)=coordinf(5,pole_axis)
          polex(4,pole_axis)=coordinf(8,pole_axis)
        elseif(trim(pole_type)=='pointaxis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          tcoord=coordinf(1,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(1,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(1,:)=pole_coord1
          endif
          tcoord=coordinf(4,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(2,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(2,:)=pole_coord1
          endif
          tcoord=coordinf(5,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(3,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(3,:)=pole_coord1
          endif
          tcoord=coordinf(8,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(4,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(4,:)=pole_coord1
          endif
        else
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
        endif
        ! X2
        coordinf(2,:)=polex(1,:)+gaminf*(coordinf(1,:)-polex(1,:))
        coordinf(3,:)=polex(2,:)+gaminf*(coordinf(4,:)-polex(2,:))
        coordinf(6,:)=polex(3,:)+gaminf*(coordinf(5,:)-polex(3,:))
        coordinf(7,:)=polex(4,:)+gaminf*(coordinf(8,:)-polex(4,:))
        ! Set near face coordinates to the pole coordinates
        ! X0
        coordinf(1,:)=polex(1,:)
        coordinf(4,:)=polex(2,:)
        coordinf(5,:)=polex(3,:)
        coordinf(8,:)=polex(4,:)
      
      elseif(i_face==3)then
        ! ymax face
        ! Set Y coordinate of the pole 
        if(trim(pole_type)=='plane')then
          polex(1,:)=coordinf(1,:)
          polex(2,:)=coordinf(2,:)
          polex(3,:)=coordinf(5,:)
          polex(4,:)=coordinf(6,:)
          polex(:,idir)=pole_coord0(idir)
        elseif(trim(pole_type)=='axis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          polex(1,pole_axis)=coordinf(1,pole_axis)
          polex(2,pole_axis)=coordinf(2,pole_axis)
          polex(3,pole_axis)=coordinf(5,pole_axis)
          polex(4,pole_axis)=coordinf(6,pole_axis)
        elseif(trim(pole_type)=='pointaxis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          tcoord=coordinf(1,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(1,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(1,:)=pole_coord1
          endif
          tcoord=coordinf(2,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(2,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(2,:)=pole_coord1
          endif
          tcoord=coordinf(5,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(3,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(3,:)=pole_coord1
          endif
          tcoord=coordinf(6,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(4,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(4,:)=pole_coord1
          endif
        else
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
        endif
        ! X2
        coordinf(4,:)=polex(1,:)+gaminf*(coordinf(1,:)-polex(1,:))
        coordinf(3,:)=polex(2,:)+gaminf*(coordinf(2,:)-polex(2,:))
        coordinf(8,:)=polex(3,:)+gaminf*(coordinf(5,:)-polex(3,:))
        coordinf(7,:)=polex(4,:)+gaminf*(coordinf(6,:)-polex(4,:))
        ! Set near face coordinates to the pole coordinates
        ! X0
        coordinf(1,:)=polex(1,:)
        coordinf(2,:)=polex(2,:)
        coordinf(5,:)=polex(3,:)
        coordinf(6,:)=polex(4,:)

      elseif(i_face==4)then
        ! xmin face
        ! Set Y coordinate of the pole 
        if(trim(pole_type)=='plane')then
          polex(1,:)=coordinf(2,:)
          polex(2,:)=coordinf(3,:)
          polex(3,:)=coordinf(6,:)
          polex(4,:)=coordinf(7,:)
          polex(:,idir)=pole_coord0(idir)
        elseif(trim(pole_type)=='axis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          polex(1,pole_axis)=coordinf(2,pole_axis)
          polex(2,pole_axis)=coordinf(3,pole_axis)
          polex(3,pole_axis)=coordinf(6,pole_axis)
          polex(4,pole_axis)=coordinf(7,pole_axis)
        elseif(trim(pole_type)=='pointaxis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          tcoord=coordinf(2,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(1,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(1,:)=pole_coord1
          endif
          tcoord=coordinf(3,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(2,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(2,:)=pole_coord1
          endif
          tcoord=coordinf(6,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(3,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(3,:)=pole_coord1
          endif
          tcoord=coordinf(7,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(4,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(4,:)=pole_coord1
          endif
        else
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
        endif
        ! X2
        coordinf(1,:)=polex(1,:)+gaminf*(coordinf(2,:)-polex(1,:))
        coordinf(4,:)=polex(2,:)+gaminf*(coordinf(3,:)-polex(2,:))
        coordinf(5,:)=polex(3,:)+gaminf*(coordinf(6,:)-polex(3,:))
        coordinf(8,:)=polex(4,:)+gaminf*(coordinf(7,:)-polex(4,:))
        ! Set near face coordinates to the pole coordinates
        ! X0
        coordinf(2,:)=polex(1,:)
        coordinf(3,:)=polex(2,:)
        coordinf(6,:)=polex(3,:)
        coordinf(7,:)=polex(4,:)

      elseif(i_face==5)then
        ! zmin face
        ! Set Y coordinate of the pole 
        if(trim(pole_type)=='plane')then
          polex(1,:)=coordinf(5,:)
          polex(2,:)=coordinf(6,:)
          polex(3,:)=coordinf(8,:)
          polex(4,:)=coordinf(7,:)
          polex(:,idir)=pole_coord0(idir)
        elseif(trim(pole_type)=='axis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          polex(1,pole_axis)=coordinf(5,pole_axis)
          polex(2,pole_axis)=coordinf(6,pole_axis)
          polex(3,pole_axis)=coordinf(8,pole_axis)
          polex(4,pole_axis)=coordinf(7,pole_axis)
        elseif(trim(pole_type)=='pointaxis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          tcoord=coordinf(5,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(1,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(1,:)=pole_coord1
          endif
          tcoord=coordinf(6,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(2,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(2,:)=pole_coord1
          endif
          tcoord=coordinf(8,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(3,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(3,:)=pole_coord1
          endif
          tcoord=coordinf(7,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(4,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(4,:)=pole_coord1
          endif
        else
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
        endif
        ! X2
        coordinf(1,:)=polex(1,:)+gaminf*(coordinf(5,:)-polex(1,:))
        coordinf(2,:)=polex(2,:)+gaminf*(coordinf(6,:)-polex(2,:))
        coordinf(4,:)=polex(3,:)+gaminf*(coordinf(8,:)-polex(3,:))
        coordinf(3,:)=polex(4,:)+gaminf*(coordinf(7,:)-polex(4,:))
        ! Set near face coordinates to the pole coordinates
        ! X0
        coordinf(5,:)=polex(1,:)
        coordinf(6,:)=polex(2,:)
        coordinf(8,:)=polex(3,:)
        coordinf(7,:)=polex(4,:)

      elseif(i_face==6)then
        ! zmax face
        ! Set Z coordinate of the pole 
        if(trim(pole_type)=='plane')then
          polex(1,:)=coordinf(1,:)
          polex(2,:)=coordinf(2,:)
          polex(3,:)=coordinf(3,:)
          polex(4,:)=coordinf(4,:)
          polex(:,idir)=pole_coord0(idir)
        elseif(trim(pole_type)=='axis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          polex(1,pole_axis)=coordinf(1,pole_axis)
          polex(2,pole_axis)=coordinf(2,pole_axis)
          polex(3,pole_axis)=coordinf(3,pole_axis)
          polex(4,pole_axis)=coordinf(4,pole_axis)
        elseif(trim(pole_type)=='pointaxis')then
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
          tcoord=coordinf(1,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(1,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(1,:)=pole_coord1
          endif
          tcoord=coordinf(2,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(2,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(2,:)=pole_coord1
          endif
          tcoord=coordinf(3,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(3,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(3,:)=pole_coord1
          endif
          tcoord=coordinf(4,pole_axis)
          if(tcoord.gt.axis_range(1).and.tcoord.lt.axis_range(2))then
            polex(4,pole_axis)=tcoord
          endif
          if(tcoord.ge.axis_range(2))then
            polex(4,:)=pole_coord1
          endif
        else
          polex(1,:)=pole_coord0
          polex(2,:)=pole_coord0
          polex(3,:)=pole_coord0
          polex(4,:)=pole_coord0
        endif
        ! X2
        coordinf(5,:)=polex(1,:)+gaminf*(coordinf(1,:)-polex(1,:))
        coordinf(6,:)=polex(2,:)+gaminf*(coordinf(2,:)-polex(2,:))
        coordinf(7,:)=polex(3,:)+gaminf*(coordinf(3,:)-polex(3,:))
        coordinf(8,:)=polex(4,:)+gaminf*(coordinf(4,:)-polex(4,:))
        ! Set near face coordinates to the pole coordinates
        ! X0
        coordinf(1,:)=polex(1,:)
        coordinf(2,:)=polex(2,:)
        coordinf(3,:)=polex(3,:)
        coordinf(4,:)=polex(4,:)

      endif

    enddo ! i_face

    ! convert gnode ordering to indicial order to match the infinite shape
    ! functions ordering
    ! if we use map2exodus_hex8 on exodus order it becomes indicial order
    coordinf=coordinf(map2exodus_hex8,:)
    nip=nipinf !ngll ! we use GLL or GLL-Radau quadrature
    
    ! Radau or Gauss quadrature
    call shape_function_infiniteGLHEX8ZW(infquad,ngllx,nglly,ngllz,    &
    ngll,nip,isfaces,shape_infinite,dshape_infinite,lagrange_gl,       &
    dlagrange_gl,GLw)
    
  endif

  egdof=gdof_elmt(:,i_elmt)
    
  kmat=zero
  eload=zero
  do i=1,nip
    if(ISDISP_DOF)then
      call compute_cmat_elastic(bulkmod_elmt(i,ielmt),shearmod_elmt(i,ielmt),  &
      cmat)
    endif

    if(isinf)then
      ! infinite element
      interpf=lagrange_gl(i,:)
      
      jac=matmul(dshape_infinite(:,i,:),coordinf)
      detjac=determinant(jac)
      if(detjac.le.zero.and.myrank==0)then
        write(*,*)'ERROR: zero or negative jacobian in infinite element!'
        write(*,*)'HINT: check "pole_type" and "infquad"!'
        write(*,*)'HINT: make sure that coordinates units are consistent!'
        write(*,*)myrank,i_elmt,i,nip,detjac
        write(*,*)isfaces
        write(*,*)infquad
        write(*,*)coordinf
        stop
      endif
      call invert(jac)
      deriv=matmul(jac,dlagrange_gl(:,i,:))

      ! set derivative constraint
      if(isdxval)deriv(1,:)=ZERO
      if(isdyval)deriv(2,:)=ZERO
      if(isdzval)deriv(3,:)=ZERO
      
      storederiv(:,:,i,ielmt)=deriv 
      jacw=detjac*GLw(i)
      storejw(i,ielmt)=jacw
      
      if(ISDISP_DOF)then
        ! compute only for nonempty elements
        if(.not.isempty_blk(imat))then
          call compute_bmat_stress(deriv,bmatu)
          kmat(edofu,edofu)=kmat(edofu,edofu)+matmul(matmul(transpose(bmatu),cmat),bmatu)*jacw

          if(ISGRAV0)then
            g0=eg0(i,:)
            ! compute dg0=\nabla g0
            dgmat=matmul(deriv,eg0)
            dg0(1)=dgmat(1,1)
            dg0(2)=dgmat(2,2)
            dg0(3)=dgmat(3,3)
            dg0(4)=dgmat(1,2)
            dg0(5)=dgmat(1,3)
            dg0(6)=dgmat(2,3)

            call compute_rmat_term1(massdens_elmt(i,ielmt),interpf,rmat_term1)
            call compute_wmat_term1(g0,dg0,interpf,deriv,wmat_term1)
            call compute_rmat_term2(g0,dg0,interpf,deriv,rmat_term2)
            call compute_wmat_term2(massdens_elmt(i,ielmt),interpf,wmat_term2)
            kmat(edofu,edofu)=kmat(edofu,edofu)-HALF*(matmul(wmat_term1,rmat_term1)+ &
            matmul(wmat_term2,rmat_term2))*jacw
  
            call compute_rmat_term3(massdens_elmt(i,ielmt),g0,interpf,rmat_term3)
            call compute_wmat_term3(deriv,wmat_term3)
            call compute_rmat_term4(massdens_elmt(i,ielmt),deriv,rmat_term4)
            call compute_wmat_term4(g0,interpf,wmat_term4)
            kmat(edofu,edofu)=kmat(edofu,edofu)+HALF*(matmul(wmat_term3,rmat_term3)+ &
            matmul(wmat_term4,rmat_term4))*jacw
          endif
          
          if(ISPOT_DOF)then
            ! w.rho*grad(phi)
            call compute_wmat_gradphi(interpf,wmat_gradphi)
            call compute_rmat_gradphi(massdens_elmt(i,ielmt),deriv,rmat_gradphi)          
            kmat(edofu,edofphi)=kmat(edofu,edofphi)+matmul(wmat_gradphi,rmat_gradphi)*jacw
            ! grad(w).rho*s
            ! this term is transpose of the previous. Therefore it will be
            ! computed later
            call compute_wmat_sPE(deriv,wmat_sPE)
            call compute_rmat_sPE(massdens_elmt(i,ielmt),interpf,rmat_sPE)
            kmat(edofphi,edofu)=kmat(edofphi,edofu)+matmul(wmat_sPE,rmat_sPE)*jacw
          endif
        endif
      endif

      if(ISPOT_DOF)then
        kmat(edofphi,edofphi)=kmat(edofphi,edofphi)+matmul(transpose(deriv),deriv)*jacw
        if(.not.ISDISP_DOF)then
          if(POT_TYPE==PMAGNETIC)then
            if(ismagnet_blk(imat))then
              ! first compute divegence of M
              ! \nabla.M=dxMx+dyMy+dzMz
              dxMx=dot_product(deriv(1,:),Mgll(1,:))
              dyMy=dot_product(deriv(2,:),Mgll(2,:))
              dzMz=dot_product(deriv(3,:),Mgll(3,:))
              divM=dxMx+dyMy+dzMz
              if(maxval(abs(M)).gt.ZERO)then
                eload(edofphi)=eload(edofphi)+lagrange_gl(i,:)*divM*jacw
              endif
            endif
          else
            eload(edofphi)=eload(edofphi)+lagrange_gl(i,:)*massdens_elmt(i,ielmt)*jacw
          endif
        endif
      endif

    else ! (isinf) 
      ! standard element
      interpf=lagrange_gll(i,:)
    
      jac=matmul(dshape_hex8(:,:,i),coord)
      detjac=determinant(jac)
      call invert(jac)
      deriv=matmul(jac,dlagrange_gll(:,i,:))
      ! set derivative constraint
      if(isdxval)deriv(1,:)=ZERO
      if(isdyval)deriv(2,:)=ZERO
      if(isdzval)deriv(3,:)=ZERO
      
      storederiv(:,:,i,ielmt)=deriv 
      jacw=detjac*gll_weights(i)
      storejw(i,ielmt)=jacw
     
      if(ISDISP_DOF)then
        ! compute only for nonempty elements
        if(.not.isempty_blk(imat))then
          call compute_bmat_stress(deriv,bmatu)
          kmat(edofu,edofu)=kmat(edofu,edofu)+matmul(matmul(transpose(bmatu),cmat),bmatu)*jacw
           
          if(ISGRAV0)then
            g0=eg0(i,:)
            ! compute dg0=\nabla g0
            dgmat=matmul(deriv,eg0)
            dg0(1)=dgmat(1,1)
            dg0(2)=dgmat(2,2)
            dg0(3)=dgmat(3,3)
            dg0(4)=dgmat(1,2)
            dg0(5)=dgmat(1,3)
            dg0(6)=dgmat(2,3)

            call compute_rmat_term1(massdens_elmt(i,ielmt),interpf,rmat_term1)
            call compute_wmat_term1(g0,dg0,interpf,deriv,wmat_term1)
            call compute_rmat_term2(g0,dg0,interpf,deriv,rmat_term2)
            call compute_wmat_term2(massdens_elmt(i,ielmt),interpf,wmat_term2)
            kmat(edofu,edofu)=kmat(edofu,edofu)-HALF*(matmul(wmat_term1,rmat_term1)+ &
            matmul(wmat_term2,rmat_term2))*jacw
  
            call compute_rmat_term3(massdens_elmt(i,ielmt),g0,interpf,rmat_term3)
            call compute_wmat_term3(deriv,wmat_term3)
            call compute_rmat_term4(massdens_elmt(i,ielmt),deriv,rmat_term4)
            call compute_wmat_term4(g0,interpf,wmat_term4)
            kmat(edofu,edofu)=kmat(edofu,edofu)+HALF*(matmul(wmat_term3,rmat_term3)+ &
            matmul(wmat_term4,rmat_term4))*jacw
          endif
          if(ISPOT_DOF)then
            ! w.rho*grad(phi)
            call compute_wmat_gradphi(interpf,wmat_gradphi)
            call compute_rmat_gradphi(massdens_elmt(i,ielmt),deriv,rmat_gradphi)          
            kmat(edofu,edofphi)=kmat(edofu,edofphi)+matmul(wmat_gradphi,rmat_gradphi)*jacw
            ! grad(w).rho*s
            ! This term is a transpose of the previous. Therefore, it will be
            ! computed later
            call compute_wmat_sPE(deriv,wmat_sPE)
            call compute_rmat_sPE(massdens_elmt(i,ielmt),interpf,rmat_sPE)
            kmat(edofphi,edofu)=kmat(edofphi,edofu)+matmul(wmat_sPE,rmat_sPE)*jacw
          endif
        endif
      endif

      if(ISPOT_DOF)then
        kmat(edofphi,edofphi)=kmat(edofphi,edofphi)+matmul(transpose(deriv),deriv)*jacw
        if(.not.ISDISP_DOF)then
          if(POT_TYPE==PMAGNETIC)then
            if(ismagnet_blk(imat))then
              ! first compute the divergence of M
              ! \nabla.M=dxMx+dyMy+dzMz
              dxMx=dot_product(deriv(1,:),Mgll(1,:))
              dyMy=dot_product(deriv(2,:),Mgll(2,:))
              dzMz=dot_product(deriv(3,:),Mgll(3,:))
              divM=dxMx+dyMy+dzMz
              if(maxval(abs(M)).gt.ZERO)then
                eload(edofphi)=eload(edofphi)+lagrange_gll(i,:)*divM*jacw
              endif
            endif
          else
            eload(edofphi)=eload(edofphi)+lagrange_gll(i,:)*massdens_elmt(i,ielmt)*jacw
          endif
        endif
      endif

    endif
  enddo
  ! kmat terms for grad(w).rho*s
  !kmat(edofphi,edofu)=transpose(kmat(edofu,edofphi))
  if(ISDISP_DOF.and.ISPOT_DOF)then
    if(devel_nondim)then
      ! Note: PI*G is nondimensionalized
      kmat(edofphi,edofphi)=0.25_kreal*kmat(edofphi,edofphi)
    else
      kmat(edofphi,edofphi)=FOUR_PI_G_INV*kmat(edofphi,edofphi)
    endif
  endif
  if(.not.issymmetric(kmat))then
    write(*,*)'ERROR: matrix is unsymmetric!'
    stop
  endif
  !do i=1,nedof
  !  do j=1,nedof
  !    xval=kmat(i,j)                                                  
  !    if(ieee_is_nan(xval).or. .not.ieee_is_finite(xval))then                    
  !      write(*,*)'ERROR: stiffness matrix has nonfinite value/s!',myrank,ielmt,isinf,&
  !      mat_id(ielmt),xval,minval(abs(kmat)),maxval(abs(kmat))         
  !      stop                                                                     
  !    endif     
  !  enddo
  !enddo
  storekmat(:,:,ielmt)=kmat
  if(.not.ISDISP_DOF .and. ISPOT_DOF)then
    rhoload(egdof)=rhoload(egdof)+eload
  endif
enddo ! i_elmt

! rhoload is computed if only the ISPOT_DOF is TRUE
! multiply rhoload by 4*PI*G
if(.not.ISDISP_DOF.and.ISPOT_DOF)then
  if(.not.devel_nondim)then
    if(POT_TYPE==PMAGNETIC)then
      !rhoload=rhoload
    else
      rhoload=FOUR_PI_G*rhoload
    endif
  else
    if(POT_TYPE==PMAGNETIC)then
      !rhoload=rhoload
    else
      rhoload=FOUR*rhoload
      ! Note: PI*G is nondimensionalized
    endif
  endif
  rhoload(0)=ZERO
endif
deallocate(shape_infinite,dshape_infinite)
deallocate(lagrange_gl,dlagrange_gl)
deallocate(GLw)

end subroutine compute_stiffness_elastic
!===============================================================================

! This subroutine computes the free surafce contribution
! on the stiffness matrix.
subroutine compute_surface_stiffness(storekmat,errcode,errtag)
use global
use math_constants
use element,only:hexface,hexface_sign
use integration,only:dshape_quad4_xy,dshape_quad4_yz,dshape_quad4_zx,          &
                     gll_weights_xy,gll_weights_yz,gll_weights_zx,             &
                     lagrange_gll_xy,lagrange_gll_yz,lagrange_gll_zx
use dof,only:set_face_vecdof,set_face_scaldof
use weakform,only:compute_rmat_sn
implicit none
real(kind=kreal),intent(inout) :: storekmat(:,:,:)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_face,i_gll,ios
integer :: ielmt,iface,nface,inode
integer :: num(nenode)
integer :: nfdofu,nfdofphi,nfdof,nfgll
!face nodal dof, face gll points
integer :: tractype,count_trac
logical :: trac_stat
real(kind=kreal) :: coord(NDIM,4)
real(kind=kreal) :: detjac
real(kind=kreal),dimension(NDIM) :: face_normal,dx_dxi,dx_deta

integer,allocatable :: edof(:),imapuf(:),imapphif(:)
real(kind=kreal),dimension(:,:,:),allocatable :: dshape_quad4
real(kind=kreal),dimension(:),allocatable :: gll_weights
real(kind=kreal),dimension(:,:),allocatable :: lagrange_gll
real(kind=kreal),dimension(:,:,:),allocatable :: dlagrange_gll
real(kind=kreal),allocatable :: interpf(:)                                       
real(kind=kreal),allocatable :: wmat_sn(:,:),rmat_sn(:,:),kmat(:,:),rhofgll(:)

character(len=80) :: fname
character(len=80) :: data_path

errtag="ERROR: unknown!"
errcode=-1
errsrc=trim(myfname)//' => compute_surface_stiffness'

! set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

fname=trim(data_path)//trim(trfile)//trim(ptail_inp)
open(unit=11,file=trim(fname),status='old',action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "',trim(fname),'" cannot be opened!'
  return
endif

allocate(dshape_quad4(2,4,maxngll2d))
allocate(gll_weights(maxngll2d),lagrange_gll(maxngll2d,maxngll2d),             &
dlagrange_gll(2,maxngll2d,maxngll2d))

trac_stat=.true. ! necessary for empty trfile
count_trac=0
traction: do
  read(11,*,iostat=ios)tractype
  if(ios/=0)exit traction
  count_trac=count_trac+1
  trac_stat=.false.

  read(11,*) ! skip this line, we do not need
  read(11,*)nface
  do i_face=1,nface
    read(11,*)ielmt,iface
    if(iface==1 .or. iface==3)then
      ! ZX faces
      nfgll=ngllzx
      lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_zx
      gll_weights(1:nfgll)=gll_weights_zx
      dshape_quad4(:,:,1:nfgll)=dshape_quad4_zx
    elseif(iface==2 .or. iface==4)then
      ! YZ faces
      nfgll=ngllyz
      lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_yz
      gll_weights(1:nfgll)=gll_weights_yz
      dshape_quad4(:,:,1:nfgll)=dshape_quad4_yz
    elseif(iface==5 .or. iface==6)then
      ! XY faces
      nfgll=ngllxy
      lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
      gll_weights(1:nfgll)=gll_weights_xy
      dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
    else
      write(errtag,'(a)')'ERROR: wrong face ID for traction!'
      exit traction
    endif
    nfdof=nfgll*nndof
    nfdofu=nfgll*NNDOFU                                                            
    nfdofphi=nfgll*NNDOFPHI                       
    allocate(edof(nfgll))
    allocate(rhofgll(nfgll))
    allocate(interpf(nfgll))
    allocate(kmat(nfdof,nfdof))                                                   
    allocate(imapuf(nfdofu),imapphif(nfdofphi))                                        
    allocate(wmat_sn(nfdofphi,1),rmat_sn(1,nfdofu))                                  
    call set_face_vecdof(nfgll,3,1,imapuf) !ndofphi=1
    call set_face_scaldof(nfgll,4,0,imapphif) !idofphi=4
  
    print*,'imapuf:',imapuf
    print*,'imapphif:',imapphif
    print*,'iface:',iface
    print*,'edof:',hexface(iface)%edof

    edof=hexface(iface)%edof
    num=g_num(:,ielmt)
    coord=g_coord(:,num(hexface(iface)%gnode))
    rhofgll=massdens_elmt(hexface(iface)%node,ielmt)
    kmat=zero
    ! compute numerical integration
    do i_gll=1,nfgll

      ! compute two vectors dx_dxi and dx_deta
      dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
      dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))

      ! Normal = (dx_dxi x dx_deta)
      face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
      face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
      face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

      detjac=sqrt(dot_product(face_normal,face_normal))
      face_normal=hexface_sign(iface)*face_normal/detjac
      print*,face_normal
      interpf=lagrange_gll(i_gll,:)
      wmat_sn(:,1)=interpf
      ! -w_\phi \rho s.n                                                                      
      call compute_rmat_sn(nfgll,face_normal,interpf,rmat_sn)                           
      kmat(imapphif,imapuf)=kmat(imapphif,imapuf)-matmul(wmat_sn,rmat_sn)* &
                            rhofgll(i_gll)*detjac*gll_weights(i_gll)
    enddo ! i_gll
    !kmat(imapuf,imapphif)=transpose(kmat(imapphif,imapuf))
    storekmat(edof,edof,ielmt)= storekmat(edof,edof,ielmt)+kmat!*1e0
    print*,minval(abs(kmat)),maxval(abs(kmat)),minval(rhofgll),maxval(rhofgll)
    deallocate(edof)                                                   
    deallocate(rhofgll)                                                   
    deallocate(interpf)
    deallocate(kmat)                                                   
    deallocate(imapuf,imapphif)                                        
    deallocate(wmat_sn,rmat_sn)                                  
  enddo ! i_face
  trac_stat=.true.
enddo traction

close(11)
deallocate(dshape_quad4)
deallocate(gll_weights,lagrange_gll,dlagrange_gll)
!deallocate(kmat,wmat_sn,rmat_sn)
!deallocate(imapuf,imapphif)
if(.not.trac_stat)then
  write(errtag,'(a)')'ERROR: all tractions cannot be read!'
  return
endif

errcode=0

return
end subroutine compute_surface_stiffness
!===============================================================================

! This subrotine computes the stiffness matrix for the viscoelastic elements.
! WARNING: it has to be modified for gravity perturbation.
subroutine compute_stiffness_viscoelastic(nelmt_viscoelas,eid_viscoelas,       &
           dt,relaxtime,storekmat,                              &
           errcode,errtag)
use set_precision
use global,only:myrank,NDIM,nst,ngll,nedof,nedofu,nedofphi,nenode,ngnode,      &
ngllx,nglly,ngllz,ngll,g_coord,gdof_elmt,g_num,mat_domain,mat_id,massdens_elmt,&
bulkmod_elmt,shearmod_elmt,rho_blk,ym_blk,imat_to_imatve, &
ISDISP_DOF,ISPOT_DOF,storederiv,        &
POT_TYPE,PGRAVITY,PMAGNETIC,    &
storejw,devel_nondim,isdxval,isdyval,isdzval, &
edofu,edofphi,grav0_nodal,dgrav0_elmt,ISGRAV0, &
muratio_blk,visco_model,VISCO_MAXWELL,VISCO_ZENER,VISCO_GENMAXWELL,nmaxwell
use element,only:hex8_gnode
use viscoelastic,only:compute_cmat_maxwell,compute_cmat_zener, &
compute_cmat_genmaxwell
use math_constants,only:HALF,ONE,ZERO,FOUR,GRAV_CONS,PI
use math_library,only:determinant,invert,issymmetric
use weakform
use shape_library
use gll_library
use infinite_element
implicit none
integer,intent(in) :: nelmt_viscoelas
integer,intent(in) :: eid_viscoelas(:)
real(kind=kreal),intent(in) :: dt,relaxtime(:,:)
real(kind=kreal),intent(inout) :: storekmat(:,:,:)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: i
integer :: i_elmt,ielmt,imat
integer :: imatve,mdomain
integer :: num(nenode)
real(kind=kreal) :: cmat(nst,nst)
real(kind=kreal) :: jacw !jacobian*weight
real(kind=kreal) :: eload(nedof)
real(kind=kreal) :: dgmat(NDIM,NDIM),eg0(ngll,NDIM),g0(NDIM),dg0(6)

real(kind=kreal) :: kmatu(nedofu,nedofu)

real(kind=kreal) :: interpf(NGLL),deriv(NDIM,nenode)
real(kind=kreal) :: Mvec(NDIM)
real(kind=kreal) :: bmatu(NST,NEDOFU),wmat_gradphi(NEDOFU,NDIM),      &
rmat_gradphi(NDIM,NEDOFPHI),wmat_sPE(NEDOFPHI,NDIM),rmat_sPE(NDIM,NEDOFU), &
rmat_term1(NDIM,NEDOFU),wmat_term1(NEDOFU,NDIM),                           &
rmat_term2(NDIM,NEDOFU),wmat_term2(NEDOFU,NDIM),                           &
rmat_term3(1,NEDOFU),wmat_term3(NEDOFU,1),                                 &
rmat_term4(1,NEDOFU),wmat_term4(NEDOFU,1)

real(kind=kreal) :: muratio(nmaxwell),tratio(nmaxwell)

errtag="ERROR: unknown!"
errcode=-1
errsrc=trim(myfname)//' => compute_stiffness_viscoelastic'

! Viscoelastic part of storekmat will be replaced. 
! kmat(edofu,edofu) is the only viscoelastic portion.
! kmat(edofphi,edofphi) and other offdiagonal terms must remain the same.

! If statement here isn't necessary since viscoelasticity matters only if 
! if there are displacement DOFs.
if(.not.ISDISP_DOF)return

! Viscoelastic elements
do i_elmt=1,nelmt_viscoelas
  ielmt=eid_viscoelas(i_elmt)
  imat=mat_id(ielmt)
  imatve=imat_to_imatve(imat)
  muratio=muratio_blk(:,imatve)
  ! For a more gneneral case tratio can also be variable within an element
  tratio=dt/relaxtime(:,imatve)
  
  kmatu=zero
  do i=1,ngll
    deriv=storederiv(:,:,i,ielmt)
    jacw=storejw(i,ielmt)
   
    call compute_cmat_genmaxwell(bulkmod_elmt(i,ielmt),shearmod_elmt(i,ielmt),&
    tratio,muratio,cmat)

    call compute_bmat_stress(deriv,bmatu)
    kmatu=kmatu+matmul(matmul(transpose(bmatu),cmat),bmatu)*jacw
  enddo !i
  storekmat(edofu,edofu,ielmt)=kmatu

enddo ! i_elmt

end subroutine compute_stiffness_viscoelastic
!===============================================================================

! This subrotine computes the stiffness matrix, diagonal preconditioner, and
! body loads (gravity and pseudostatic loads) optionally
! TODO: optional precoditioner,optional assembly of stiffness
subroutine stiffness_bodyload(nelmt,neq,gnod,g_num,gdof_elmt,mat_id,gam,       &
storkm,dprecon,extload,gravity,pseudoeq)
use set_precision
use global,only:NDIM,nst,ngll,nedof,nenode,ngnode,g_coord,eqkx,eqky,eqkz,      &
nmatblk,bulkmod_elmt,shearmod_elmt
use elastic,only:compute_cmat_elastic
use math_library,only:determinant,invert
use weakform,only:compute_bmat_stress
use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights
implicit none
integer,intent(in) :: nelmt,neq,gnod(8) ! nelmt (only intact elements)
integer,intent(in) :: g_num(nenode,nelmt),gdof_elmt(nedof,nelmt),mat_id(nelmt)
! only intact elements
real(kind=kreal),intent(in) :: gam(nmatblk)
real(kind=kreal),intent(out) :: storkm(nedof,nedof,nelmt),dprecon(0:neq)
real(kind=kreal),intent(inout),optional :: extload(0:neq)
logical,intent(in),optional :: gravity,pseudoeq

real(kind=kreal) :: detjac,zero=0.0_kreal
real(kind=kreal) :: cmat(nst,nst),coord(ngnode,NDIM),jac(NDIM,NDIM),           &
deriv(NDIM,nenode),bmatu(nst,nedof),eld(nedof),eqload(nedof),km(nedof,nedof)
integer :: egdof(nedof),num(nenode)
integer :: i,i_elmt,k

if(present(extload).and.(.not.present(gravity) .or. .not.present(pseudoeq)))then
  write(*,'(a)')'ERROR: both "gravity" and "pseudoeq" must be defined for &
  &"extload"!'
  stop
endif

! compute stiffness matrices
storkm=zero; dprecon=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))
  egdof=gdof_elmt(:,i_elmt)
  km=zero; eld=zero; eqload=zero
  do i=1,ngll
    call compute_cmat_elastic(bulkmod_elmt(i,i_elmt),shearmod_elmt(i,i_elmt),  &
    cmat)
    ! compute Jacobian at GLL point using 20 noded element
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)
    call invert(jac)

    deriv=matmul(jac,dlagrange_gll(:,i,:))
    call compute_bmat_stress(deriv,bmatu)
    km=km+matmul(matmul(transpose(bmatu),cmat),bmatu)*detjac*gll_weights(i)
    eld(3:nedof:3)=eld(3:nedof:3)+lagrange_gll(i,:)*detjac*gll_weights(i)
    !eld(2:nedof-1:3)=eld(2:nedof-1:3)+fun(:)*detjac*weights(i)
  enddo ! i=1,ngll
  storkm(:,:,i_elmt)=km
  do k=1,nedof
    dprecon(egdof(k))=dprecon(egdof(k))+km(k,k)
  enddo

  if(.not.present(extload))cycle
  ! compute body loads
  ! gravity load and add to extload
  ! WARNING: this must be changed to massdens_elmt
  if(gravity)extload(egdof)=extload(egdof)-eld*gam(mat_id(i_elmt))
  if(pseudoeq)then
    ! compute pseudostatic earthquake loads and add to extload
    eqload(1:nedof:3)=eqkx*eld(3:nedof:3)
    eqload(2:nedof:3)=eqky*eld(3:nedof:3)
    eqload(3:nedof:3)=eqkz*eld(3:nedof:3)
    extload(egdof)=extload(egdof)+eqload*gam(mat_id(i_elmt)) ! KN
  endif
enddo ! i_elmt=1,nelmt
!write(*,*)'complete!'
dprecon(0)=zero
if(present(extload))extload(0)=zero
end subroutine stiffness_bodyload
!===============================================================================

end module preprocess
!===============================================================================
