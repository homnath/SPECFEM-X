module mtraction
contains
!-------------------------------------------------------------------------------

! This subroutine routine read and applies the magnetic traction specified in the 
! mtraction file
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
subroutine apply_mtraction(load,errcode,errtag)
use global
use math_constants
use element,only:hexface,hexface_sign
use integration,only:dshape_quad4_xy,dshape_quad4_yz,dshape_quad4_zx,          &
                     gll_weights_xy,gll_weights_yz,gll_weights_zx,             &
                     lagrange_gll_xy,lagrange_gll_yz,lagrange_gll_zx
use math_library,only : angle,magnetic_unitvec
use conversion_constants,only:DEG2RAD
implicit none
real(kind=kreal),intent(inout) :: load(0:neq)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: lndofphi,lngdofu(ndim)
integer :: i_face,i_gll,ios,iaxis
integer :: ielmt,iface,nface,inode
integer :: num(nenode)
integer :: nfdofphi,nfgll !face nodal dof, face gll pts
integer,allocatable :: fgdof(:) ! face global degrees of freedom
integer :: tractype,count_trac
logical :: trac_stat

real(kind=kreal) :: coord(ndim,4)
real(kind=kreal) :: x,x1,x2,detjac
real(kind=kreal),dimension(ndim) :: face_normal,dx_dxi,dx_deta
real(kind=kreal),allocatable :: ftracload(:) ! face mtraction load

! Magnetization
! magnetization properties type for the traction surface
! 'inherit': inherit the magnetization form the parent element
! 'define': defien the magnetization
character(len=10) :: mag_type
! Inclination (0) or latitude (1). If latitude, the inclination is obtained
! using the relation: inclination = ATAN(2*TAN(latitude))
integer :: incORlat
integer :: imat,imatmag
! magnitude of the magnetization
real(kind=kreal) :: M0,M(NDIM)
real(kind=kreal) :: inc,dec,azim,rval
real(kind=kreal),allocatable :: Mgll(:,:)

real(kind=kreal),dimension(:,:,:),allocatable :: dshape_quad4
real(kind=kreal),dimension(:),allocatable :: gll_weights
real(kind=kreal),dimension(:,:),allocatable :: lagrange_gll
real(kind=kreal),dimension(:,:,:),allocatable :: dlagrange_gll

character(len=80) :: fname
character(len=80) :: data_path

errtag="ERROR: unknown!"
errcode=-1
! set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

fname=trim(data_path)//trim(mtrfile)//trim(ptail_inp)
open(unit=11,file=trim(fname),status='old',action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "',trim(fname),'" cannot be opened!'
  return
endif

allocate(fgdof(nndof*maxngll2d),ftracload(nndof*maxngll2d))
allocate(dshape_quad4(2,4,maxngll2d))
allocate(gll_weights(maxngll2d),lagrange_gll(maxngll2d,maxngll2d),             &
dlagrange_gll(2,maxngll2d,maxngll2d))

!read(11,*)ntrac
trac_stat=.true. ! necessary for empty trfile
count_trac=0
mtraction: do
  read(11,*,iostat=ios)tractype
  if(ios/=0)exit mtraction
  count_trac=count_trac+1
  trac_stat=.false.

  if(tractype==0)then
    ! inherit magnetization properties from the parent element
    mag_type='inherit'
  elseif(tractype==1)then
    ! define the uniform magnetization properties
    mag_type='uniform'
    read(11,*)M0,incORlat,rval,dec,azim
    if(incORlat==0)then
      inc=rval*DEG2RAD
    elseif(incORlat==1)then
      inc=atan(TWO*tan(rval*DEG2RAD))
    else
      write(errtag,'(a,i0)')'ERROR: incORlat must be either 0 or 1 but has a &
      &value: ',incORlat
      return
    endif
    dec=dec*DEG2RAD
    azim=azim*DEG2RAD
    M=M0*magnetic_unitvec(inc,dec,azim)
  else
    write(errtag,'(a)')'ERROR: mtraction type ',tractype,' not supported!'
    return
  endif
    
  read(11,*)nface
  do i_face=1,nface
    read(11,*)ielmt,iface
    imat=mat_id(ielmt)
    if(.not.ismagnet_blk(imat))then
      write(errtag,'(a)')'ERROR: mtraction element is unmagnetic!'
      return
    endif
    !if(trim(mag_type).eq.'inherit')then
    !  imatmag=imat_to_imatmag(imat)
    !  M=magnetization_blk(:,imatmag)
    !endif
    if(iface==1 .or. iface==3)then
      nfgll=ngllzx
      lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_zx
      gll_weights(1:nfgll)=gll_weights_zx
      dshape_quad4(:,:,1:nfgll)=dshape_quad4_zx
    elseif(iface==2 .or. iface==4)then
      nfgll=ngllyz
      lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_yz
      gll_weights(1:nfgll)=gll_weights_yz
      dshape_quad4(:,:,1:nfgll)=dshape_quad4_yz
    elseif(iface==5 .or. iface==6)then
      nfgll=ngllxy
      lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
      gll_weights(1:nfgll)=gll_weights_xy
      dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
    else
      write(errtag,'(a)')'ERROR: wrong face ID for mtraction!'
      return
    endif
    allocate(Mgll(NDIM,nfgll))
    Mgll=ZERO
    if(trim(mag_type).eq.'uniform')then
      do i_gll=1,nfgll
        Mgll(:,i_gll)=M
      enddo
    elseif(trim(mag_type).eq.'inherit')then
      do i_gll=1,nfgll
        Mgll(:,i_gll)=magnetization_elmt(:,hexface(iface)%node(i_gll),ielmt)
      enddo
    else
      write(errtag,'(a)')'ERROR: unsupported mag_type ',trim(mag_type),'!'
      return
    endif
    nfdofphi=nfgll*nndofphi

    num=g_num(:,ielmt)
    coord=g_coord(:,num(hexface(iface)%gnode))
    fgdof(1:nfdofphi)=reshape(gdof(idofphi,g_num(hexface(iface)%node,ielmt)),(/nfdofphi/))
    ftracload=zero
    ! compute numerical integration
    do i_gll=1,nfgll
      ! compute d(area)
      dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
      dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))
      ! Normal
      face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
      face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
      face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

      detjac=sqrt(dot_product(face_normal,face_normal))
      face_normal=hexface_sign(iface)*face_normal/detjac
      ! TODO:for constant q this can be computed only once!!
      M=Mgll(:,i_gll)
      ftracload(1:nfdofphi)=ftracload(1:nfdofphi)- &
      dot_product(M,face_normal)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll)
    enddo ! i_gll
    load(fgdof(1:nfdofphi))=load(fgdof(1:nfdofphi))+ftracload(1:nfdofphi)
    deallocate(Mgll)
  enddo ! i_face
  trac_stat=.true.
enddo mtraction
close(11)
deallocate(dshape_quad4)
deallocate(gll_weights,lagrange_gll,dlagrange_gll)
if(.not.trac_stat)then
  write(errtag,'(a)')'ERROR: all tractions cannot be read!'
  return
endif

errcode=0

return

end subroutine apply_mtraction
!===============================================================================

end module mtraction
!===============================================================================
