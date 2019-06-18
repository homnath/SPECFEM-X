module traction
implicit none
character(len=250),private :: myfname=' => apply_traction.f90'
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

! This routine reads and applies the traction specified in the traction file.
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
subroutine apply_traction(load,errcode,errtag)
use global
use math_constants
use element,only:hexface,hexface_sign,hex8_gnode
use integration,only:dshape_quad4_xy,dshape_quad4_yz,dshape_quad4_zx,          &
                     gll_weights_xy,gll_weights_yz,gll_weights_zx,             &
                     lagrange_gll_xy,lagrange_gll_yz,lagrange_gll_zx
use math_library,only : angle
use free_surface,only:nelmt_fs,nnode_fs,gnode_fs,gnum_fs,rgnum_fs,iface_fs
implicit none
real(kind=kreal),intent(inout) :: load(0:neq)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: lngdof(ndim)
integer :: i_face,i_gll,ios,iaxis
integer :: ielmt,iface,nface,inode,ifnode
integer :: num(nenode)
integer,allocatable :: numf(:)
integer :: nfdof,nfgll !face nodal dof, face gll pts
integer,allocatable :: fgdof(:) ! face global degrees of freedom
integer :: tractype,count_trac
logical :: trac_stat

real(kind=kreal) :: coord(ndim,4)
real(kind=kreal) :: x,x1,x2,detjac
real(kind=kreal),dimension(ndim) :: face_normal,dx_dxi,dx_deta,dq_dx,q,q1,q2
real(kind=kreal),allocatable :: ftracload(:) ! face traction load

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal

real(kind=kreal),dimension(:,:,:),allocatable :: dshape_quad4
real(kind=kreal),dimension(:),allocatable :: gll_weights
real(kind=kreal),dimension(:,:),allocatable :: lagrange_gll
real(kind=kreal),dimension(:,:,:),allocatable :: dlagrange_gll

real(kind=kreal) :: r0,r_cyl,y,z,z_cyl,theta_cyl,Tq,Jt,sigma_thetaz,G,theta_comp

character(len=80) :: fname
character(len=80) :: data_path

character(len=80) :: char80
integer :: dumi
! ufs must be real because it is read from the Ensight file
real,allocatable :: ufs(:,:),ufs0(:,:)

errtag="ERROR: unknown!"
errcode=-1
! set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

allocate(numf(maxngll2d))
allocate(fgdof(nndof*maxngll2d),ftracload(nndof*maxngll2d))
allocate(dshape_quad4(2,4,maxngll2d))
allocate(gll_weights(maxngll2d),lagrange_gll(maxngll2d,maxngll2d),             &
dlagrange_gll(2,maxngll2d,maxngll2d))

trac_stat=.true. ! Necessary for empty trfile
if(istraction)then
  fname=trim(data_path)//trim(trfile)//trim(ptail_inp)
  open(unit=11,file=trim(fname),status='old',action='read',iostat=ios)
  if (ios /= 0)then
    write(errtag,'(a)')'ERROR: input file "',trim(fname),'" cannot be opened!'
    return
  endif

  !read(11,*)ntrac
  count_trac=0
  traction: do
    read(11,*,iostat=ios)tractype
    if(ios/=0)exit traction
    count_trac=count_trac+1
    trac_stat=.false.

    if(tractype==0)then ! point loading
      read(11,*)q ! vector
      read(11,*)nface ! number of points
      do i_face=1,nface

        read(11,*)ielmt,inode
        lngdof=gdof(idofu,g_num(hex8_gnode(inode),ielmt))
        load(lngdof)=load(lngdof)+q
      enddo
      trac_stat=.true.
    elseif(tractype==1)then ! uniform loading
      read(11,*)q ! vector
      read(11,*)nface
      do i_face=1,nface
        read(11,*)ielmt,iface
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
          write(errtag,'(a)')'ERROR: wrong face ID for traction!'
          exit traction
        endif
        nfdof=nfgll*ndim

        num=g_num(:,ielmt)
        coord=g_coord(:,num(hexface(iface)%gnode))
        fgdof(1:nfdof)=reshape(gdof(idofu,g_num(hexface(iface)%node,ielmt)),(/nfdof/))

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
          ftracload(1:nfdof:3)=ftracload(1:nfdof:3)+ &
          q(1)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(1) !only in X direction
          ftracload(2:nfdof:3)=ftracload(2:nfdof:3)+ &
          q(2)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(2) !only in Y direction
          ftracload(3:nfdof:3)=ftracload(3:nfdof:3)+ &
          q(3)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(3) !only in Z direction
        enddo
        load(fgdof(1:nfdof))=load(fgdof(1:nfdof))+ftracload(1:nfdof)
      enddo
      trac_stat=.true.
    elseif(tractype==2)then ! linearly distributed loading
      read(11,*)iaxis,x1,x2,q1,q2 ! q1 and q2 are vectors, x1 and x2 can be any coordinates
      dq_dx=(q2-q1)/(x2-x1)
      read(11,*)nface
      do i_face=1,nface
        read(11,*)ielmt,iface
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
          nfgll=ngllzx
          lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
          gll_weights(1:nfgll)=gll_weights_xy
          dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
        else
          write(errtag,'(a)')'ERROR: wrong face ID for traction!'
          exit traction
        endif
        nfdof=nfgll*ndim

        num=g_num(:,ielmt)
        coord=g_coord(:,num(hexface(iface)%gnode))
        fgdof(1:nfdof)=reshape(gdof(idofu,g_num(hexface(iface)%node,ielmt)),(/nfdof/))
        ftracload=zero
        ! compute numerical integration
        do i_gll=1,nfgll
          x=g_coord(iaxis,num(hexface(iface)%node(i_gll)))
          q=q1+dq_dx*(x-x1) ! vector of nodal values

          ! compute two vectors dx_dxi and dx_deta
          dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
          dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))

          ! Normal = (dx_dxi x dx_deta)
          face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
          face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
          face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

          detjac=sqrt(dot_product(face_normal,face_normal))
          face_normal=hexface_sign(iface)*face_normal/detjac

          ftracload(1:nfdof:3)=ftracload(1:nfdof:3)+ &
          q(1)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(1) !only in X direction
          ftracload(2:nfdof:3)=ftracload(2:nfdof:3)+ &
          q(2)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(2) !only in Y direction
          ftracload(3:nfdof:3)=ftracload(3:nfdof:3)+ &
          q(3)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(3) !only in Z direction
       enddo
       load(fgdof(1:nfdof))=load(fgdof(1:nfdof))+ftracload(1:nfdof)
     enddo
     trac_stat=.true.
    elseif(tractype==23)then ! torsion
      ! WARNING: this is implemented ONLY for simple particular torsion of a cylinder
      ! compute polar moment of inertia Jt

      read(11,*)r0,Tq ! Torque
      Jt=half*pi*r0*r0*r0*r0
      G=half*ym_blk(1)/(1+nu_blk(1))
      theta_comp=Tq*0.2/(G*Jt)
      print*,'Comuted twisting angle:',half,ym_blk(1),nu_blk(1),theta_comp

      read(11,*)nface
      do i_face=1,nface
        read(11,*)ielmt,iface
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
          nfgll=ngllzx
          lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
          gll_weights(1:nfgll)=gll_weights_xy
          dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
        else
          write(errtag,'(a)')'ERROR: wrong face ID for traction!'
          exit traction
        endif
        nfdof=nfgll*ndim

        num=g_num(:,ielmt)
        coord=g_coord(:,num(hexface(iface)%gnode))
        fgdof(1:nfdof)=reshape(gdof(idofu,g_num(hexface(iface)%node,ielmt)),(/nfdof/))
        ftracload=zero
        ! compute numerical integration
        do i_gll=1,nfgll
          ! x,y,z coordinates
          x=g_coord(1,num(hexface(iface)%node(i_gll)))
          y=g_coord(2,num(hexface(iface)%node(i_gll)))
          z=g_coord(3,num(hexface(iface)%node(i_gll)))

          ! compute r, \theta, z
          r_cyl=sqrt(x*x+y*y)

          theta_cyl=angle(x,y)
          z_cyl=z

          ! compute \sigma_{\theta z}
          sigma_thetaz=Tq*r_cyl/Jt
          q=zero ! vector of nodal values

          ! compute two vectors dx_dxi and dx_deta
          dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
          dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))

          ! Normal = (dx_dxi x dx_deta)
          face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
          face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
          face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

          detjac=sqrt(dot_product(face_normal,face_normal))
          face_normal=hexface_sign(iface)*face_normal/detjac

          ! compute \sigma_{ij}.n_j
          q(1)=-sin(theta_cyl)*sigma_thetaz*face_normal(3)
          q(2)=cos(theta_cyl)*sigma_thetaz*face_normal(3)

          ftracload(1:nfdof:3)=ftracload(1:nfdof:3)+ &
          q(1)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(1) !only in X direction
          ftracload(2:nfdof:3)=ftracload(2:nfdof:3)+ &
          q(2)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(2) !only in Y direction
          ftracload(3:nfdof:3)=ftracload(3:nfdof:3)+ &
          q(3)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) ! *face_normal(3) !only in Z direction
       enddo
       load(fgdof(1:nfdof))=load(fgdof(1:nfdof))+ftracload(1:nfdof)
     enddo
     trac_stat=.true.
    else
      write(errtag,'(a)')'ERROR: traction type ',tractype,' not supported!'
      exit traction
    endif
  enddo traction
  close(11)
endif !(istraction)

if(.not.trac_stat)then
  write(errtag,'(a)')'ERROR: all tractions cannot be read!'
  return
endif

! Traction defined on the free surface on SEM points.
! WARNING: valid only for the FREE SURFACE already determined in the program.
if(isfstraction.and.nnode_fs>0)then
  fname=trim(trfspath)//trim(trfsfile)//trim(ptail_inp)//'.dis'
  open(unit=11,file=trim(fname),status='old',action='read',access='stream', &
  form='unformatted',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  
  allocate(ufs(nnode_fs,nndofu))
  read(11)char80
  read(11)char80
  read(11)dumi
  read(11)char80
  read(11)ufs
  close(11)
  if(isfstraction0)then
    fname=trim(trfspath0)//trim(trfsfile0)//trim(ptail_inp)//'.dis'
    open(unit=11,file=trim(fname),status='old',action='read',access='stream', &
    form='unformatted',iostat = ios)
    if( ios /= 0 ) then
      write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
      return
    endif
    
    allocate(ufs0(nnode_fs,nndofu))
    read(11)char80
    read(11)char80
    read(11)dumi
    read(11)char80
    read(11)ufs0
    close(11)
    ufs=ufs-ufs0
  endif

  do i_face=1,nelmt_fs
    iface=iface_fs(i_face)
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
      nfgll=ngllzx
      lagrange_gll(1:nfgll,1:nfgll)=lagrange_gll_xy
      gll_weights(1:nfgll)=gll_weights_xy
      dshape_quad4(:,:,1:nfgll)=dshape_quad4_xy
    else
      write(errtag,'(a)')'ERROR: wrong face ID for traction!'
      return
    endif
    nfdof=nfgll*ndim

    numf=gnum_fs(:,i_face)
    coord=g_coord(:,numf)
    fgdof(1:nfdof)=reshape(gdof(idofu,numf),(/nfdof/))
    ftracload=zero
    ! compute numerical integration
    do i_gll=1,nfgll
      ifnode=rgnum_fs(i_gll,i_face)
      q=ufs(ifnode,:) ! vector of nodal values

      ! compute two vectors dx_dxi and dx_deta
      dx_dxi=matmul(coord,dshape_quad4(1,:,i_gll))
      dx_deta=matmul(coord,dshape_quad4(2,:,i_gll))

      ! Normal = (dx_dxi x dx_deta)
      face_normal(1)=dx_dxi(2)*dx_deta(3)-dx_deta(2)*dx_dxi(3)
      face_normal(2)=dx_deta(1)*dx_dxi(3)-dx_dxi(1)*dx_deta(3)
      face_normal(3)=dx_dxi(1)*dx_deta(2)-dx_deta(1)*dx_dxi(2)

      detjac=sqrt(dot_product(face_normal,face_normal))
      face_normal=hexface_sign(iface)*face_normal/detjac

      ftracload(1:nfdof:3)=ftracload(1:nfdof:3)+ &
      q(1)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) !*face_normal(1) !only in X direction
      ftracload(2:nfdof:3)=ftracload(2:nfdof:3)+ &
      q(2)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) !*face_normal(2) !only in Y direction
      ftracload(3:nfdof:3)=ftracload(3:nfdof:3)+ &
      q(3)*lagrange_gll(i_gll,:)*detjac*gll_weights(i_gll) !*face_normal(3) !only in Z direction
      
      !ftracload(1:nfdof:3)=q(1)
      !ftracload(2:nfdof:3)=q(2)
      !ftracload(3:nfdof:3)=q(3)
    enddo
    load(fgdof(1:nfdof))=load(fgdof(1:nfdof))+ftracload(1:nfdof)
  enddo
  
  deallocate(ufs)

endif !(isfstraction)

deallocate(numf)
deallocate(fgdof,ftracload)
deallocate(dshape_quad4)
deallocate(gll_weights,lagrange_gll,dlagrange_gll)

errcode=0

return

end subroutine apply_traction
!===============================================================================

end module traction
!===============================================================================
