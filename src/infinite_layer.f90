! DESCRIPTION
!  This module contains the infinite-element routines which assist to create
!  infinite element layer.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Apr 11,2012; HNG, Jul 12,2011; HNG, Apr 09,2010
! TODO
!  - 
module infinite_layer
use set_precision
use global,only : ndim,pole_coord0,pole_coord1
implicit none
character(len=250),private :: myfname=' => infinite_layer'
character(len=500),private :: errsrc

integer,parameter,private :: xdir=1,ydir=2,zdir=3
integer,allocatable,private :: elmt_face_pair(:,:)
integer,private :: nelmt_inf,nface_inf
contains
!-------------------------------------------------------------------------------

! This subroutine add a layer of infinite mesh outside the model given the
! reference surface and infinite element information.  
! REVISION
!   HNG, Jul 12,2011; ; HNG, Apr 09,2010
subroutine add_infinite_elements(errcode,errtag)
use global
use math_constants,only:zero
use math_library,only:distance,i_uniinv
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: i,ios
integer :: i_elmt,ielmt,iface,nelpart,i_elpart
! local node numbers in each face
integer,dimension(6,4) :: node_face

real(kind=kreal) :: gaminf,x0(ndim)
integer :: g_numOLD(8,nelmt),mat_idOLD(nelmt),nelmtOLD,nnodeOLD
real(kind=kreal) :: r1,g_coordOLD(ndim,nnode)
real(kind=kreal),allocatable :: xs(:,:),mirxs(:,:)

integer :: n1,n2,nelmtTRINF,nelmtINF,nsnode,nsnode_all
integer,allocatable :: nodelist(:),inode_order(:),g_numinf(:),iface_elmt(:)
! material DOMAINs of infinite surface elements
integer,allocatable :: mat_domainINFS(:)
! material IDs of infinite surface elements
integer :: i_mat,imat,imat_vis,istat,mdomain
integer :: nmatblk_inherit,nmatblkOLD,nmatblk_viscoelasOLD,nvis_inherit
logical,allocatable :: ismat(:),ismat_vis(:)
real(kind=kreal) :: rho_or_gam
integer,allocatable :: imat_inherit(:),imat_INFS(:)
character(len=250) :: inherit_matfile

integer,allocatable :: mat_domainOLD(:),type_blkOLD(:)
real(kind=kreal),allocatable :: gam_blkOLD(:),rho_blkOLD(:),ym_blkOLD(:), &
nu_blkOLD(:),coh_blkOLD(:),phi_blkOLD(:),psi_blkOLD(:)
character(len=60),allocatable :: mfile_blkOLD(:)

real(kind=kreal),allocatable :: muratio_blkOLD(:,:),viscosity_blkOLD(:,:)
integer,allocatable :: imat_to_imatveOLD(:),imatve_to_imatOLD(:)

logical,allocatable :: isnode(:)
logical :: reset_pole
character(len=250) :: fname
character(len=150) :: data_path

errtag="ERROR: unknown!"
errcode=-1
! set data path
if(ismpi .and. nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

! local node numbering in each face CUBIT/EXODUS convention
node_face(1,:)=(/1,2,6,5/) ! counterclockwise w.r.t. outer normal
node_face(2,:)=(/2,3,7,6/) ! counterclockwise w.r.t. outer normal
node_face(3,:)=(/4,3,7,8/) ! clockwise w.r.t. outer normal
node_face(4,:)=(/1,4,8,5/) ! clockwise w.r.t. outer normal
node_face(5,:)=(/1,2,3,4/) ! clockwise w.r.t. outer normal
node_face(6,:)=(/5,6,7,8/) ! counterclockwise w.r.t. outer normal

fname=trim(data_path)//trim(infrfile)//trim(ptail_inp)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if(ios.ne.0) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif

read(11,*,iostat=ios)nelpart
if(ios.ne.0)then
  write(errtag,*)'ERROR: bad infrfile!'
  return
endif
nelmtINF=nelpart
write(logunit,*)'Total number of infinite elements:',nelmtINF
flush(logunit)
allocate(iface_elmt(nelmtINF))
allocate(mat_domainINFS(nelmtINF),imat_INFS(nelmtINF))
nsnode_all=4*nelmtINF
allocate(nodelist(nsnode_all),inode_order(nsnode_all))
allocate(ismat(nmatblk),ismat_vis(nmatblk))
ismat=.false.
ismat_vis=.false.
n1=1; n2=4
do i_elpart=1,nelpart
  ! This will read a line and proceed to next line
  read(11,*)ielmt,iface
  imat=mat_id(ielmt)
  iface_elmt(i_elpart)=iface
  nodelist(n1:n2)=g_num(node_face(iface,:),ielmt)
  n1=n2+1; n2=n1+3
  mat_domainINFS(i_elpart)=mat_domain(imat)
  imat_INFS(i_elpart)=imat
  ismat(imat)=.true.
  if(mat_domain(imat).eq.VISCOELASTIC_DOMAIN .or. &
     mat_domain(imat).eq.VISCOELASTIC_TRINFDOMAIN .or. &
     mat_domain(imat).eq.VISCOELASTIC_INFDOMAIN)then
    ismat_vis(imat)=.true.
  endif
enddo
close(11)

nmatblk_inherit=0
! if the material block IDs have to be inherited, prepare inheritance
! information and create inherited material list file 
if(trim(matinf_type).eq.'inherit')then
  nmatblk_inherit=count(ismat)
  allocate(imat_inherit(nmatblk))

  ! reassign material blocks
  allocate(mat_domainOLD(nmatblk),type_blkOLD(nmatblk),gam_blkOLD(nmatblk),    &
  rho_blkOLD(nmatblk),ym_blkOLD(nmatblk),coh_blkOLD(nmatblk), &
  nu_blkOLD(nmatblk),phi_blkOLD(nmatblk),psi_blkOLD(nmatblk))
  
  mat_domainOLD=mat_domain
  type_blkOLD=type_blk
  mfile_blkOLD=mfile_blk
  gam_blkOLD=gam_blk
  rho_blkOLD=rho_blk
  ym_blkOLD=ym_blk
  nu_blkOLD=nu_blk
  coh_blkOLD=coh_blk
  phi_blkOLD=phi_blk
  psi_blkOLD=psi_blk

  nmatblkOLD=nmatblk
  
  deallocate(mat_domain,type_blk,gam_blk,             &
  rho_blk,ym_blk,coh_blk,nu_blk,           &
  phi_blk,psi_blk)
  deallocate(mfile_blk)

  nmatblk=nmatblkOLD+nmatblk_inherit
  
  allocate(mat_domain(nmatblk),type_blk(nmatblk),gam_blk(nmatblk),             &
  rho_blk(nmatblk),ym_blk(nmatblk),coh_blk(nmatblk),nu_blk(nmatblk),           &
  phi_blk(nmatblk),psi_blk(nmatblk))
  allocate(mfile_blk(nmatblk))
  
  mat_domain(1:nmatblkOLD)=mat_domainOLD
  type_blk(1:nmatblkOLD)=type_blkOLD
  mfile_blk(1:nmatblkOLD)=mfile_blkOLD
  gam_blk(1:nmatblkOLD)=gam_blkOLD
  rho_blk(1:nmatblkOLD)=rho_blkOLD
  ym_blk(1:nmatblkOLD)=ym_blkOLD
  nu_blk(1:nmatblkOLD)=nu_blkOLD
  coh_blk(1:nmatblkOLD)=coh_blkOLD
  phi_blk(1:nmatblkOLD)=phi_blkOLD
  psi_blk(1:nmatblkOLD)=psi_blkOLD

  nmatblk_viscoelasOLD=nmatblk_viscoelas
  nvis_inherit=count(ismat_vis)
  if(nvis_inherit.gt.0)then
    allocate(muratio_blkOLD(nmaxwell,nmatblk_viscoelasOLD), &
    viscosity_blkOLD(nmaxwell,nmatblk_viscoelasOLD))
    allocate(imat_to_imatveOLD(nmatblkOLD),imatve_to_imatOLD(nmatblk_viscoelasOLD))
    muratio_blkOLD=muratio_blk
    viscosity_blkOLD=viscosity_blk
    imat_to_imatveOLD=imat_to_imatve
    imatve_to_imatOLD=imatve_to_imat
    
    deallocate(muratio_blk,viscosity_blk)
    deallocate(imat_to_imatve,imatve_to_imat)

    nmatblk_viscoelas=nmatblk_viscoelasOLD+nvis_inherit

    allocate(muratio_blk(nmaxwell,nmatblk_viscoelas), &
    viscosity_blk(nmaxwell,nmatblk_viscoelas))
    allocate(imat_to_imatve(nmatblk),imatve_to_imat(nmatblk_viscoelas))
    
    imat_to_imatve(1:nmatblkOLD)=imat_to_imatveOLD
    imatve_to_imat(1:nmatblk_viscoelasOLD)=imatve_to_imatOLD
    muratio_blk(:,1:nmatblk_viscoelasOLD)=muratio_blkOLD
    viscosity_blk(:,1:nmatblk_viscoelasOLD)=viscosity_blkOLD
  endif

  imat_inherit=-9999
  imat=nmatblkOLD
  imat_vis=nmatblk_viscoelasOLD
  do i_mat=1,nmatblkOLD
    if(ismat(i_mat))then
      imat=imat+1
      imat_inherit(i_mat)=imat
      mdomain=mat_domainOLD(i_mat)
      ! make mat_domain infinite
      if(mdomain.lt.ELASTIC_TRINFDOMAIN)then
        mdomain=mdomain*1000
      elseif(mdomain.eq.ELASTIC_TRINFDOMAIN)then
        mdomain=mdomain*10
      else
        write(stdout,*)'ERROR: wrong material ID for inhertance!'
        stop
      endif
      ! inherit properties
      mat_domain(imat)=mdomain
      type_blk(imat)=type_blkOLD(i_mat)
      mfile_blk(imat)=mfile_blkOLD(i_mat)
      gam_blk(imat)=gam_blkOLD(i_mat)
      rho_blk(imat)=rho_blkOLD(i_mat)
      ym_blk(imat)=ym_blkOLD(i_mat)
      nu_blk(imat)=nu_blkOLD(i_mat)
      coh_blk(imat)=coh_blkOLD(i_mat)
      phi_blk(imat)=phi_blkOLD(i_mat)
      psi_blk(imat)=psi_blkOLD(i_mat)
    endif

    ! viscelastic blocks
    if(ismat_vis(i_mat))then
      imat_vis=imat_vis+1
      imat_to_imatve(imat)=imat_vis
      imatve_to_imat(imat_vis)=imat
      muratio_blk(:,imat_vis)=muratio_blk(:,i_mat)
      viscosity_blk(:,imat_vis)=viscosity_blk(:,i_mat)
    endif
  enddo
  ! Write inherited material list file.
  ! Save infinite element file.
  inherit_matfile=trim(inp_path)//trim(matfile)//'_inherit'
  open(unit=17,file=trim(inherit_matfile),status='replace',action='write', &
  iostat=istat)
  if(istat/=0) then
    print*,'ERROR: cannot open inherited material file:',trim(inherit_matfile)
    stop
  endif
  write(17,'(a)')'# material properties (id,domain,type,gamma,ym,nu,phi,coh,psi)'
  write(17,'(i0)')nmatblk
  do i_mat=1,nmatblk
    if(isdensity)then
      rho_or_gam=rho_blk(i_mat)
    else
      rho_or_gam=gam_blk(i_mat)
    endif
    write(17,'(i0,a,i0,a,i0,a,g0.6,a,es14.6,5(a,g0.6))')i_mat,', ', &
    mat_domain(i_mat),', ',type_blk(i_mat),', ',rho_or_gam,', ', &
    ym_blk(i_mat),', ',nu_blk(i_mat),', ',phi_blk(i_mat),', ', &
    coh_blk(i_mat),', ',psi_blk(i_mat)
  enddo
  ! write viscoelastic properties
  ! TODO: WARNING: this segment must be modified for nmaxwell>1 and  
  ! rheologies other than general maxwell family
  if(nmatblk_viscoelas.gt.0)then
    write(17,'(a)')trim(visco_rheology)
    do i_mat=1,nmatblk_viscoelas
      imat=imatve_to_imat(i_mat)
      write(17,'(i0,a,g0.6,a,es14.6)')imat,', ', &
      muratio_blk(:,i_mat),', ',viscosity_blk(:,i_mat)
    enddo
  endif
  ! write water properties 
  if(iswater)then
    write(17,'(i0)')nwmat
    do i=1,nwmat
      write(17,'(i0)')waterid(i)
    enddo
  endif
  ! magnetization
  if(nmatblk_magnet.gt.0)then
    write(17,'(a)')'magnetization:'
    do i_mat=1,nmatblk_magnet
      imat=imatmag_to_imat(i_mat)
      write(17,'(i0,a,g0.6,a,es14.6)')imat,', ', &
      muratio_blk(:,i_mat),', ',viscosity_blk(:,i_mat)
    enddo
  endif
  close(17)
  deallocate(ismat)
  deallocate(mat_domainOLD,type_blkOLD,gam_blkOLD,             &
  rho_blkOLD,ym_blkOLD,coh_blkOLD,nu_blkOLD,           &
  phi_blkOLD,psi_blkOLD)
  deallocate(mfile_blkOLD)
  if(nvis_inherit.gt.0)then
    deallocate(muratio_blkOLD,viscosity_blkOLD)
    deallocate(imat_to_imatveOLD,imatve_to_imatOLD)
  endif
endif

call i_uniinv(nodelist,inode_order)

nsnode=maxval(inode_order)
allocate(isnode(nsnode),xs(ndim,nsnode),mirxs(ndim,nsnode),g_numinf(nsnode))

isnode=.false.
! assign xs
xs(:,inode_order(1))=g_coord(:,nodelist(1))
isnode(inode_order(1))=.true.
do i=2,nsnode_all
  if(.not.isnode(inode_order(i)))then
     xs(:,inode_order(i))=g_coord(:,nodelist(i))
     isnode(inode_order(i))=.true.
  endif
enddo
deallocate(isnode)

! compute mirror nodes
do i=1,nsnode
  ! pole specifies the reference point for the decaying functions
  x0=pole_coord0
  if(trim(pole_type)=='pointaxis')then
    if(xs(pole_axis,i).ge.axis_range(2))x0=pole_coord1
  endif

  ! set pole axis if option provided
  reset_pole=.false.
  ! reset all poles to axis
  if(trim(pole_type)=='axis')then
    reset_pole=.true.
  endif
  ! reset some poles to axis
  if(trim(pole_type)=='pointaxis')then
    if(xs(pole_axis,i).gt.axis_range(1).and.xs(pole_axis,i).lt.axis_range(2))reset_pole=.true.
  endif
  if(reset_pole)then
    x0(pole_axis)=xs(pole_axis,i)
  endif
  
  r1=distance(x0,xs(:,i),ndim)
  if(rinf.le.r1)then
    write(errtag,*)'ERROR: reference infinite radius is smaller than the model!'
    return
  endif
  gaminf=r1/(rinf-r1)
  ! division formula
  mirxs(:,i)=((gaminf+one)*xs(:,i)-x0)/gaminf
  g_numinf(i)=nnode+i
enddo
!stop
deallocate(xs)
g_numOLD=g_num;
g_coordOLD=g_coord;
mat_idOLD=mat_id
deallocate(g_num,g_coord,mat_id)

nelmtOLD=nelmt; nnodeOLD=nnode

nelmt=nelmtOLD+nelmtINF
nnode=nnodeOLD+nsnode

! reallocate global node - and element-arrays
allocate(g_num(8,nelmt),g_coord(ndim,nnode),mat_id(nelmt))

! update connectivity, coordinates, and material IDs
g_num(:,1:nelmtOLD)=g_numOLD
g_coord(:,1:nnodeOLD)=g_coordOLD
mat_id(1:nelmtOLD)=mat_idOLD
!! add to material list
!mat_id(nelmtOLD+1:nelmt)=imat_inf
! add to global node list
g_coord(:,nnodeOLD+1:nnode)=mirxs
deallocate(mirxs)

! add to global element list
! to have face 5 of new element always clockwise w.r.t. outer normal 
! face 3 or 4 or 5 of reference element has to be reordered
nelmtTRINF=0
n1=1; n2=4
do i=1,nelmtINF
  ielmt=nelmtOLD+i
  if(iface_elmt(i).eq.3.or.iface_elmt(i).eq.4.or.iface_elmt(i).eq.5)then
    g_num(1:4,ielmt)=nodelist(n2:n1:-1)
    g_num(5:8,ielmt)=g_numinf(inode_order(n2:n1:-1))
  else
    g_num(1:4,ielmt)=nodelist(n1:n2)
    g_num(5:8,ielmt)=g_numinf(inode_order(n1:n2))
  endif
  ! assign material ID
  if(trim(matinf_type).eq.'inherit')then
  ! inherit material block ID/s from the parent elements
    imat=imat_INFS(i)
    mat_id(ielmt)=imat_inherit(imat)
  else
  ! set specified material block ID/s
    if(mat_domainINFS(i).eq.ELASTIC_TRINFDOMAIN)then
      ! adjacent element is the transition element
      if(imat_trinf.lt.1)then
        write(*,'(a,i0)')'ERROR: invalid material ID for transition infinite elements! ',imat_trinf
        write(*,'(a)')'HINT: Make sure that the material_list file is properly defined!'
        write(*,'(a)')'      If there is a transition layer, make sure that the'
        write(*,'(a)')'      paramater "imat_trinf" is defined for "bc:" line in the main input file.'
        stop
      endif
      mat_id(ielmt)=imat_trinf
      nelmtTRINF=nelmtTRINF+1
    else
      mat_id(ielmt)=imat_inf
    endif
  endif
  n1=n2+1; n2=n1+3
enddo
deallocate(nodelist,inode_order,g_numinf,iface_elmt, &
mat_domainINFS,imat_INFS)
if(allocated(imat_inherit))deallocate(imat_inherit)
write(logunit,*)'Number of transition infinite elements:',nelmtTRINF
flush(logunit)

ielmtINF1=nelmtOLD+1; ielmtINF2=nelmt
ifaceINF=6

! if infinite elements are internally added, each infinite elements always has 
! exactly one infinite face
allocate(elmt_face_pair(3,nelmtINF)) !elementID, faceID, infDIR
elmt_face_pair=-9999

! save infinite elements
! each elements has one infinite face
nface_inf=nelmtINF
elmt_face_pair(1,:)=(/ (i_elmt,i_elmt=ielmtINF1,ielmtINF2) /)
elmt_face_pair(2,:)=ifaceINF
! infinite face is always 6 or the top face, i.e., z direction
elmt_face_pair(3,:)=zdir

call save_infinite_elements()

! compute nodal to global
errcode=0
return

end subroutine add_infinite_elements
!===============================================================================

!-------------------------------------------------------------------------------
! This subroutine reorder infinite elements such that the infinite faces are
! always located on the outer directions. This subroutine basically changes the 
! connectivity of the infinite elements if necessary.
! Note: Face ordering in EXODUS/CUBIT convention:
!--------------------------
! Face      Face ID
!--------------------------
! xmin         4
! xmax         2
! ymin         1
! ymax         3
! zmin         5
! zmax         6
!--------------------------
subroutine reorder_infinite_elements(errcode,errtag)
use global
use math_constants,only:zero
use math_library,only:sort_array_coord
use element,only:map2exodus_hex8,get_faceidHEX8
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
integer :: i,i_elmt,imat,inum,istat
integer :: ielmt
integer,parameter :: ignodes(8)=(/ (i,i=1,8) /)
integer :: ncount,nelmt_ing0,neface,xnelmt_inf
integer :: gnum(ngnode)
integer,allocatable :: ielmt_inf(:),g_numinf(:,:)

integer :: nxmin,nxmax,nymin,nymax,nzmin,nzmax
real(kind=kreal) :: xmin,xmax,ymin,ymax,zmin,zmax
real(kind=kreal) :: xp0(ngnode),yp0(ngnode),zp0(ngnode)
real(kind=kreal) :: xp(ngnode),yp(ngnode),zp(ngnode)

integer :: iface
integer :: i_node,inode,inodes(4)
logical :: isinf
logical :: isxmin(8),isxmax(8),isymin(8),isymax(8),iszmin(8),iszmax(8)
character(len=250) :: inf_file

! Format string
integer :: ewidth
character(len=20) :: format_ef

errtag="ERROR: unknown!"
errcode=-1

!nelmt_ing0=count(mat_id.eq.imat_inf)
nelmt_ing0=count(mat_domain(mat_id).ge.ELASTIC_INFDOMAIN)

write(*,'(a,i0)')'number of original infinite elements:',nelmt_ing0
allocate(ielmt_inf(nelmt_ing0),g_numinf(ngnode,nelmt_ing0))
allocate(elmt_face_pair(3,nelmt_ing0*3)) !elementID, faceID, infDIR
elmt_face_pair=-9999
! Identify infinite elements
inum=0
do i_elmt=1,nelmt
  !if(mat_id(i_elmt).eq.imat_inf)then
  imat=mat_id(i_elmt)
  if(mat_domain(imat).ge.ELASTIC_INFDOMAIN)then
    inum=inum+1
    ielmt_inf(inum)=i_elmt
    g_numinf(:,inum)=g_num(:,i_elmt)
  endif
enddo
! Check the number of infinite elements
if(inum.ne.nelmt_ing0)then
  write(errtag,*)'ERROR: inconsistent number of infinite elements!'
  return
endif

! Compute domain range
! This is in SERIAL mode
xmin=minval(g_coord(1,:))
xmax=maxval(g_coord(1,:))
ymin=minval(g_coord(2,:))
ymax=maxval(g_coord(2,:))
zmin=minval(g_coord(3,:))
zmax=maxval(g_coord(3,:))

!print*,xmin,xmax,ymin,ymax,zmin,zmax
nface_inf=0
nelmt_inf=0
do i_elmt=1,nelmt_ing0
  ielmt=ielmt_inf(i_elmt)
  gnum=g_numinf(:,i_elmt)
  xp=g_coord(1,gnum)
  yp=g_coord(2,gnum)
  zp=g_coord(3,gnum)

  xp0=xp
  yp0=yp
  zp0=zp

  isxmin=xp.eq.xmin
  isxmax=xp.eq.xmax
  isymin=yp.eq.ymin
  isymax=yp.eq.ymax
  iszmin=zp.eq.zmin
  iszmax=zp.eq.zmax
  
  nxmin=count(isxmin)
  nxmax=count(isxmax)
  nymin=count(isymin)
  nymax=count(isymax)
  nzmin=count(iszmin)
  nzmax=count(iszmax)

  neface=0
  isinf=.false.
  ! XMIN face
  if(nxmin.eq.4)then
    neface=neface+1

    ! extract active nodes
    inode=0
    do i_node=1,8
      if(isxmin(i_node))then
        inode=inode+1
        inodes(inode)=i_node
      endif
    enddo
    if(inode.ne.4)then
      write(errtag,*)'ERROR: number of face nodes must be 4!',inode
      return
    endif
   
    ! get face ID
    !iface=4 !get_faceidHEX8(inodes,ignodes)
    iface=get_faceidHEX8(inodes,ignodes)

    !! Reverse x direction only for reordering
    !xp=-xp
    nface_inf=nface_inf+1
    ! This face becomes xmax in reordering
    elmt_face_pair(1,nface_inf)=ielmt
    elmt_face_pair(2,nface_inf)=iface ! see note on face ordering above
    elmt_face_pair(3,nface_inf)=xdir
    isinf=.true.
  elseif(nxmin.ne.0)then
    write(errtag,*)'ERROR: invalid number of face nodes on infinite element!', &
    nxmin
    return
  endif
  ! XMAX face
  if(nxmax.eq.4)then
    neface=neface+1

    ! extract active nodes
    inode=0
    do i_node=1,8
      if(isxmax(i_node))then
        inode=inode+1
        inodes(inode)=i_node
      endif
    enddo
    if(inode.ne.4)then
      write(errtag,*)'ERROR: number of face nodes must be 4!',inode
      return
    endif
   
    ! get face ID
    !iface=2 !get_faceidHEX8(inodes,ignodes)
    iface=get_faceidHEX8(inodes,ignodes)

    ! DO NOT reverse x direction
    nface_inf=nface_inf+1
    ! This remians xmax in reordering
    elmt_face_pair(1,nface_inf)=ielmt
    elmt_face_pair(2,nface_inf)=iface ! see note on face ordering above
    elmt_face_pair(3,nface_inf)=xdir
    isinf=.true.
  elseif(nxmax.ne.0)then
    write(errtag,*)'ERROR: invalid number of face nodes on infinite element!', &
    nxmax
    return
  endif
  
  ! YMIN face
  if(nymin.eq.4)then
    neface=neface+1

    ! extract active nodes
    inode=0
    do i_node=1,8
      if(isymin(i_node))then
        inode=inode+1
        inodes(inode)=i_node
      endif
    enddo
    if(inode.ne.4)then
      write(errtag,*)'ERROR: number of face nodes must be 4!',inode
      return
    endif
   
    ! get face ID
    !iface=1 !get_faceidHEX8(inodes,ignodes)
    iface=get_faceidHEX8(inodes,ignodes)

    !! Reverse y direction only for reordering
    !yp=-yp
    nface_inf=nface_inf+1
    ! This face becomes ymax in reordering
    elmt_face_pair(1,nface_inf)=ielmt
    elmt_face_pair(2,nface_inf)=iface ! see note on face ordering above
    elmt_face_pair(3,nface_inf)=ydir
    isinf=.true.
  elseif(nymin.ne.0)then
    write(errtag,*)'ERROR: invalid number of face nodes on infinite element!', &
    nymin
    return
  endif
  ! YMAX face
  if(nymax.eq.4)then
    neface=neface+1

    ! extract active nodes
    inode=0
    do i_node=1,8
      if(isymax(i_node))then
        inode=inode+1
        inodes(inode)=i_node
      endif
    enddo
    if(inode.ne.4)then
      write(errtag,*)'ERROR: number of face nodes must be 4!',inode
      return
    endif
   
    ! get face ID
    !iface=3 !get_faceidHEX8(inodes,ignodes)
    iface=get_faceidHEX8(inodes,ignodes)

    ! DO NOT reverse y direction
    nface_inf=nface_inf+1
    ! This remains ymax in reordering
    elmt_face_pair(1,nface_inf)=ielmt
    elmt_face_pair(2,nface_inf)=iface ! see note on face ordering above
    elmt_face_pair(3,nface_inf)=ydir
    isinf=.true.
  elseif(nymax.ne.0)then
    write(errtag,*)'ERROR: invalid number of face nodes on infinite element!', &
    nymax
    return
  endif
  
  ! ZMIN face
  if(nzmin.eq.4)then
    neface=neface+1

    ! extract active nodes
    inode=0
    do i_node=1,8
      if(iszmin(i_node))then
        inode=inode+1
        inodes(inode)=i_node
      endif
    enddo
    if(inode.ne.4)then
      write(errtag,*)'ERROR: number of face nodes must be 4!',inode
      return
    endif
   
    ! get face ID
    !iface=5 !get_faceidHEX8(inodes,ignodes)
    iface=get_faceidHEX8(inodes,ignodes)

    !! Reverse z direction only for reordering
    !zp=-zp
    nface_inf=nface_inf+1
    ! This face becomes zmax in reordering
    elmt_face_pair(1,nface_inf)=ielmt
    elmt_face_pair(2,nface_inf)=iface ! see note on face ordering above
    elmt_face_pair(3,nface_inf)=zdir
    isinf=.true.
  elseif(nzmin.ne.0)then
    write(errtag,*)'ERROR: invalid number of face nodes on infinite element!', &
    nzmin
    return
  endif
  ! ZMAX face
  if(zmaxinf)then
    if(nzmax.eq.4)then
      neface=neface+1

      ! extract active nodes
      inode=0
      do i_node=1,8
        if(iszmax(i_node))then
          inode=inode+1
          inodes(inode)=i_node
        endif
      enddo
      if(inode.ne.4)then
        write(errtag,*)'ERROR: number of face nodes must be 4!',inode
        return
      endif
   
      ! get face ID
      !iface=6 !get_faceidHEX8(inodes,ignodes)
      iface=get_faceidHEX8(inodes,ignodes)

      ! DO NOT reverse z direction
      nface_inf=nface_inf+1
      ! This face remains zmax in reordering
      elmt_face_pair(1,nface_inf)=ielmt
      elmt_face_pair(2,nface_inf)=iface ! see note on face ordering above
      elmt_face_pair(3,nface_inf)=zdir
      isinf=.true.
    elseif(nzmax.ne.0)then
      write(errtag,*)'ERROR: nvalid number of face nodes on infinite element!',&
      nzmax
      return
    endif
  endif

  if(neface==3)print*,'nface INF:',neface
  
  ! count assigned infinite elements
  if(isinf)then
    nelmt_inf=nelmt_inf+1
  endif 
  if(neface.gt.3)then
    write(errtag,*)'ERROR: NOT all faces of an element can be infinite!',neface
    return

  endif
  

!  print*,'old',gnum
!  write(*,*)(xp0(i),yp0(i),zp0(i),i=1,ngnode)
!
!  call sort_array_coord(ndim,ngnode,zp,yp,xp,gnum,ncount)
!  if(ngnode.ne.ncount)then
!    write(errtag,*)'ERROR: number of gnodes mismatched after sorting!',        &
!    ngnode,ncount,xmin,xmax,ymin,ymax,zmin,zmax
!    return
!  endif
!  print*,'New',(gnum(map2exodus_hex8(i)),i=1,8)
!  write(*,*)(xp(i),yp(i),zp(i),i=1,ngnode)
!  !!stop
!  !!print*,'-----------------------------------------------' 
!  ! Set new connectivity after coverting to EXODUS/CUBIT convention
!
!  g_num(:,ielmt)=gnum(map2exodus_hex8)
enddo

deallocate(ielmt_inf,g_numinf)

xnelmt_inf=nelmt_ing0-nelmt_inf
write(stdout,'(a,i0)')'total number of assigned infinite elements:',nelmt_inf
write(stdout,'(a,i0)')'total number of discarded infinite elements:',xnelmt_inf
write(stdout,'(a,i0)')'total number of assigned infinite faces:',nface_inf
if(xnelmt_inf.gt.0)then
  write(*,'(a,i0,a)')'WARNING: ',xnelmt_inf,' originally assigned &
  &infinite element found NOT on the boundary.'
  write(*,*)'These elements will become transition elements!'
endif

call save_infinite_elements()

!! format for elmt_face_pair
!ewidth=ceiling(log10(real(nelmt)+1.))
!write(format_ef,*)ewidth
!! faceID and infDIR are single digits
!format_ef='(i'//trim(adjustl(format_ef))//',1x,i1,1x,i2)'
!
!! Save infinite element file
!open(unit=17,file=trim(inp_path)//trim(inffile),status='replace',              &
!action='write',iostat=istat)
!if(istat/=0) then
!  print*,'ERROR: cannot open file:',trim(inffile)
!  stop
!endif
!
!write(17,'(i0,1x,g0.4)')200,0.0
!write(17,'(i0)')nface_inf
!write(17,format_ef)elmt_face_pair(:,1:nface_inf)
!close(17)
!deallocate(elmt_face_pair)

errcode=0
return

end subroutine reorder_infinite_elements
!===============================================================================

subroutine save_infinite_elements()
use global,only:nelmt,inffile,inp_path
implicit none
integer :: istat
! Format string
integer :: ewidth
character(len=20) :: format_ef

! format for elmt_face_pair
ewidth=ceiling(log10(real(nelmt)+1.))
write(format_ef,*)ewidth
! faceID and infDIR are single digits
format_ef='(i'//trim(adjustl(format_ef))//',1x,i1,1x,i2)'

! Save infinite element file
open(unit=17,file=trim(inp_path)//trim(inffile),status='replace',              &
action='write',iostat=istat)
if(istat/=0) then
  print*,'ERROR: cannot open file:',trim(inffile)
  stop
endif
write(17,'(i0,1x,g0.4)')200,0.0
write(17,'(i0)')nface_inf
write(17,format_ef)elmt_face_pair(:,1:nface_inf)
close(17)
deallocate(elmt_face_pair)

return

end subroutine save_infinite_elements
!===============================================================================

end module infinite_layer
!===============================================================================
