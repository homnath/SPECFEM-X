!WARNING: valid only for single fault
! DESCRIPTION
!  This module contain the fault preprocessing routines for split-node method.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO
!  - taper rom the edge of the element
module fault
use set_precision
use global,only:ndim
implicit none
character(len=250),private :: myfname=' => fault.f90'
character(len=500),private :: errsrc

! Fault surface: Plus side
! Fault element IDs
integer,allocatable :: elmt_pfault(:)
! Number of (hex)element faces which are quads
integer :: pfault_nface
! Element IDs
integer,allocatable :: pfault_ielmt(:)
! Face IDs
integer,allocatable :: pfault_iface(:)
! Edge IDs on the boundary. 0 indicates that the edge is not on the fault
! boundary
integer,allocatable :: pfault_iedge(:,:)
! slip vector
real(kind=kreal) :: pfault_svec(ndim)

! Following variables are determined in activaet_dof function in apply_bc.f90
integer :: pfault_nnode
integer,allocatable :: pfault_inode(:)

! Fault surface: Minus side
! Fault element IDs
integer,allocatable :: elmt_mfault(:)
! Number of (hex)element faces which are quads
integer :: mfault_nface
! Element IDs
integer,allocatable :: mfault_ielmt(:)
! Face IDs
integer,allocatable :: mfault_iface(:)
! Edge IDs on the boundary. 0 indicates that the edge is not on the fault
! boundary
integer,allocatable :: mfault_iedge(:,:)
! slip vector
real(kind=kreal) :: mfault_svec(ndim)

! Following variables are determined in activaet_dof function in apply_bc.f90
integer :: mfault_nnode
integer,allocatable :: mfault_inode(:)

contains
!-------------------------------------------------------------------------------

! This subroutine activates degrees of freedoms.
subroutine prepare_fault(errcode,errtag)
use global
use dimensionless,only:NONDIM_L
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: ios
integer :: i_face
integer :: count_fsurf
logical :: fsurf_stat,isempty
character(len=250) :: fname
character(len=250) :: fplus_file,fminus_file
character(len=150) :: data_path

errtag="ERROR: unknown!"
errcode=-1
! Set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

! Fault surface: Plus side
! Allocate and initialize fault element IDs
allocate(elmt_pfault(nelmt))
elmt_pfault=-1

fname=trim(data_path)//trim(faultslipfile_plus)//'_recreated'//trim(ptail_inp)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
! This is necessary for empty faultslipfile.
fsurf_stat=.true.
count_fsurf=0
isempty=.true.
fsurface_plus: do
  read(11,*,iostat=ios)pfault_svec
  if(ios/=0)exit fsurface_plus
  ! Nondimensionalise
  pfault_svec=NONDIM_L*pfault_svec
  count_fsurf=count_fsurf+1
  fsurf_stat=.false.

  read(11,*)pfault_nface
  allocate(pfault_ielmt(pfault_nface),pfault_iface(pfault_nface),              &
  pfault_iedge(4,pfault_nface))
  do i_face=1,pfault_nface
    read(11,*)pfault_ielmt(i_face),pfault_iface(i_face),pfault_iedge(:,i_face)
    elmt_pfault(pfault_ielmt(i_face))=i_face
  enddo
  fsurf_stat=.true.
  isempty=.false.
enddo fsurface_plus
close(11)
if(.not.fsurf_stat)then
  write(errtag,'(a)')'ERROR: some fault surfaces cannot be read for PLUS side!'
  return
endif

! Plot VTK file.
! File names for plotting fault slip VTK format.
fplus_file=trim(out_path)//trim(file_head)//'_fault_plus'//trim(ptail)//'.vtk'
call plot_fault_slip_vtk(pfault_nface,pfault_ielmt,pfault_iface,pfault_iedge,&
pfault_svec,fplus_file)

! Fault surface: Minus side
! Allocate and initialize fault element IDs
allocate(elmt_mfault(nelmt))
elmt_mfault=-1

fname=trim(data_path)//trim(faultslipfile_minus)//'_recreated'//trim(ptail_inp)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
! This is necessary for empty faultslipfile.
fsurf_stat=.true.
count_fsurf=0
isempty=.true.
fsurface_minus: do
  read(11,*,iostat=ios)mfault_svec
  if(ios/=0)exit fsurface_minus
  ! Nondimensionalise
  mfault_svec=NONDIM_L*mfault_svec
  count_fsurf=count_fsurf+1
  fsurf_stat=.false.

  read(11,*)mfault_nface
  allocate(mfault_ielmt(mfault_nface),mfault_iface(mfault_nface),              &
  mfault_iedge(4,mfault_nface))
  do i_face=1,mfault_nface
    read(11,*)mfault_ielmt(i_face),mfault_iface(i_face),mfault_iedge(:,i_face)
    elmt_mfault(mfault_ielmt(i_face))=i_face
  enddo
  fsurf_stat=.true.
  isempty=.false.
enddo fsurface_minus
close(11)
if(.not.fsurf_stat)then
  write(errtag,'(a)')'ERROR: some fault surfaces cannot be read for PLUS side!'
  return
endif

! Plot VTK file.
! File names for plotting fault slip VTK format.
fminus_file=trim(out_path)//trim(file_head)//'_fault_minus'//trim(ptail)//'.vtk'
call plot_fault_slip_vtk(mfault_nface,mfault_ielmt,mfault_iface,mfault_iedge,&
mfault_svec,fminus_file)

errcode=0

end subroutine prepare_fault
!===============================================================================

subroutine plot_fault_slip_vtk(nface,fault_ielmt,fault_iface,         &
fault_iedge,slip_vec,vtkout)
use global,only:myrank,ndim,maxngll2d,ngllx,nglly,g_num,g_coord,itaper_slip
use dimensionless,only:DIM_L
use math_library,only:i_uniinv
use element,only:hexface,hexface_edge
use math_constants,only:INFTOL,ZERO
implicit none
integer,intent(in) :: nface
integer,intent(in) :: fault_ielmt(:),fault_iface(:),fault_iedge(:,:)
real(kind=kreal),intent(in) :: slip_vec(ndim)
character(250),intent(in) :: vtkout

integer :: i,i_edge,i_face,i_gll,j
integer :: ielmt,iface,iedge(4)
integer :: nfgll
integer :: n1,n2
integer :: funit,nface_small,nsnode,nsnode_all
integer :: node_quad4(4)
integer,allocatable :: f_num(:,:),inode_order(:),nodelist(:)
real(kind=kreal),allocatable :: f_coord(:,:),f_slip(:,:),slip_face(:,:,:)
real(kind=kreal),allocatable :: slip_gll(:,:)

logical,allocatable :: isnode(:)

! WARNING: it works only for NGLLX=NGLLY=NGLLZ!!!
nfgll=maxngll2d

funit=101
open(funit,file=trim(vtkout),action='write',status='replace')
! allocate for nodelist and order
nsnode_all=maxngll2d*nface
allocate(nodelist(nsnode_all),inode_order(nsnode_all))
allocate(f_num(maxngll2d,nface),slip_face(NDIM,maxngll2d,nface))
! This has to be inside the loop but for ngllx=nglly=ngllz, we can put it here
allocate(slip_gll(NDIM,nfgll))

n1=1; n2=maxngll2d
do i_face=1,nface
  ielmt=fault_ielmt(i_face)
  iface=fault_iface(i_face)
  iedge=fault_iedge(:,i_face)
  
  ! Address for new connectivity
  f_num(:,i_face)= (/ (i, i=n1,n2) /)
  nodelist(n1:n2)=g_num(hexface(iface)%node,ielmt)
 
  ! Do NOT modify these parameters after the line below
  n1=n2+1; n2=n1+maxngll2d-1

  ! set slip on the GLL points  
  do i_gll=1,nfgll
    slip_gll(:,i_gll)=slip_vec
  enddo

  ! Set slip on the fault boundary nodes to ZERO
  if(itaper_slip.eq.1)then
    do i_edge=1,4
      if(iedge(i_edge).gt.0)then
        slip_gll(:,hexface_edge(iface,i_edge)%fnode)=ZERO
      endif
    enddo
  endif 

  slip_face(:,:,i_face)=slip_gll

enddo ! i_face
deallocate(slip_gll)

call i_uniinv(nodelist,inode_order)

nsnode=0

if(nsnode_all.gt.0)then
  nsnode=maxval(inode_order)
endif
! new connectivity
allocate(f_coord(ndim,nsnode),f_slip(NDIM,nsnode))
f_slip=INFTOL
do i_face=1,nface
  f_num(:,i_face)=inode_order(f_num(:,i_face))
  f_slip(:,f_num(:,i_face))=slip_face(:,:,i_face)
enddo
deallocate(slip_face) 

allocate(isnode(nsnode))
isnode=.false.

! assign f_coord
if(nsnode_all.gt.0)then
  f_coord(:,inode_order(1))=g_coord(:,nodelist(1))
  isnode(inode_order(1))=.true.
endif
do i=2,nsnode_all
  if(.not.isnode(inode_order(i)))then
     if(inode_order(i).le.0)then
       print*,'ERROR: "inode_order" value beyond range!',myrank,nsnode_all,    &
       inode_order(i)
       stop
     endif
     f_coord(:,inode_order(i))=g_coord(:,nodelist(i))
     isnode(inode_order(i))=.true.
  endif
enddo
deallocate(inode_order,isnode,nodelist)

! write VTK file
write(funit,'(a)')'# vtk DataFile Version 2.0'
write(funit,'(a)')'Unstructured Grid Example'
write(funit,'(a)')'ASCII'
write(funit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(funit,'(a,i5,a)')'POINTS',nsnode,' float'
do i=1,nsnode
  write(funit,'(3(f14.6,1x))')DIM_L*f_coord(:,i)
enddo
write(funit,*)
nface_small=nface*(ngllx-1)*(nglly-1) ! ngllx = nglly
write(funit,'(a,i5,a,i5)')'CELLS',nface_small,' ',5*nface_small
do i_face=1,nface
  do j=1,nglly-1
    do i=1,ngllx-1
      node_quad4(1)=(j-1)*ngllx+i
      node_quad4(2)=node_quad4(1)+1

      node_quad4(4)=node_quad4(1)+ngllx ! G4 is located NGLLX above G1
      node_quad4(3)=node_quad4(4)+1     ! G3 is located 1 of G4

      ! VTK format indexing starts from 0
      write(funit,'(5(i5,1x))')4,f_num(node_quad4,i_face)-1
    enddo
  enddo
enddo
write(funit,*)
write(funit,'(a,i5)')'CELL_TYPES',nface_small
do i=1,nface_small
  write(funit,'(i1)')9
enddo
write(funit,*)
write(funit,'(a,1x,i5)')'POINT_DATA',nsnode
write(funit,'(a)')'VECTORS slip float'
do i=1,nsnode
  write(funit,'(3(e13.6,1x))')DIM_L*f_slip(:,i)
enddo
close(funit)
deallocate(f_coord,f_num,f_slip)

close(101)
end subroutine plot_fault_slip_vtk
!===============================================================================

! This routine read and applies the fault slip specified in the faultslip file
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
subroutine compute_fault_slip_load(plus_or_minus,sfac,storekmat,load,          &
                                     errcode,errtag)
use global
use math_constants
use element,only:hexface,hexface_edge!,hexface_sign
!use preprocess
!use dimensionless,only:NONDIM_L,DIM_L
implicit none
integer,intent(in) :: plus_or_minus
real(kind=kreal),intent(in) :: sfac
real(kind=kreal),intent(in) :: storekmat(:,:,:)
real(kind=kreal),intent(inout) :: load(0:neq)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i,i_face,i_gll
integer :: i_edge,iedge(4),iedof,ielmt,iface
integer,allocatable :: egdof(:)

integer :: fault_nface
real(kind=kreal) :: slip_vec(NDIM)
real(kind=kreal),allocatable :: fault_ielmt(:),fault_iface(:),fault_iedge(:,:)

real(kind=kreal),allocatable :: slip_gll(:,:)
real(kind=kreal),allocatable :: kmat(:,:)

errtag="ERROR: unknown!"
errcode=-1

! Assign appropriate fault parameters
if(plus_or_minus.gt.0)then
  slip_vec=pfault_svec
  fault_nface=pfault_nface
else
  slip_vec=mfault_svec
  fault_nface=mfault_nface
endif
allocate(fault_ielmt(fault_nface))
allocate(fault_iface(fault_nface))
allocate(fault_iedge(4,fault_nface))
if(fault_nface.gt.0)then
  if(plus_or_minus.gt.0)then
    fault_ielmt=pfault_ielmt
    fault_iface=pfault_iface
    fault_iedge=pfault_iedge
  else
    fault_ielmt=mfault_ielmt
    fault_iface=mfault_iface
    fault_iedge=mfault_iedge
  endif
endif

!print*,'slip_vec:',slip_vec
!allocate(kmat(nedof,nedof))
allocate(kmat(nedofu,nedofu))
allocate(slip_gll(NDIM,ngll))
!allocate(egdof(nedof))
allocate(egdof(nedofu))
kmat=ZERO
do i_face=1,fault_nface
  ielmt=fault_ielmt(i_face)
  iface=fault_iface(i_face)
  iedge=fault_iedge(:,i_face)
  
  ! set slip on the GLL points  
  slip_gll=ZERO
  do i_gll=1,maxngll2d != nfgll=ngllxy=ngllyz=ngllzx
    slip_gll(:,hexface(iface)%node(i_gll))=slip_vec
  enddo
  slip_gll=sfac*slip_gll

  ! Following will set the slip to ZERO on the fault boundary
  ! set slip on the fault boundary nodes to ZERO
  if(itaper_slip.eq.1)then
    do i_edge=1,4
      if(iedge(i_edge).gt.0)then
        slip_gll(:,hexface_edge(iface,i_edge)%node)=ZERO
      endif
    enddo
  endif
   
  !egdof=gdof_elmt(:,ielmt)
  egdof=gdof_elmt(edofu,ielmt)
 
  !kmat=storekmat(:,:,ielmt)
  kmat=storekmat(edofu,edofu,ielmt)
  iedof=0
  do i_gll=1,ngll
    do i=1,nndofu
      iedof=iedof+1
      if(slip_gll(i,i_gll)/=zero)then
        load(egdof)=load(egdof)-plus_or_minus*kmat(:,iedof)*slip_gll(i,i_gll)
      endif
    enddo
  enddo
enddo

deallocate(kmat)
deallocate(slip_gll)
deallocate(egdof)
if(allocated(fault_ielmt))deallocate(fault_ielmt)
if(allocated(fault_iface))deallocate(fault_iface)
if(allocated(fault_iedge))deallocate(fault_iedge)

errcode=0

end subroutine compute_fault_slip_load
!===============================================================================

subroutine cleanup_fault
implicit none
if(allocated(elmt_pfault ))deallocate(elmt_pfault )
if(allocated(pfault_ielmt))deallocate(pfault_ielmt)
if(allocated(pfault_iface))deallocate(pfault_iface)
if(allocated(pfault_iedge))deallocate(pfault_iedge)
if(allocated(pfault_inode))deallocate(pfault_inode)
if(allocated(elmt_mfault ))deallocate(elmt_mfault )
if(allocated(mfault_ielmt))deallocate(mfault_ielmt)
if(allocated(mfault_iface))deallocate(mfault_iface)
if(allocated(mfault_iedge))deallocate(mfault_iedge)
if(allocated(mfault_inode))deallocate(mfault_inode)

end subroutine cleanup_fault
!===============================================================================

end module fault
!===============================================================================
