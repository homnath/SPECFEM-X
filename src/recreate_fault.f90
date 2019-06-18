! This module contain the routines related to earthquake fault.
! REVISION
!   HNG, Oct 12,2017
module recreate_fault
use set_precision
implicit none

contains
!-------------------------------------------------------------------------------

! This subroutine geometrically splits fault nodes. 
! This subroutine is used by a single/serial processor.
! This subroutine is called from read_input if the read_input is 
! called by serial programs such as partmesh and specfemx.
subroutine split_fault_nodes(fsfile_plus,fsfile_minus)
use global,only:g_num,inp_path,NDIM
use math_library,only:i_uniinv
implicit none
character(len=250),intent(in) :: fsfile_plus,fsfile_minus
character(len=250) :: newfsfile
character(len=60) :: format_str1

integer :: fs_nelmt ! number of BC elements
! fist row = element ID, second row = entity ID
integer,dimension(:,:),allocatable :: fs_elmt
integer :: ios,i_edge,i_elmt,nfs,istat

integer :: ielmt,iface
integer :: n1,n2,nsnode,nsnode_all,sum_eval
integer,dimension(4) :: inodes
integer,dimension(6,4) :: node_face ! local node numbers in each face
integer,dimension(4,2) :: node_fedge ! local node numbers in each face edge

integer,allocatable :: evalence(:),inode_felmt(:,:),inode_order(:),nodelist(:)
integer,allocatable :: isbound_edge(:,:)
! Fault slip vector
real(kind=kreal) :: fsvec(NDIM) ! replace this with NDIM
integer,allocatable :: temp_mat(:,:)

! Local node numbering in each face CUBIT/EXODUS convention
node_face(1,:)=(/1,2,6,5/) ! counterclockwise w.r.t. outer normal
node_face(2,:)=(/2,3,7,6/) ! counterclockwise w.r.t. outer normal
node_face(3,:)=(/4,3,7,8/) ! clockwise w.r.t. outer normal
node_face(4,:)=(/1,4,8,5/) ! clockwise w.r.t. outer normal
node_face(5,:)=(/1,2,3,4/) ! clockwise w.r.t. outer normal
node_face(6,:)=(/5,6,7,8/) ! counterclockwise w.r.t. outer normal

! Local node numbering in each face edge CUBIT/EXODUS convention
node_fedge(1,:)=(/1,2/)
node_fedge(2,:)=(/2,3/)
node_fedge(3,:)=(/3,4/)
node_fedge(4,:)=(/4,1/)

! Open fault surface file
open(unit=16,file=trim(inp_path)//trim(fsfile_plus),status='old',action='read',     &
iostat=istat)
if(istat/=0) then
  print*,'ERROR: cannot open file:',trim(fsfile_plus)
  stop
endif

! New faultslip file
newfsfile=trim(inp_path)//trim(fsfile_plus)//'_recreated'
open(unit=17,file=trim(newfsfile),status='replace',action='write',&
iostat=istat)
if(istat/=0) then
  print*,'ERROR: cannot open file:',trim(newfsfile)
  stop
endif

nfs=0
fsurface: do ! nfs=1,nfs
  read(16,*,iostat=ios)fsvec
  if(ios/=0)then
    exit fsurface
  endif
  nfs=nfs+1
  read(16,*)fs_nelmt
  if(fs_nelmt>0)then
    allocate(fs_elmt(2,fs_nelmt))

    read(16,*,iostat=istat)fs_elmt
    if(istat/=0)then
      write(*,*)'ERROR: cannot read complete fault slip file: ',trim(fsfile_plus),'!'
      stop
    endif
  endif

  !-----------------------------------------------------------------------------
  ! Determine boundary nodes on the fault.
  ! Allocate for nodelist and order.
  nsnode_all=4*fs_nelmt
  allocate(nodelist(nsnode_all),inode_order(nsnode_all))
  nodelist=-1
  n1=1; n2=4
  do i_elmt=1,fs_nelmt
    ielmt=fs_elmt(1,i_elmt)
    iface=fs_elmt(2,i_elmt)
    ! NOTE: surprisingly intel Fortran 19 compiler gives the strange error 
    ! for a following statement.
    ! nodelist(n1:n2)=g_num(node_face(iface,:),ielmt)
    ! Therefore, we use two statements.
    inodes=node_face(iface,:)
    nodelist(n1:n2)=g_num(inodes,ielmt)
    n1=n2+1; n2=n1+3
  enddo
  if(any(nodelist.lt.0))then
    write(*,*)'ERROR: nodelist array is not filled!'
    stop
  endif
  call i_uniinv(nodelist,inode_order)
  deallocate(nodelist)

  nsnode=maxval(inode_order)

  allocate(evalence(nsnode),inode_felmt(4,fs_nelmt))

  evalence=0
  
  n1=1; n2=4
  do i_elmt=1,fs_nelmt
    ielmt=fs_elmt(1,i_elmt)
    iface=fs_elmt(2,i_elmt)
    inode_felmt(:,i_elmt)=inode_order(n1:n2)
    evalence(inode_felmt(:,i_elmt))=evalence(inode_felmt(:,i_elmt))+1
    n1=n2+1; n2=n1+3
  enddo
  deallocate(inode_order)
  print*,'element valency of fault node MIN MAX:',minval(evalence),maxval(evalence)

  ! Boundary edges
  allocate(isbound_edge(4,fs_nelmt))
  isbound_edge=0 !.false.
  do i_elmt=1,fs_nelmt
    ielmt=fs_elmt(1,i_elmt)
    iface=fs_elmt(2,i_elmt)
    do i_edge=1,4
      sum_eval=evalence(inode_felmt(node_fedge(i_edge,1),i_elmt)) + &
               evalence(inode_felmt(node_fedge(i_edge,2),i_elmt))
      if(sum_eval.lt.6)isbound_edge(i_edge,i_elmt)=1 !.true.
    enddo
  enddo

  deallocate(evalence,inode_felmt)
  !-----------------------------------------------------------------------------

  ! Format string for element ID and face ID
  write(format_str1,*)ceiling(log10(real(maxval(fs_elmt(1,:)))+1.))
  format_str1='(i'//trim(adjustl(format_str1))//',1x,I2,1x,I2,1x,I2,1x,I2,1x,I2)'
  ! i2 is sufficient for face id or node id

  ! Write new fault information including boundary edges
  write(17,*)fsvec
  write(17,*)fs_nelmt
  allocate(temp_mat(6,fs_nelmt))
  temp_mat(1:2,:)=fs_elmt
  temp_mat(3:6,:)=isbound_edge
  write(17,format_str1)temp_mat!fs_elmt(:,mpart(i_part)%iloc)
  deallocate(temp_mat)
  deallocate(fs_elmt,isbound_edge)
enddo fsurface
close(16)
close(17)
if(nfs.le.0)then
  write(*,*)'ERROR: invalid or empty faultslip file!'
  stop
endif

end subroutine split_fault_nodes
!===============================================================================

! This subroutine recreates fault slip file with added information on the face
! edges on the fault boundary.
! This subroutine is used by a single/serial processor.
! This subroutine is called from read_input if the read_input called by serial
! program such as partmesh library and semvisco.
subroutine recreate_faultslip_file(fsfile)
use global,only:g_num,inp_path,NDIM
use math_library,only:i_uniinv
implicit none
character(len=250),intent(in) :: fsfile
character(len=250) :: newfsfile
character(len=60) :: format_str1

integer :: fs_nelmt ! number of BC elements
! fist row = element ID, second row = entity ID
integer,dimension(:,:),allocatable :: fs_elmt
integer :: ios,i_edge,i_elmt,nfs,istat

integer :: ielmt,iface
integer :: n1,n2,nsnode,nsnode_all,sum_eval
integer,dimension(4) :: inodes
integer,dimension(6,4) :: node_face ! local node numbers in each face
integer,dimension(4,2) :: node_fedge ! local node numbers in each face edge

integer,allocatable :: evalence(:),inode_felmt(:,:),inode_order(:),nodelist(:)
integer,allocatable :: isbound_edge(:,:)
! fault slip vector
real(kind=kreal) :: fsvec(NDIM) ! replace this with NDIM
integer,allocatable :: temp_mat(:,:)

! local node numbering in each face CUBIT/EXODUS convention
node_face(1,:)=(/1,2,6,5/) ! counterclockwise w.r.t. outer normal
node_face(2,:)=(/2,3,7,6/) ! counterclockwise w.r.t. outer normal
node_face(3,:)=(/4,3,7,8/) ! clockwise w.r.t. outer normal
node_face(4,:)=(/1,4,8,5/) ! clockwise w.r.t. outer normal
node_face(5,:)=(/1,2,3,4/) ! clockwise w.r.t. outer normal
node_face(6,:)=(/5,6,7,8/) ! counterclockwise w.r.t. outer normal

! local node numbering in each face edge CUBIT/EXODUS convention
node_fedge(1,:)=(/1,2/)
node_fedge(2,:)=(/2,3/)
node_fedge(3,:)=(/3,4/)
node_fedge(4,:)=(/4,1/)

! open fs file
open(unit=16,file=trim(inp_path)//trim(fsfile),status='old',action='read',     &
iostat=istat)
if(istat/=0) then
  print*,'ERROR: cannot open file:',trim(fsfile)
  stop
endif

! new faultslip file
newfsfile=trim(inp_path)//trim(fsfile)//'_recreated'
open(unit=17,file=trim(newfsfile),status='replace',action='write',&
iostat=istat)
if(istat/=0) then
  print*,'ERROR: cannot open file:',trim(newfsfile)
  stop
endif

nfs=0
fsurface: do ! nfs=1,nfs
  read(16,*,iostat=ios)fsvec
  if(ios/=0)then
    exit fsurface
  endif
  nfs=nfs+1
  read(16,*)fs_nelmt
  if(fs_nelmt>0)then
    allocate(fs_elmt(2,fs_nelmt))

    read(16,*,iostat=istat)fs_elmt
    if(istat/=0)then
      write(*,*)'ERROR: cannot read complete fault slip file: ',trim(fsfile),'!'
      stop
    endif
  endif

  !-----------------------------------------------------------------------------
  ! Determine boundary nodes on the fault.
  ! Allocate for nodelist and order.
  nsnode_all=4*fs_nelmt
  allocate(nodelist(nsnode_all),inode_order(nsnode_all))
  nodelist=-1
  n1=1; n2=4
  do i_elmt=1,fs_nelmt
    ielmt=fs_elmt(1,i_elmt)
    iface=fs_elmt(2,i_elmt)
    ! NOTE: surprisingly intel Fortran 19 compiler gives the strange error 
    ! for a following statement.
    !nodelist(n1:n2)=g_num(node_face(iface,:),ielmt)
    ! Therefor, we use two statements.
    inodes=node_face(iface,:)
    nodelist(n1:n2)=g_num(inodes,ielmt)
    n1=n2+1; n2=n1+3
  enddo
  if(any(nodelist.lt.0))then
    write(*,*)'ERROR: nodelist array is not filled!'
    stop
  endif
  call i_uniinv(nodelist,inode_order)
  deallocate(nodelist)

  nsnode=maxval(inode_order)

  allocate(evalence(nsnode),inode_felmt(4,fs_nelmt))

  evalence=0
  
  n1=1; n2=4
  do i_elmt=1,fs_nelmt
    ielmt=fs_elmt(1,i_elmt)
    iface=fs_elmt(2,i_elmt)
    inode_felmt(:,i_elmt)=inode_order(n1:n2)
    evalence(inode_felmt(:,i_elmt))=evalence(inode_felmt(:,i_elmt))+1
    n1=n2+1; n2=n1+3
  enddo
  deallocate(inode_order)
  print*,'element valency of fault node MIN MAX:',minval(evalence),maxval(evalence)

  ! Boundary edges
  allocate(isbound_edge(4,fs_nelmt))
  isbound_edge=0 !.false.
  do i_elmt=1,fs_nelmt
    ielmt=fs_elmt(1,i_elmt)
    iface=fs_elmt(2,i_elmt)
    do i_edge=1,4
      sum_eval=evalence(inode_felmt(node_fedge(i_edge,1),i_elmt)) + &
               evalence(inode_felmt(node_fedge(i_edge,2),i_elmt))
      if(sum_eval.lt.6)isbound_edge(i_edge,i_elmt)=1 !.true.
    enddo
  enddo

  deallocate(evalence,inode_felmt)
  !-----------------------------------------------------------------------------

  ! Format string for element ID and face ID
  write(format_str1,*)ceiling(log10(real(maxval(fs_elmt(1,:)))+1.))
  format_str1='(i'//trim(adjustl(format_str1))//',1x,I2,1x,I2,1x,I2,1x,I2,1x,I2)'
  ! i2 is sufficient for face id or node id

  ! Write new fault information including boundary edges
  write(17,*)fsvec
  write(17,*)fs_nelmt
  allocate(temp_mat(6,fs_nelmt))
  temp_mat(1:2,:)=fs_elmt
  temp_mat(3:6,:)=isbound_edge
  write(17,format_str1)temp_mat!fs_elmt(:,mpart(i_part)%iloc)
  deallocate(temp_mat)
  deallocate(fs_elmt,isbound_edge)
enddo fsurface
close(16)
close(17)
if(nfs.le.0)then
  write(*,*)'ERROR: invalid or empty faultslip file!'
  stop
endif

end subroutine recreate_faultslip_file
!===============================================================================

end module recreate_fault
!===============================================================================
