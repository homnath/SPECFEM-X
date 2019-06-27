! this module contains the routines to process parallel communications across
! the ghost partitions
! TODO:
!   - better to avoid using allocatable component of derived type variable which
!     is fortran 95 standard? How much influence will be on the performance with
!     repeated allocate and deallocate
! REVISION:
!   HNG, APR 15,2011; HNG, Apr 09,2010
module ghost_library_mpi
use set_precision_mpi
use mpi_library
use global,only:absmaxcoord,myrank,nproc
integer :: ngpart,maxngnode
! ghost partitions
type ghost_partition
  integer :: rank,nnode
  integer,dimension(:),allocatable :: node,gdof
  logical,dimension(:),allocatable :: isnode
  ! whether the node on the interface is intact (.true.) or void (.false.)
end type ghost_partition
type(ghost_partition),dimension(:),allocatable :: gpart

contains

!-----------------------------------------------------------
subroutine prepare_ghost()
use element
use global,only:ndim,nnode,nndof,ngllx,nglly,ngllz,g_num0,g_num,g_coord,gfile, &
part_path,proc_str,file_head,stdout
use math_library, only : iquick_sort,sort_array_coord
use math_library_mpi, only : maxscal

implicit none
integer :: istat

integer,dimension(8) :: ign,jgn,kgn ! ith, jth, and kth GLL indices of node
integer :: ig0,ig1,jg0,jg1,kg0,kg1
integer :: i_g,j_g,k_g
!integer,dimension(6,4) :: node_face ! local node numbers in each face
!integer,dimension(12,2) :: node_edge ! local node numbers in each edge

character(len=20) :: format_str
character(len=250) :: fname,ofname

integer :: mrank,maxngpart ! partition ID
integer :: i_elmt,i_gpart,grank,melmt,ngelmt,igllp,inode,maxngelmt
integer :: etype,eid
integer :: ncount,new_ncount,ngllxy,maxngll_face
logical,dimension(nnode) :: switch_node
integer,dimension(:),allocatable :: itmp_array ! temporary integer array
double precision,dimension(:),allocatable :: xp,yp,zp
integer,dimension(3) :: ngll_vec ! (/ngllx,nglly,ngllz/) in ascending order
integer :: i
integer :: inode1,inode2(2),inode4(4),inodes(4)
integer :: v2(2),v4(4),gnode8(8)

character(len=60) :: sline
! sline must be 60 characeters long
character(len=250) :: errtag
integer :: errcode

errtag=""; errcode=-1

! find maximum ngll points in a face
! this is needed only for memory allocation
ngll_vec=(/ ngllx,nglly,ngllz /)
ngll_vec=iquick_sort(ngll_vec,NDIM)
maxngll_face=ngll_vec(NDIM)*ngll_vec(NDIM-1)

!! local node numbering in each face CUBIT/EXODUS convention
!node_face(1,:)=(/1,2,6,5/) ! front
!node_face(2,:)=(/2,3,7,6/) ! right
!node_face(3,:)=(/4,3,7,8/) ! back
!node_face(4,:)=(/1,4,8,5/) ! left
!node_face(5,:)=(/1,2,3,4/) ! bottom
!node_face(6,:)=(/5,6,7,8/) ! top
!
!! local node numbering in each edge CUBIT/EXODUS convention
!! bottom edges
!node_edge(1,:)=(/1,2/);  node_edge(2,:)=(/2,3/)
!node_edge(3,:)=(/3,4/);  node_edge(4,:)=(/4,1/)
!! side edges
!node_edge(5,:)=(/1,5/);  node_edge(6,:)=(/2,6/)
!node_edge(7,:)=(/3,7/);  node_edge(8,:)=(/4,8/)
!! top edges
!node_edge(9,:)=(/5,6/);  node_edge(10,:)=(/6,7/)
!node_edge(11,:)=(/7,8/); node_edge(12,:)=(/8,5/)

! GLL indices for nodes of a hexahedral mesh
! ith
ign(1)=1; ign(4)=1; ign(5)=1; ign(8)=1
ign(2)=ngllx; ign(3)=ngllx; ign(6)=ngllx; ign(7)=ngllx
! jth
jgn(1)=1; jgn(2)=1; jgn(5)=1; jgn(6)=1
jgn(4)=nglly; jgn(3)=nglly; jgn(8)=nglly; jgn(7)=nglly
! kth
kgn(1)=1; kgn(2)=1; kgn(4)=1; kgn(3)=1
kgn(5)=ngllz; kgn(6)=ngllz; kgn(8)=ngllz; kgn(7)=ngllz

ngllxy=ngllx*nglly

! open appropriate ghost file
write(format_str,*)ceiling(log10(real(nproc)+1))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
write(fname, fmt=format_str)trim(part_path)//trim(gfile)//'_proc',myrank
open(unit=11,file=trim(fname),access='stream',form='unformatted',    &
status='old',action='read',iostat=istat)
if (istat /= 0)then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  call control_error(errcode,errtag,stdout,myrank)
endif
read(11)sline
read(11)mrank ! master partition ID

if(mrank/=myrank)then
  write(errtag,*)'ERROR: wrong gpart file partition ',mrank,' !'
  call control_error(errcode,errtag,stdout,myrank)
endif

read(11)sline
read(11)ngpart

! allocate gpart
allocate(gpart(ngpart))
gpart(1:ngpart)%nnode=0
gpart(1:ngpart)%rank=-1
ofname='tmp/'//trim(file_head)//'_partitioninfo'//trim(adjustl(proc_str))
open(unit=22,file=trim(ofname),access='stream',form='unformatted',             &
status='replace',action='write',iostat=istat)
if (istat /= 0)then
  write(*,'(a)')'ERROR: output file "'//trim(fname)//'" cannot be opened!'
  stop
endif
write(22)nnode
write(22)ngpart

do i_gpart=1,ngpart ! ghost partitions loop
  read(11)sline
  read(11)grank
  write(22)grank

  read(11)sline
  read(11)ngelmt
  read(11)sline
  allocate(itmp_array(ngelmt*maxngll_face))
  itmp_array=-1
  switch_node=.false.
  ncount=0
  do i_elmt=1,ngelmt ! ghost elements loop
    read(11)melmt,etype,inodes
    ! Determine the entity ID
    ! We need to use old connectivity because the node ID in the ghost files are
    ! for the old connectivity.
    gnode8=g_num0(:,melmt)
    if(etype==1)then
      inode1=inodes(1)
      eid=get_nodeidHEX8(inode1,gnode8)
    elseif(etype==2)then
      inode2=inodes(1:etype)
      eid=get_edgeidHEX8(inode2,gnode8)
    elseif(etype==4)then
      inode4=inodes(1:etype)
      eid=get_faceidHEX8(inode4,gnode8)
    else
      write(errtag,*)'ERROR: wrong etype:',etype,' for ghost partition ',mrank,'!'
      call control_error(errcode,errtag,stdout,myrank)
    endif

    ! initialize
    ig0=-1; ig1=-1
    jg0=-1; jg1=-1
    kg0=-1; kg1=-1
    ! find range of GLL indices
    if(etype==1)then ! node
      ig0=ign(eid); ig1=ign(eid)
      jg0=jgn(eid); jg1=jgn(eid)
      kg0=kgn(eid); kg1=kgn(eid)
    elseif(etype==2)then ! edge
      ig0=minval(ign(node_edge(:,eid))); ig1=maxval(ign(node_edge(:,eid)))
      jg0=minval(jgn(node_edge(:,eid))); jg1=maxval(jgn(node_edge(:,eid)))
      kg0=minval(kgn(node_edge(:,eid))); kg1=maxval(kgn(node_edge(:,eid)))
    elseif(etype==4)then ! face
      ig0=minval(ign(node_face(:,eid))); ig1=maxval(ign(node_face(:,eid)))
      jg0=minval(jgn(node_face(:,eid))); jg1=maxval(jgn(node_face(:,eid)))
      kg0=minval(kgn(node_face(:,eid))); kg1=maxval(kgn(node_face(:,eid)))
    else
      write(errtag,*)'ERROR: wrong etype:',etype,' for ghost partition ',mrank,'!'
      call control_error(errcode,errtag,stdout,myrank)
    endif

    do k_g=kg0,kg1
      do j_g=jg0,jg1
        do i_g=ig0,ig1
          igllp=(k_g-1)*ngllxy+(j_g-1)*ngllx+i_g
          inode=g_num(igllp,melmt)
          if(.not.switch_node(inode))then
            ncount=ncount+1
            itmp_array(ncount)=inode
            switch_node(inode)=.true.
          endif
        enddo
      enddo
    enddo
  enddo ! do i_elmt

  gpart(i_gpart)%nnode=ncount
  gpart(i_gpart)%rank=grank
  allocate(gpart(i_gpart)%node(ncount),gpart(i_gpart)%isnode(ncount))
  gpart(i_gpart)%isnode=.true. ! all nodes are active
  allocate(gpart(i_gpart)%gdof(ncount*nndof))

  gpart(i_gpart)%node=itmp_array(1:ncount)

  !order nodal array to match with ghost partitions
  !extract coordinates
  allocate(xp(ncount),yp(ncount),zp(ncount))

  !xp=g_coord(1,itmp_array(1:ncount))
  !yp=g_coord(2,itmp_array(1:ncount))
  !zp=g_coord(3,itmp_array(1:ncount))
  do i=1,ncount
    xp(i)=g_coord(1,itmp_array(i))
    yp(i)=g_coord(2,itmp_array(i))
    zp(i)=g_coord(3,itmp_array(i))
  enddo
  deallocate(itmp_array)

  !sort nodes
  call sort_array_coord(ndim,ncount,xp,yp,zp,gpart(i_gpart)%node,new_ncount)
  deallocate(xp,yp,zp)
  if(ncount/=new_ncount)then
    write(errtag,*)'ERROR: number of ghost nodes mismatched after sorting!'
    call control_error(errcode,errtag,stdout,myrank)
  endif

  ! find ghost gdof
  !gpart(i_gpart)%gdof=reshape(gdof(:,gpart(i_gpart)%node),(/ncount*nndof/))
  write(22)gpart(i_gpart)%nnode
  write(22)gpart(i_gpart)%node

enddo ! do i_gpart
close(22)
close(11)
!stop
maxngnode=maxscal(maxval(gpart(1:ngpart)%nnode))

errcode=0

call sync_process
return
end subroutine prepare_ghost
!========i======================================================================

! this subroutine prepare GDOFs. note that the variable gpart(i_gpart)%gdof is
! allocated in the previous routine
subroutine prepare_ghost_gdof()
use global,only:gdof,nnode,nndof

implicit none

integer :: i_gpart

do i_gpart=1,ngpart ! ghost partitions loop
  ! find ghost gdof
  gpart(i_gpart)%gdof=reshape(gdof(:,gpart(i_gpart)%node),(/gpart(i_gpart)%nnode*nndof/))
enddo ! do i_gpart
return
end subroutine prepare_ghost_gdof
!===============================================================================

! modify ghost gdof based on the modified gdof
subroutine modify_ghost(isnode)
use global,only:gdof,nnode,nndof

implicit none
logical,intent(in) :: isnode(nnode)

integer :: i,i_gpart,ncount

! modify ghost gdof based on the modified gdof
do i_gpart=1,ngpart
  ncount=gpart(i_gpart)%nnode
  gpart(i_gpart)%gdof=reshape(gdof(:,gpart(i_gpart)%node),(/ncount*nndof/))
  do i=1,ncount
    gpart(i_gpart)%isnode(i)=isnode(gpart(i_gpart)%node(i))
  enddo
enddo ! do i_gpart
return
end subroutine modify_ghost

!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine send_ghosts_for_gpu_mode(nndof,neq,send_array,recv_array)
!use math_library, only : maxscal_par
use mpi
implicit none
integer,intent(in) :: nndof,neq
real(kind=kreal),dimension(nndof*maxngnode,ngpart), intent(in) :: send_array
real(kind=kreal),dimension(nndof*maxngnode,ngpart), intent(out) :: recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
real(kind=kreal),parameter :: zero=0.0_kreal
integer :: ierr,i_gpart,ncount

do i_gpart=1,ngpart
  ncount=gpart(i_gpart)%nnode*nndof
  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo
! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine send_ghosts_for_GPU_mode


!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine assemble_ghosts(nndof,neq,array,array_g)
!use math_library, only : maxscal_par
use mpi
implicit none
integer,intent(in) :: nndof,neq
real(kind=kreal),dimension(0:neq),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g
real(kind=kreal),dimension(nndof*maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
real(kind=kreal),parameter :: zero=0.0_kreal

integer :: ierr,i_gpart,ncount

array_g=array
send_array=zero; recv_array=zero
do i_gpart=1,ngpart
  ncount=gpart(i_gpart)%nnode*nndof
  ! array to send
  send_array(1:ncount,i_gpart)=array(gpart(i_gpart)%gdof)
  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ncount=gpart(i_gpart)%nnode*nndof
  array_g(gpart(i_gpart)%gdof)=array_g(gpart(i_gpart)%gdof)+ &
  recv_array(1:ncount,i_gpart)
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
array_g(0)=zero
return
end subroutine assemble_ghosts
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations
subroutine assemble_ghosts_nodal(nndof,array,array_g)
use mpi
use global,only:nnode
implicit none
integer,intent(in) :: nndof
real(kind=kreal),dimension(nndof,nnode),intent(in) :: array
real(kind=kreal),dimension(nndof,nnode),intent(out) :: array_g

!logical,dimension(maxngnode,ngpart) :: lsend_array,lrecv_array
real(kind=kreal),dimension(nndof*maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
!integer,dimension(nnode) :: ngpart_node ! number of ghost partition for a node
real(kind=kreal),parameter :: zero=0.0_kreal
integer :: ierr,i_gpart,ncount,ngnode

array_g=array
send_array=zero; recv_array=zero
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  ! store array-to-send in a garray
  send_array(1:ncount,i_gpart)=reshape(array(:,gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  array_g(:,gpart(i_gpart)%node)=array_g(:,gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/nndof,ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine assemble_ghosts_nodal
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations for the n-component vector.
subroutine assemble_ghosts_nodal_vector(array,array_g)
use mpi
use global,only:NDIM,nnode
implicit none
real(kind=kreal),dimension(NDIM,nnode),intent(in) :: array
real(kind=kreal),dimension(NDIM,nnode),intent(out) :: array_g

!logical,dimension(maxngnode,ngpart) :: lsend_array,lrecv_array
real(kind=kreal),dimension(NDIM*maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
!integer,dimension(nnode) :: ngpart_node ! number of ghost partition for a node
real(kind=kreal),parameter :: zero=0.0_kreal
integer :: ierr,i_gpart,ncount,ngnode

array_g=array
send_array=zero; recv_array=zero
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*NDIM
  ! store array-to-send in a garray
  send_array(1:ncount,i_gpart)=reshape(array(:,gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*NDIM
  array_g(:,gpart(i_gpart)%node)=array_g(:,gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/NDIM,ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine assemble_ghosts_nodal_vector
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations for the n-component vector.
subroutine assemble_ghosts_nodal_vectorn(ncomp,array,array_g)
use mpi
use global,only:nnode
implicit none
integer,intent(in) :: ncomp
real(kind=kreal),dimension(ncomp,nnode),intent(in) :: array
real(kind=kreal),dimension(ncomp,nnode),intent(out) :: array_g

!logical,dimension(maxngnode,ngpart) :: lsend_array,lrecv_array
real(kind=kreal),dimension(ncomp*maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
!integer,dimension(nnode) :: ngpart_node ! number of ghost partition for a node
real(kind=kreal),parameter :: zero=0.0_kreal
integer :: ierr,i_gpart,ncount,ngnode

array_g=array
send_array=zero; recv_array=zero
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*ncomp
  ! store array-to-send in a garray
  send_array(1:ncount,i_gpart)=reshape(array(:,gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*ncomp
  array_g(:,gpart(i_gpart)%node)=array_g(:,gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/ncomp,ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine assemble_ghosts_nodal_vectorn
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations for the integer scalar.
subroutine assemble_ghosts_nodal_iscalar(array,array_g)
use mpi
use global,only:nnode
implicit none
integer,dimension(nnode),intent(in) :: array
integer,dimension(nnode),intent(out) :: array_g

integer,dimension(maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
integer :: ierr,i_gpart,ncount,ngnode

array_g=array
send_array=0; recv_array=0
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode
  ! store array-to-send in a garray
  send_array(1:ncount,i_gpart)=reshape(array(gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KINT,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KINT,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode
  array_g(gpart(i_gpart)%node)=array_g(gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine assemble_ghosts_nodal_iscalar
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations for the float scalar.
subroutine assemble_ghosts_nodal_fscalar(array,array_g)
use mpi
use global,only:nnode
implicit none
real(kind=kreal),dimension(nnode),intent(in) :: array
real(kind=kreal),dimension(nnode),intent(out) :: array_g

real(kind=kreal),dimension(maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
integer :: ierr,i_gpart,ncount,ngnode

array_g=array
send_array=0; recv_array=0
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode
  ! store array-to-send in a garray
  send_array(1:ncount,i_gpart)=reshape(array(gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode
  array_g(gpart(i_gpart)%node)=array_g(gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine assemble_ghosts_nodal_fscalar
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine assemble_ghosts_gdof(nndof,array,array_g)
use mpi
use global,only:nnode
implicit none
integer,intent(in) :: nndof
integer,dimension(nndof,nnode),intent(in) :: array
integer,dimension(nndof,nnode),intent(out) :: array_g

integer,dimension(nndof*maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
integer :: ierr,i_gpart,ncount,ngnode

array_g=array
send_array=0; recv_array=0
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  ! store array-to-send in a garray
  send_array(1:ncount,i_gpart)=reshape(array(:,gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KINT,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KINT,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  array_g(:,gpart(i_gpart)%node)=array_g(:,gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/nndof,ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine assemble_ghosts_gdof
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine assemble_ghosts_gdof_min(nndof,array,array_g)
use mpi
use global,only:nnode
implicit none
integer,intent(in) :: nndof
integer,dimension(nndof,nnode),intent(in) :: array
integer,dimension(nndof,nnode),intent(out) :: array_g

integer,dimension(nndof*maxngnode,ngpart) :: send_array,recv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
integer :: ierr,i_gpart,ncount,ngnode

array_g=array
send_array=0; recv_array=0
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  ! store array-to-send in a garray
  send_array(1:ncount,i_gpart)=reshape(array(:,gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KINT,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KINT,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  array_g(:,gpart(i_gpart)%node)=array_g(:,gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/nndof,ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine assemble_ghosts_gdof_min
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at gdof locations.
subroutine undo_unmatching_displacementBC(bcnodalu)
use mpi
use global,only:gdof,nnode,nndof,nndofu
implicit none
real(kind=kreal),dimension(nndofu,nnode),intent(inout) :: bcnodalu

real(kind=kreal),dimension(nndofu*maxngnode,ngpart) :: send_array,recv_array
real(kind=kreal),dimension(nndofu,maxngnode) :: garray
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
real(kind=kreal),parameter :: zero=0.0_kreal
integer :: ierr,i_dof,i_node,i_gpart,ignode,ncount,ngnode

send_array=zero; recv_array=zero
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndofu
  ! store array-to-send
  send_array(1:ncount,i_gpart)=reshape(bcnodalu(:,gpart(i_gpart)%node),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndofu
  garray(:,1:ngnode)=reshape(recv_array(1:ncount,i_gpart),(/nndofu,ngnode/))
  do i_node=1,ngnode
    ignode=gpart(i_gpart)%node(i_node)
    do i_dof=1,nndofu
      if(bcnodalu(i_dof,ignode).ne.garray(i_dof,i_node))then
        ! undo
        bcnodalu(i_dof,ignode)=zero
        gdof(i_dof,gpart(i_gpart)%node(i_node))=1
      endif
    enddo
  enddo
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()
return
end subroutine undo_unmatching_displacementBC
!===============================================================================

! this subroutine counts the active ghost partitions for each node on the
! interfaces.
! logical flag representing whether the nodes in the interfaces are intact or
! void has to be communicated across the processors
subroutine count_active_nghosts(ngpart_node)
use mpi
use global,only:nnode
implicit none
! number of active ghost partitions for a node
integer,dimension(nnode),intent(out) :: ngpart_node
! only the interfacial nodes can be saved for the storage (TODO)

logical,dimension(maxngnode,ngpart) :: lsend_array,lrecv_array
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
integer :: i,ierr,i_gpart,ignode,ngnode

ngpart_node=0
lsend_array=.true.; lrecv_array=.true.
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode

  ! store array-to-send
  lsend_array(1:ngnode,i_gpart)=gpart(i_gpart)%isnode(1:ngnode)

  ! send
  call MPI_ISSEND(lsend_array(:,i_gpart),ngnode,MPI_LOGICAL,gpart(i_gpart)%rank,&
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(lrecv_array(:,i_gpart),ngnode,MPI_LOGICAL,gpart(i_gpart)%rank,&
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! count active partitons along the interfaces
do i_gpart=1,ngpart
  do i=1,gpart(i_gpart)%nnode
    ignode=gpart(i_gpart)%node(i)
    if(lrecv_array(i,i_gpart))ngpart_node(ignode)=ngpart_node(ignode)+1
  enddo
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()

return
end subroutine count_active_nghosts
!===============================================================================

! this subroutine distributes the excavation loads discarded by a processors due
! to the special geoemtry partition. it will not distribute if the load is used
! within the partition
subroutine distribute2ghosts(gdof,nndof,neq,ngpart_node, &
array,array_g)
use mpi
use global,only:nnode
implicit none
integer,intent(in) :: nndof,neq
integer,dimension(nndof,nnode),intent(in) :: gdof ! global degree of freedom
! number of active ghost partitions for a node
integer,dimension(nnode),intent(in) :: ngpart_node
real(kind=kreal),dimension(nndof,nnode),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g

real(kind=kreal),dimension(nndof,nnode) :: tarray
!logical,dimension(maxngnode,ngpart) :: lsend_array,lrecv_array
real(kind=kreal),dimension(nndof*maxngnode,ngpart) :: send_array,recv_array
real(kind=kreal),dimension(nndof,maxngnode) :: garray
integer,parameter :: tag=0
integer, dimension(MPI_STATUS_SIZE) :: mpistatus
integer,dimension(ngpart) :: send_req,recv_req
!integer,dimension(nnode) :: ngpart_node ! number of ghost partition for a node
real(kind=kreal),parameter :: zero=0.0_kreal
integer :: i,j,ierr,i_gpart,igdof,ignode,ncount,ngnode

array_g=zero
tarray=array; send_array=zero; recv_array=zero
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  ! store array-to-send in a garray
  garray=zero
  garray(:,1:ngnode)=array(:,gpart(i_gpart)%node)
  ! make appropriate correction so that only the discarded loads will equally be
  ! distributed
  do j=1,ngnode
    ignode=gpart(i_gpart)%node(j)
    if(ngpart_node(ignode)==0)then
      ! this node is dead
      garray(:,j)=zero !stop 'strange!'
      cycle
    endif
    do i=1,nndof
      if(gdof(i,ignode)/=0)then
        ! do not distribute if the load is already used in the partition
        garray(i,j)=zero !garray(i,j)/real(ngpart_node(ignode)+1,kreal) !
        !tarray(i,ignode)=garray(i,j)
      else
        ! distribute equally to ghost partitions if the load is discarded by the
        ! partition.
        !if(ngpart_node(ignode)==0)stop 'strange!'
        garray(i,j)=garray(i,j)/real(ngpart_node(ignode),kreal)
        tarray(i,ignode)=zero
      endif
    enddo
  enddo
  send_array(1:ncount,i_gpart)=reshape(garray(:,1:ngnode),(/ncount/))

  ! send
  call MPI_ISSEND(send_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,  &
  tag,MPI_COMM_WORLD,send_req(i_gpart),ierr)
  ! receive
  call MPI_IRECV(recv_array(:,i_gpart),ncount,MPI_KREAL,gpart(i_gpart)%rank,   &
  tag,MPI_COMM_WORLD,recv_req(i_gpart),ierr)
enddo

! wait for receive-communications completion (recv)
do i_gpart=1,ngpart
  call MPI_WAIT(recv_req(i_gpart),mpistatus,ierr)
enddo

! adding contributions of all ghost neighbours
do i_gpart=1,ngpart
  ngnode=gpart(i_gpart)%nnode
  ncount=ngnode*nndof
  tarray(:,gpart(i_gpart)%node)=tarray(:,gpart(i_gpart)%node)+ &
  reshape(recv_array(1:ncount,i_gpart),(/nndof,ngnode/))
enddo

! wait for send communications completion (send)
do i_gpart=1,ngpart
  call MPI_WAIT(send_req(i_gpart),mpistatus,ierr)
enddo
call sync_process()

! store nodal values to gdof locations
do j=1,nnode
  do i=1,nndof
    igdof=gdof(i,j)
    array_g(igdof)=tarray(i,j)
  enddo
enddo
array_g(0)=zero
return
end subroutine distribute2ghosts
!===============================================================================

! deallocate ghost variables
subroutine cleanup_ghost()
implicit none
integer :: i
do i=1,ngpart
  deallocate(gpart(i)%node,gpart(i)%gdof)
enddo

if(allocated(gpart))deallocate(gpart)
return
end subroutine cleanup_ghost
!===============================================================================

end module ghost_library_mpi
!===============================================================================

