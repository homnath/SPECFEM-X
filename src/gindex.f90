! DEVELOPER
!   Hom Nath Gharti
! HISTORY
!   Nov 20,2017; Feb 07,2014; Sep 27,2013
module global_dof
contains
!-------------------------------------------------------------------------------

subroutine gindex() !(ipnode)
use set_precision
use global,only:myrank,nenode,nproc,NNDOF,IDOFU,        &
                iseqsource,eqsource_type,file_head,logunit
implicit none
! local parameters
!integer :: ipnode(NENODE)
integer :: nnode
integer :: i_proc,j_proc

integer :: i,i_dof,i_fnode,i_node,j
integer :: ifnode,ig,istat

! local
logical,allocatable :: isnode(:)

! global
integer,allocatable :: ignode(:)
integer,allocatable :: idofid(:,:),jdofid(:,:)
logical,allocatable :: isgnode(:)

integer,allocatable :: gnf(:,:),gnf_read(:,:)

integer :: ig_end,neq_ext

logical,allocatable :: isgnf(:,:),iseq(:)

integer :: igdof,gnf_end,nibool,nnode_ext
integer,allocatable :: gghost(:,:),ighost(:)
character(len=100) :: fhead
character(len=20) :: spm,spn
character(len=250) :: fname,ifname,pfname

integer,allocatable :: tmpvec(:)
integer,allocatable :: tmpmat(:,:)

integer :: pfault_nnode
integer,allocatable :: pfault_inode(:)

integer :: mfault_nnode
integer,allocatable :: mfault_inode(:)

! ghost partitions
integer :: ngpart
type ghost_partition
  integer :: id,nnode
  integer,dimension(:),allocatable :: inode
end type ghost_partition
type(ghost_partition),dimension(:),allocatable :: gpart

if(myrank.ne.0)then
  write(logunit,'(a)')'ERROR: gindex routine must be run on a single processor!'
  flush(logunit)
endif
if(myrank==0)then
  write(logunit,'(a)')'computing global DOFs ...'
  flush(logunit)
endif
!print*,'gindex rank',myrank
! initialize global indices
ig_end=0 ! global node
gnf_end=0    ! global gdof

! loop through the processors
do i_proc=0,nproc-1
  write(spm,'(i10)')i_proc
  ! read partition information
  pfname='tmp/'//trim(file_head)//'_partitioninfo'//trim(adjustl(spm))
  open(unit=22,file=trim(pfname),access='stream',form='unformatted',    &
  status='old',action='read',iostat=istat)
  if (istat /= 0)then
    write(*,'(a)')'ERROR: partition info file "'//trim(pfname)//&
    '" cannot be opened!'
    stop
  endif
  read(22)nnode
  read(22)ngpart
  allocate(gpart(ngpart))
  do i=1,ngpart
    read(22)gpart(i)%id
    read(22)gpart(i)%nnode
    allocate(gpart(i)%inode(gpart(i)%nnode))
    gpart(i)%inode=-9999 ! Initialize
    read(22)gpart(i)%inode
  enddo
  close(22,status='delete')

  !print*,'Processor:',i_proc,' Neighbours:',ngpart

  allocate(isnode(nnode),ignode(nnode),isgnode(nnode))
  isnode=.false.
  ignode=-1;    isgnode=.false.

  !============================================================
  !! read fault information 
  !if(iseqsource.and.eqsource_type.eq.3)then
  !  pfault_nnode=0
  !  pfname=trim(part_path)//'pfault'//trim(adjustl(spm))
  !  open(22,file=pfname,action='read',status='old')
  !  read(22,*)pfault_nnode
  !  allocate(pfault_inode(pfault_nnode))
  !  read(22,*)pfault_inode
  !  close(22)
  !  
  !  mfault_nnode=0
  !  pfname=trim(part_path)//'mfault'//trim(adjustl(spm))
  !  open(22,file=pfname,action='read',status='old')
  !  read(22,*)mfault_nnode
  !  allocate(mfault_inode(mfault_nnode))
  !  read(22,*)mfault_inode
  !  close(22)
  !endif

  ! global nodal indexing

  ! WARNING: is it correct to put these statements here?
  isgnode=.false.

  fhead=trim(file_head)//'_gnode'
  ! copy global indices from preceeding partitions
  do i=1,ngpart
    j_proc=gpart(i)%id
    if(j_proc<i_proc)then
      write(spn,'(i10)')j_proc
      fname='tmp/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      open(unit=10,file=trim(fname),access='stream',form='unformatted',    &
      status='old',action='read',iostat=istat)
      if (istat /= 0)then
        write(*,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
        stop
      endif
      read(10)nibool
      allocate(ighost(nibool))
      read(10)ighost(1:nibool)
      close(10,status='delete')
      isgnode(gpart(i)%inode)=.true.
      ignode(gpart(i)%inode)=ighost
      deallocate(ighost)
    endif
  enddo
  !print*,'Previous largest node ID:',ig_end

  ! indexify global nodes and store in a region array
  ! inner core
  ig=ig_end
  do i_node=1,nnode
    if(.not.isgnode(i_node))then
      ig=ig+1
      isgnode(i_node)=.true.
      ignode(i_node)=ig
    endif
  enddo

  ! save global indices for neighbouring partitions
  write(spm,'(i10)')i_proc
  do i=1,ngpart
    j_proc=gpart(i)%id
    if(j_proc>i_proc)then
      write(spn,'(i10)')j_proc
      fname='tmp/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      open(unit=10,file=trim(fname),access='stream',form='unformatted', &
      status='replace',action='write',iostat=istat)
      if (istat /= 0)then
        write(*,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
        stop
      endif
      write(10)gpart(i)%nnode
      allocate(tmpvec(gpart(i)%nnode))
      tmpvec=ignode(gpart(i)%inode)
      write(10)tmpvec
      deallocate(tmpvec)
      close(10)
    endif
  enddo
  ig_end=maxval(ignode)
  !print*,'Largest node ID:',ig_end
  !write(spm,'(i10)')i_proc
  !fname=trim(part_path)//'gnode_proc'//trim(adjustl(spm))
  !open(unit=10,file=trim(fname),access='stream',form='unformatted',        &
  !status='replace',action='write',iostat=istat)
  !if (istat /= 0)then
  !  write(*,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  !  stop
  !endif
  !write(10)nnode
  !write(10)ignode
  !close(10)

  ! global indexing of degrees of freedoms
  allocate(gnf(NNDOF,nnode),isgnf(NNDOF,nnode))
  gnf=0
  isgnf=.false.

  ! externally defined Dirichlet BCs
  write(spm,'(i10)')i_proc
  ifname='tmp/'//trim(file_head)//'_gdof'//trim(adjustl(spm))
  open(unit=22,file=trim(ifname),access='stream',form='unformatted',    &
  status='old',action='read',iostat=istat)
  if (istat /= 0)then
    write(*,'(a)')'ERROR: file "'//trim(ifname)//'" cannot be opened!'
    stop
  endif
  read(22)nnode_ext
  read(22)neq_ext
  if(nnode.ne.nnode_ext)then
    write(*,*)'ERROR: nnode & nnode_ext mismatch!'
    stop
  endif
  allocate(gnf_read(NNDOF,nnode_ext),iseq(0:neq_ext))
  read(22)gnf_read
  close(22)

  gnf=gnf_read

  !if(i_proc==1)then
  !print*,'before:',NNDOF
  !do j=1,nenode
  ! print*,j,' = ',gnf(:,ipnode(j))
  !enddo
  !endif
 ! copy global indices from preceeding partitions
  fhead=trim(file_head)//'_gdof'
  write(spm,'(i10)')i_proc
  do i=1,ngpart
    j_proc=gpart(i)%id
    if(j_proc<i_proc)then
      write(spn,'(i10)')j_proc
      fname='tmp/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      open(unit=10,file=trim(fname),access='stream',form='unformatted',    &
      status='old',action='read',iostat=istat)
      if (istat /= 0)then
        write(*,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
        stop
      endif
      read(10)nibool
      if(nibool.eq.0)print*,'nibool:',nibool
      allocate(gghost(NNDOF,nibool),idofid(NNDOF,nibool),jdofid(NNDOF,nibool))
      read(10)gghost
      close(10,status='delete')
      isgnf(:,gpart(i)%inode)=.true.
!      if(i_proc==13)then
!      do j=1,nibool
!      print*,i_proc,j_proc,'hilo0:',gnf(:,gpart(i)%inode(j))
!      print*,i_proc,j_proc,'hilo1:',gghost(:,j)
!      enddo
!      endif
      gnf(:,gpart(i)%inode)=gghost

      !! Check for split-fault
      !! Make sure that split varibles are not passed across the split surface
      !! Plus side
      !if(iseqsource.and.eqsource_type.eq.3)then
      !  if(pfault_nnode>0)then
      !    do i_node=1,nibool
      !      if(any(pfault_inode==gpart(i)%inode(i_node)))then
      !        ! Needed only for test below: "MPI interfaces mismatch....
      !        print*,i_proc,j_proc,'hilop0:',(gghost(j,i_node),j=1,NNDOF)
      !        gghost(IDOFUM,i_node)=0
      !        
      !        gnf(IDOFUM,gpart(i)%inode(i_node))=0
      !        isgnf(IDOFUM,gpart(i)%inode(i_node))=.false.
      !        print*,i_proc,j_proc,'hilop1:',(gghost(j,i_node),j=1,NNDOF)
      !      endif
      !    enddo
      !  endif

      !  ! Minus side
      !  if(mfault_nnode>0)then
      !    do i_node=1,nibool
      !      if(any(mfault_inode==gpart(i)%inode(i_node)))then
      !        ! Needed only for test below: "MPI interfaces mismatch....
      !        print*,i_proc,j_proc,'hilom0:',(gghost(j,i_node),j=1,NNDOF)
      !        gghost(IDOFU,i_node)=0
      !        
      !        gnf(IDOFU,gpart(i)%inode(i_node))=0
      !        isgnf(IDOFU,gpart(i)%inode(i_node))=.false.
      !        print*,i_proc,j_proc,'hilom1:',(gghost(j,i_node),j=1,NNDOF)
      !      endif
      !    enddo
      !  endif
      !endif
      
      idofid=gnf(:,gpart(i)%inode)
      jdofid=gghost
      where(idofid>0)idofid=1
      where(jdofid>1)jdofid=1
      if(.not.all(idofid.eq.jdofid))then
        write(*,*)'ERROR: MPI interfaces mismatch between processors:',i_proc,j_proc
        stop
      endif
      deallocate(idofid,jdofid)
      deallocate(gghost)
    endif
  enddo
  !if(i_proc==1)then
  !print*,'NNDOF:',NNDOF
  !print*,'after:'
  !do j=1,nenode
  ! print*,j,' = ',gnf(:,ipnode(j))
  !enddo
  !endif
  
  !if(i_proc.eq.1)print*,'--------------------------------------------'
  !print*,'Previous largest gnf ID:',gnf_end

  ! test only
  iseq=.false.
  do i_node=1,nnode
    do i_dof=1,NNDOF
      iseq(gnf_read(i_dof,i_node))=.true.
    enddo
  enddo
  if(count(.not.iseq).gt.1)then
    write(*,*)'ERRORSP gindex: some degrees of freedoms missing!',i_proc,neq_ext,nnode,&
    count(iseq),count(.not.iseq),neq_ext,minval(gnf_read),maxval(gnf_read)
  endif

  where(gnf_read>0)gnf_read=1
  igdof=gnf_end ! gdof
  do i_node=1,nnode
    do i_dof=1,NNDOF
      if(gnf(i_dof,i_node).gt.0 .and. .not.isgnf(i_dof,i_node))then
        isgnf(i_dof,i_node)=.true.
        igdof=igdof+1
        gnf(i_dof,i_node)=igdof
      endif
    enddo
  enddo
  deallocate(gnf_read,iseq)
  
  ! save global degrees of freedom for neighbouring partitions
  write(spm,'(i10)')i_proc
  do i=1,ngpart
    j_proc=gpart(i)%id
    if(j_proc>i_proc)then
      write(spn,'(i10)')j_proc
      fname='tmp/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      open(unit=10,file=trim(fname),access='stream',form='unformatted',        &
      status='replace',action='write',iostat=istat)
      if (istat /= 0)then
        write(*,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
        stop
      endif
      write(10)gpart(i)%nnode
      allocate(tmpmat(NNDOF,gpart(i)%nnode))
      tmpmat=gnf(:,gpart(i)%inode)
      write(10)tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  gnf_end=maxval(gnf)
  !print*,'Largest gnf ID:',gnf_end
  write(spm,'(i10)')i_proc
  fname='tmp/'//trim(file_head)//'_ggdof_proc'//trim(adjustl(spm))
  open(unit=10,file=trim(fname),access='stream',form='unformatted',        &
  status='replace',action='write',iostat=istat)
  if (istat /= 0)then
    write(*,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    stop
  endif
  write(10)nnode
  write(10)gnf
  close(10)

  deallocate(gnf,isgnf)
  do i=1,ngpart
    deallocate(gpart(i)%inode)
  enddo
  deallocate(gpart)
  deallocate(isnode,ignode,isgnode)
  
enddo !i_proc=0,nproc
if(myrank==0)then
  write(logunit,'(a,i0)')' Exact total number of nodes: ',ig_end
  write(logunit,'(a,i0)')' Exact total number of global DOFs: ',gnf_end
  write(logunit,'(a)')'complete!'
  flush(logunit)
endif

end subroutine gindex
!===============================================================================

end module global_dof
!===============================================================================
