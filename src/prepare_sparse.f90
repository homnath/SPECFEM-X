module sparse
contains
!-------------------------------------------------------------------------------

subroutine prepare_sparse()
use math_library,only:i_uniinv,i8_uniinv
use math_library_mpi,only:maxscal,minscal
use mpi_library,only:check_allocate,sync_process
use global
implicit none
integer,parameter :: kint8=selected_int_kind(13)

integer :: i,j,i_elmt,ielmt,i_count,istat,ncount
integer :: igdof,jgdof
integer :: nmax,ndof
integer :: egdof(NEDOF),gegdof(NEDOF)
integer(kind=kint8),allocatable :: ind0(:),iorder(:)
integer,allocatable :: row0(:),col0(:),grow0(:),gcol0(:)

character(len=12) :: spm
character(len=250) :: fname

integer :: gmin,gmax

! counting nonzero elements in offdiagonal portion
integer :: ig0,ig1

integer :: neq_actual
integer :: idgdof(NEDOF),idggdof(NEDOF)
integer,allocatable :: gdof_read(:,:),gnum_read(:,:)
integer :: ngz,ng0,ng1,np0,maxrank0,neq_read,nnode_read
logical,allocatable :: iseq(:)
integer,allocatable :: igorder(:)
integer :: ierr
character(len=250) :: myfname=" => prepare_sparse.f90"
character(len=500) :: errsrc

errsrc=trim(myfname)//' => prepare_sparse'

if(myrank==0) then
  write(logunit,*) 'preparing sparse matrix...'
endif

nmax=nelmt*(NEDOF*NEDOF)
allocate(col0(nmax),row0(nmax),gcol0(nmax),grow0(nmax),stat=ierr)
call check_allocate(ierr,errsrc)

allocate(ggdof(NNDOF,nnode),stat=ierr)
call check_allocate(ierr,errsrc)

if(nproc.eq.1)then
  ggdof=gdof
else
  ! read global degrees of freedoms from DATABASE files
  ! inner core
  write(spm,'(i10)')myrank
  fname='tmp/'//trim(file_head)//'_ggdof_proc'//trim(adjustl(spm))
  open(unit=10,file=trim(fname),access='stream',form='unformatted',    &
  status='old',action='read',iostat=istat)
  if (istat /= 0)then
    write(logunit,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    stop
  endif
  read(10)nnode_read
  if(nnode.ne.nnode_read)then
    ! Error
    write(logunit,*)'ERROR: nnode & nnode_read mismatch!',nnode,nnode_read
    stop
  endif
  read(10)ggdof
  close(10,status='delete')

  allocate(gdof_read(nndof,nnode),gnum_read(nenode,nelmt))
  fname='tmp/'//trim(file_head)//'_gdof'//trim(adjustl(spm))
  open(unit=10,file=trim(fname),access='stream',form='unformatted',    &
  status='old',action='read',iostat=istat)
  if (istat /= 0)then
    write(logunit,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    stop
  endif
  read(10)nnode_read
  read(10)neq_read
  read(10)gdof_read
  close(10,status='delete')
  if(any(.not.(gdof==gdof_read)))print*,myrank,'gdof NOT_EQUAL_TO gdof_read!'

  fname='tmp/'//trim(file_head)//'_gnum'//trim(adjustl(spm))
  open(unit=10,file=trim(fname),access='stream',form='unformatted',    &
  status='old',action='read',iostat=istat)
  if (istat /= 0)then
    write(logunit,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    stop
  endif
  read(10)nnode_read
  read(10)gnum_read
  close(10,status='delete')
  if(any(.not.(g_num==gnum_read)))print*,myrank, &
    'gnum NOT_EQUAL_TO gnum_read!',count(g_num-gnum_read.ne.0)
endif

! total degrees of freedoms
ngdof=maxscal(maxval(ggdof))

if(myrank==0)write(logunit,'(a,i0)')'Total global degrees of freedom: ',ngdof
flush(logunit)

! precompute ownership range OR partion layout
ng1=ngdof/nproc
ng0=ceiling(real(ngdof)/real(nproc))

np0=ngdof-nproc*ng1

if(np0.eq.0)then
! all processors have equal gdofs
  ngz=ng0
  ig0=myrank*ng0 ! 0-based index
  ig1=ig0+ng0-1
elseif(np0.gt.0)then
! first np0 processors have ng0 gdofs each and remainging processors have ng1
! gdofs each
  maxrank0=np0-1 ! myrank is 0-based
  if(myrank.le.maxrank0)then
    ngz=ng0
    ig0=myrank*ng0 ! 0-based index
    ig1=ig0+ng0-1
  else !myrank.gt.maxrank0
    ngz=ng1
    ig0=np0*ng0+(myrank-np0)*ng1 ! 0-based index
    ig1=ig0+ng1-1
  endif
else
! Error
  write(logunit,*)'ERROR: illegal value of "np0"!'
  stop
endif

allocate(iseq(0:neq))
iseq=.false.
! stage 0: store all elements
ncount=0

do i_elmt=1,nelmt
  ielmt=i_elmt
  egdof=reshape(gdof(:,g_num(:,ielmt)),(/NEDOF/))
  gegdof=reshape(ggdof(:,g_num(:,ielmt)),(/NEDOF/))
  iseq(egdof)=.true.
  idgdof=egdof; idggdof=gegdof
  where(idgdof.gt.0)idgdof=1
  where(idggdof.gt.0)idggdof=1
  if(any((idgdof-idggdof).ne.0))then
    print*,myrank,' gdof mismatch!',ielmt
    print*,'OK:',idgdof-idggdof
    print*,'gdof:',egdof
    print*,'ggdof:',gegdof
    print*,'x:',minval(g_coord(1,g_num(:,ielmt))),maxval(g_coord(1,g_num(:,ielmt)))
    print*,'y:',minval(g_coord(2,g_num(:,ielmt))),maxval(g_coord(2,g_num(:,ielmt)))
    print*,'z:',minval(g_coord(3,g_num(:,ielmt))),maxval(g_coord(3,g_num(:,ielmt)))
    print*,'matID:',mat_id(ielmt)
    stop
  endif
  do i=1,NEDOF
    do j=1,NEDOF
      igdof=egdof(i)
      jgdof=egdof(j)
      if(igdof.gt.0.and.jgdof.gt.0)then
        ncount=ncount+1
        row0(ncount)=igdof
        col0(ncount)=jgdof
        grow0(ncount)=gegdof(i)
        gcol0(ncount)=gegdof(j)
      endif
    enddo
  enddo
enddo
call sync_process

if(count(.not.iseq).gt.1)then
  write(logunit,*)'ERRORSP: some degrees of freedoms missing!',myrank,count(.not.iseq),maxval(gdof)
endif
deallocate(iseq)
neq_actual=maxval(gdof)
call sync_process
! stage 1: assemble duplicates
! sort global indices
allocate(ind0(ncount),iorder(ncount),stat=ierr)
call check_allocate(ierr,errsrc)
ind0=0
do i=1,ncount
  ind0(i)=int(neq,kint8)*(int(row0(i),kint8)-1_kint8)+int(col0(i),kint8)
  if(ind0(i).lt.0)print*,'IMPOSSIBLE:',myrank,neq,row0(i),col0(i),ind0(i)
enddo

call i8_uniinv(ind0,iorder)
nsparse=maxval(iorder)
if(myrank==0)write(logunit,'(a,i0,a,i0)')' neq_local: ',neq,' nsparse_local: ',nsparse
call sync_process
allocate(krow_sparse(nsparse),kcol_sparse(nsparse),stat=ierr)
call check_allocate(ierr,errsrc)
allocate(kgrow_sparse(nsparse),kgcol_sparse(nsparse),stat=ierr)
call check_allocate(ierr,errsrc)

krow_sparse=-1
kcol_sparse=-1
kgrow_sparse=-1
kgcol_sparse=-1
do i_count=1,ncount!nmax
  krow_sparse(iorder(i_count))=row0(i_count)
  kcol_sparse(iorder(i_count))=col0(i_count)
  kgrow_sparse(iorder(i_count))=grow0(i_count)
  kgcol_sparse(iorder(i_count))=gcol0(i_count)
enddo
if(minval(krow_sparse).lt.1.or.minval(kcol_sparse).lt.1.or.                    &
minval(kgrow_sparse).lt.1.or.minval(kgcol_sparse).lt.1)then
  write(logunit,*)'ERROR: local and global indices are less than 1!',                &
  minval(krow_sparse),minval(kcol_sparse),minval(kgrow_sparse),                &
  minval(kgcol_sparse)
  flush(logunit)
  stop
endif

deallocate(row0,col0,grow0,gcol0,ind0,iorder)

! local DOF to global DOF mapping
allocate(l2gdof(0:neq),stat=ierr)
call check_allocate(ierr,errsrc)
l2gdof=-1
ndof=nnode*NNDOF
l2gdof(reshape(gdof, (/ndof/)))=reshape(ggdof, (/ndof/))
if(myrank==0)then
  do i=1,nsparse
    if(kgrow_sparse(i).ne.l2gdof(krow_sparse(i)).or.kgcol_sparse(i).ne.l2gdof(kcol_sparse(i)))then
      print*,'VERY STRANGE!!!!!'
      stop
    endif
  enddo
endif
l2gdof=l2gdof-1 ! PETSC uses 0 indexing
gmin=minscal(minval(l2gdof(1:)))
gmax=maxscal(maxval(l2gdof(1:)))
if(myrank==0)then
  write(logunit,'(a,i0,1x,i0)')' l2gdof range: ',gmin,gmax
  flush(logunit)
endif
call sync_process
if(minval(l2gdof(1:)).lt.0)then
  write(logunit,*)'ERROR: local-to-global indices are less than 1!'
  flush(logunit)
  stop
endif

allocate(igorder(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call i_uniinv(l2gdof(1:),igorder)

if(myrank==0)then
  write(logunit,'(a)')'complete!'
  flush(logunit)
endif
end subroutine prepare_sparse
!===============================================================================

end module sparse
!===============================================================================
