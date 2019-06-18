! routines requied for MPI proecsses
! REVISION:
!   HNG, Jul 11,2011; Apr 09,2010
module mpi_library
use mpi
contains

! start MPI processes
subroutine start_process()
use global,only:ismpi,myrank,nproc,stdout
implicit none
integer :: errcode
ismpi=.true. ! parallel
call MPI_INIT(errcode)
if(errcode /= 0) call mpierror('ERROR: cannot initialize MPI!',errcode,stdout)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,errcode)
if(errcode /= 0) call mpierror('ERROR: cannot find processor ID (rank)!',errcode,stdout)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,errcode)
if(errcode /= 0) call mpierror('ERROR: cannot find number of processors!',errcode,stdout)
return
end subroutine start_process
!=======================================================

! close all MPI processes
subroutine close_process()
implicit none
integer :: errcode
call MPI_FINALIZE(errcode)
stop
return
end subroutine close_process
!=======================================================

! MPI error
subroutine mpierror(errtag,errcode,stdout)
implicit none
character(len=*) :: errtag
integer :: errcode,stdout
write(stdout,'(a,a,i4,a)')errtag,' (MPI ERROR code:',errcode,')!'
stop
end subroutine mpierror
!=======================================================

! syncronize all MPI processes
subroutine sync_process()
implicit none
integer :: errcode

call MPI_BARRIER(MPI_COMM_WORLD,errcode)

end subroutine sync_process
!=======================================================

subroutine check_allocate(ierr,errsrc)
implicit none
integer,intent(in) :: ierr
character(len=500),intent(in) :: errsrc
if(ierr.ne.0)then
    write(*,*)'ERROR: insufficient memory!'//trim(errsrc)
    stop
endif
end subroutine check_allocate
!===========================================================

! write error and stop
subroutine control_error(errcode,errtag,stdout,myrank)
implicit none
integer,intent(in) :: errcode
character(*),intent(in) :: errtag
integer,intent(in) :: stdout,myrank
integer :: ierr

if(errcode.eq.0)return
! print error message and stop execution
if(myrank==0)write(stdout,'(a)')trim(errtag)
flush(stdout)
call close_process
! stop all the MPI processes, and exit
write(stdout,'(a)')'aborting MPI...'
call MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
stop
end subroutine control_error
!=======================================================

! get processor tag
function proc_tag() result(ptag)
use global,only:myrank,nproc
implicit none
character(len=20) :: format_str,ptag

write(format_str,*)ceiling(log10(real(nproc)+1.))
format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'

write(ptag,fmt=format_str)'_proc',myrank

return
end function
!=======================================================

end module mpi_library
