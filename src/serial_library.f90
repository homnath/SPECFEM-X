! Routines requied to mimic similar MPI routines found in mpi_library.f90 and
! others
! REVISION:
!   HNG, Jul 11,2011; Apr 09,2010
module serial_library
use set_precision
integer :: ngpart
contains

subroutine start_process()
use global,only:ismpi,myrank,nproc
implicit none
ismpi=.false. ! serial
myrank=0
nproc=1
return
end subroutine start_process
!===============================================================================

subroutine close_process()
implicit none
stop
return
end subroutine close_process
!===============================================================================

subroutine sync_process()
implicit none
return
end subroutine sync_process
!===============================================================================

! write error and stop
subroutine control_error(errcode,errtag,stdout,myrank)
implicit none
integer,intent(in) :: errcode
character(*),intent(in) :: errtag
integer,intent(in) :: stdout,myrank

if(errcode.eq.0)return

! print error message and stop execution
if(myrank==0)write(stdout,'(a)')errtag
stop
end subroutine control_error
!===============================================================================

! get processor tag
function proc_tag() result(ptag)
!use global,only:myrank,nproc
implicit none
character(len=20) :: ptag

ptag=''

return
end function
!===============================================================================

subroutine prepare_ghost()
use global,only:nnode,nndof
implicit none
return
end subroutine prepare_ghost
!===============================================================================

subroutine prepare_ghost_gdof()

implicit none

return
end subroutine prepare_ghost_gdof
!===============================================================================

subroutine modify_ghost(isnode)
use global,only:nnode
implicit none
logical,intent(in) :: isnode(nnode)
return
end subroutine modify_ghost
!===============================================================================

subroutine assemble_ghosts(nndof,neq,array,array_g)
implicit none
integer,intent(in) :: nndof,neq
real(kind=kreal),dimension(0:neq),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g
array_g=array
return
end subroutine assemble_ghosts
!===============================================================================

! this subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine assemble_ghosts_nodal(nndof,array,array_g)
use global,only:nnode
implicit none
integer,intent(in) :: nndof
real(kind=kreal),dimension(nndof,nnode),intent(in) :: array
real(kind=kreal),dimension(nndof,nnode),intent(out) :: array_g

array_g=array
return
end subroutine assemble_ghosts_nodal
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations for the n-component vector.
subroutine assemble_ghosts_nodal_vectorn(ncomp,array,array_g)
use global,only:nnode
implicit none
integer,intent(in) :: ncomp
real(kind=kreal),dimension(ncomp,nnode),intent(in) :: array
real(kind=kreal),dimension(ncomp,nnode),intent(out) :: array_g


array_g=array
return
end subroutine assemble_ghosts_nodal_vectorn
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations
subroutine assemble_ghosts_nodal_vector(array,array_g)
use global,only:NDIM,nnode
implicit none
real(kind=kreal),dimension(NDIM,nnode),intent(in) :: array
real(kind=kreal),dimension(NDIM,nnode),intent(out) :: array_g

array_g=array
return
end subroutine assemble_ghosts_nodal_vector
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations
subroutine assemble_ghosts_nodal_iscalar(array,array_g)
use global,only:nnode
implicit none
integer,dimension(nnode),intent(in) :: array
integer,dimension(nnode),intent(out) :: array_g

array_g=array
return
end subroutine assemble_ghosts_nodal_iscalar
!===============================================================================

! This subroutine assembles the contributions of all ghost partitions
! at nodal locations
subroutine assemble_ghosts_nodal_fscalar(array,array_g)
use global,only:nnode
implicit none
real(kind=kreal),dimension(nnode),intent(in) :: array
real(kind=kreal),dimension(nnode),intent(out) :: array_g

array_g=array
return
end subroutine assemble_ghosts_nodal_fscalar
!===============================================================================

! this subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine assemble_ghosts_gdof(nndof,array,array_g)
use global,only:nnode
implicit none
integer,intent(in) :: nndof
integer,dimension(nndof,nnode),intent(in) :: array
integer,dimension(nndof,nnode),intent(out) :: array_g

array_g=array
return
end subroutine assemble_ghosts_gdof
!===============================================================================

! this subroutine assembles the contributions of all ghost partitions
! at gdof locations
subroutine undo_unmatching_displacementBC(bcnodalu)
use global,only:gdof,nnode,nndofu
implicit none
real(kind=kreal),dimension(nndofu,nnode),intent(inout) :: bcnodalu
! for serial this is done in apply_bc()
return
end subroutine undo_unmatching_displacementBC
!===============================================================================

! this subroutine counts the active ghost partitions for each node on the
! interfaces.
! logical flag representing whether the nodes in the interfaces are intact or
! void has to be communicated across the processors
subroutine count_active_nghosts(nndof,ngpart_node)
use global,only:nnode
implicit none
integer,intent(in) :: nndof
! number of active ghost partitions for a node
integer,dimension(nnode),intent(out) :: ngpart_node
! only the interfacial nodes can be saved for the storage (TODO)
ngpart_node=0
return
end subroutine count_active_nghosts
!===============================================================================

! this subroutine distributes the excavation loads discarded by a processors due
! to the special geoemtry partition. it will not distribute if the load is used
! within the partition
subroutine distribute2ghosts(gdof,nndof,neq,ngpart_node, &
array,array_g)
use global,only:nnode
implicit none
integer,intent(in) :: nndof,neq
integer,dimension(nndof,nnode),intent(in) :: gdof ! global degree of freedom
! number of active ghost partitions for a node
integer,dimension(nnode),intent(in) :: ngpart_node
real(kind=kreal),dimension(nndof,nnode),intent(in) :: array
real(kind=kreal),dimension(0:neq),intent(out) :: array_g

integer :: i,j,igdof
real(kind=kreal),parameter :: zero=0.0_kreal
! store nodal values to gdof locations
do j=1,nnode
  do i=1,nndof
    igdof=gdof(i,j)
    array_g(igdof)=array(i,j)
  enddo
enddo
array_g(0)=zero
return
end subroutine distribute2ghosts
!===============================================================================

! deallocate ghost variables
subroutine cleanup_ghost()
implicit none
return
end subroutine cleanup_ghost
!===============================================================================

end module serial_library
!===============================================================================
