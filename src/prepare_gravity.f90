module gravity
use set_precision
contains
!-------------------------------------------------------------------------------

subroutine prepare_gravity()
use dimensionless,only:NONDIM_ACCEL
use global,only:agrav,logunit,myrank,g_num,ndim,ngll,nelmt,nnode, &
grav0_nodal,dgrav0_elmt,mat_id,mat_domain,devel_nondim
!use global,only:storederiv,dgrav0_elmt
use math_constants,only:ZERO
implicit none

integer :: i,i_elmt

!real(kind=kreal) :: deriv(ndim,ngll),eg0(ngll,ndim),dgrav0(ndim,ndim)
!integer :: num(ngll)

character(len=250) :: myfname=" => prepare_gravity.f90"
character(len=500) :: errsrc

errsrc=trim(myfname)//' => prepare_gravity'

if(myrank==0) then
  write(logunit,'(a)',advance='no') 'preparing gravity...'
endif

allocate(grav0_nodal(ndim,nnode))!,dgrav0_elmt(6,ngll,nelmt))
grav0_nodal=ZERO
! set background gravity
if(devel_nondim)then
  ! gravity acts downward (-Z direction)
  grav0_nodal(3,:)=-NONDIM_ACCEL*agrav
else
  grav0_nodal(3,:)=-agrav
endif

do i_elmt=1,nelmt
  ! coorec formulation is to compute the gravity on the infinite/transition 
  ! elements where there is no displacement DOF
  ! but these gravity values are used to compute the convetion terms which are
  ! ZERO on those elements due to ZERO displacement
  if(mat_domain(mat_id(i_elmt)).eq.1000)grav0_nodal=ZERO
enddo

!dgrav0_elmt=ZERO
!do i_elmt=1,nelmt
!  num=g_num(:,i_elmt)
!  eg0=transpose(grav0_nodal(:,num))
!
!  do i=1,ngll ! loop over integration points
!    deriv=storederiv(:,:,i,i_elmt)
!    dgrav0=matmul(deriv,eg0)
!    ! set to linear elemental indexing
!    dgrav0_elmt(1,i,i_elmt)=dgrav0(1,1)
!    dgrav0_elmt(2,i,i_elmt)=dgrav0(2,2)
!    dgrav0_elmt(3,i,i_elmt)=dgrav0(3,3)
!    dgrav0_elmt(4,i,i_elmt)=dgrav0(1,2)
!    dgrav0_elmt(5,i,i_elmt)=dgrav0(1,3)
!    dgrav0_elmt(6,i,i_elmt)=dgrav0(2,3)
!  enddo ! i GLL
!enddo ! i_elmt

if(myrank==0) then
  write(logunit,'(a)') 'complete!'
endif
return

end subroutine prepare_gravity
!===============================================================================

end module gravity
!===============================================================================
