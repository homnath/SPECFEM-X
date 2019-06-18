module benchmark
implicit none
contains
!-------------------------------------------------------------------------------

subroutine compute_okada_solution(iserror,nodalu)
use math_constants,only:ZERO,ONE,TWO,HALF
use global,only:NDIM,nelmt,nnode,nndofu,g_coord,g_num,gdof_elmt, &
                ym_blk,nu_blk,savedata,stdout, &
                okada_aw1,okada_aw2,okada_al1,okada_al2,okada_origin, &
                okada_depth,okada_dip,okada_disl1,okada_disl2,okada_disl3
use dimensionless,only:DIM_L
use element,only:hex8_gnode
use okada_solution
use postprocess
#if (USE_MPI)
use math_library_mpi
#else
use math_library_serial
#endif
implicit none
logical,intent(in),optional :: iserror
real(kind=kreal),intent(in),optional :: nodalu(:,:)
!nodalu_okada: okada displacement
real(kind=kreal),allocatable :: nodalu_okada(:,:)
!diff_u: difference between okada, sem displacement
real(kind=kreal),allocatable :: diff_u(:,:)
!coords: holds a single node's coordinates (x,y,z)
!fault coords: coords that lie along fault for Okada, these yeild singular
real(kind=kreal),allocatable :: coords(:),fault_coords(:)
!error norm for benchmark calculation
real(kind=kreal) :: errnorm,errnorm_g
!volume computed by numerical quadrature
real(kind=kreal) :: vol,vol_g
!params for okada
real(kind=kreal) :: alpha,lambda,mu
!'dummy' values to hold derivatives from okada. the result is not used. LL 9/16
real(kind=kreal) :: uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
!okada return value, index for fault coords
integer :: iret,i_f_coord

integer :: i_dof,i_node,istat

if(myrank==0)write(stdout,*)'computing okada displacement'

! Check for optional parameters
if(present(iserror).and.iserror.and. .not.present(nodalu))then
  write(*,*)'ERROR: for benchmark error "nodalu" must be provided!'
  stop
endif

! calculate alpha. assumes uniform material properties
!lambda=ym_blk(1)*nu_blk(1)/((1.+nu_blk(1))*(1.-2.*nu_blk(1)))
!mu=ym_blk(1)/(2.*(1.+nu_blk(1)))
!alpha=(lambda+mu)/(lambda+2.*mu) ! This depends only on Poisson's ratio
alpha=HALF/(ONE-TWO*nu_blk(1))

allocate(coords(ndim),nodalu_okada(nndofu,nnode),     &
fault_coords(nnode),stat=istat)
if(present(iserror).and.iserror)allocate(diff_u(nndofu,nnode),stat=istat)

i_f_coord=1 !index for fault coords

nodalu_okada=ZERO
do i_node=1,nnode
  ! compute nodalu_okada
  coords=g_coord(:,i_node)
  ! correct for origin
  coords(1)=coords(1)-okada_origin(1)
  coords(2)=coords(2)-okada_origin(2)
  ! okada_origin(3) is always ZERO
  if(coords(3).gt.zero)coords(3)=zero
  call dc3d(alpha,coords(1),coords(2),coords(3),okada_depth,okada_dip,      &
            okada_al1,okada_al2,okada_aw1,okada_aw2,okada_disl1,okada_disl2,&
            okada_disl3,nodalu_okada(1,i_node),nodalu_okada(2,i_node),      &
            nodalu_okada(3,i_node),uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
  if(iret==1)then
      ! In this case, we need to set the slip values depending on which
      ! side of the fault the point lies.
      !write(stdout,*)'WARNING: Okada result yeilds singular value!',     &
      !coords(1),coords(2),coords(3)
      fault_coords(i_f_coord)=i_node !save global index of fault node
      i_f_coord = i_f_coord+1
  else if(iret==2)then
      write(stdout,*)'ERROR: z coord must be negative for okada benchmark!'
      write(stdout,*)coords
      stop
  endif
enddo ! i_node

! write Okada displacement
if(savedata%disp)call write_vector_to_file(nnode,DIM_L*nodalu_okada,&
ext='dis',stag='okada') 

if(present(iserror).and.iserror)then
  diff_u=nodalu_okada-nodalu
  call error_norm(nelmt,hex8_gnode,g_num,gdof_elmt,diff_u,nodalu_okada,      &
  errnorm,vol,i_f_coord,fault_coords(1:i_f_coord))
  deallocate(diff_u)
  ! sum across all processors
  errnorm_g=sumscal(errnorm)
  vol_g=sumscal(vol)
  ! relative error norm
  errnorm_g=errnorm_g/vol_g
  if(myrank==0)then
    write(stdout,*)'finished benchmark calculation'
    write(stdout,*)'errnorm_local:',errnorm,'vol_local:',vol
    write(stdout,*)'errnorm:', errnorm_g,'vol:',vol_g
  endif
endif
deallocate(coords,nodalu_okada,fault_coords)

if(myrank==0)then
  write(stdout,*)'finished benchmark calculation'
endif
end subroutine compute_okada_solution
!===============================================================================

end module benchmark
!===============================================================================
