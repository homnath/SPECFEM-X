!AUTHORS:
!Hom Nath Gharti
!Stefano Zhampini
!REFERENCE:
!PETSC documentation
! All routines below are empty. these have to be modified for petsc serial solver
! if available. Serial version curerntly use builtin solver.
! -----------------------------------------------------------------------
module solver_petsc
use ksp_constants
use global
implicit none
contains
!=======================================================
subroutine petsc_initialize()
implicit none

end subroutine petsc_initialize
!=======================================================

subroutine petsc_create_vector
implicit none

end subroutine petsc_create_vector
!=======================================================

subroutine petsc_matrix_preallocate_size()

end subroutine petsc_matrix_preallocate_size
!=======================================================

subroutine petsc_create_matrix()

end subroutine petsc_create_matrix
!=======================================================

subroutine petsc_create_solver()

end subroutine petsc_create_solver
!=======================================================

subroutine petsc_set_stiffness_matrix(storekmat)
use ieee_arithmetic
implicit none
real(kind=kreal),intent(in) :: storekmat(:,:,:)

end subroutine petsc_set_stiffness_matrix
!=======================================================

subroutine petsc_set_vector(rload)
!use global,only:l2gdof,nelmt,NEDOF
use ieee_arithmetic
implicit none
!PetscScalar,intent(in) :: rload(0:)
real(kind=kreal),intent(in) :: rload(0:)

end subroutine petsc_set_vector
!=======================================================

subroutine petsc_set_ksp_operator(reuse_pc)
implicit none
logical,intent(in) :: reuse_pc
end subroutine petsc_set_ksp_operator
!=======================================================

subroutine petsc_set_initialguess(rload)
!use global,only:l2gdof,NEDOF
implicit none
!PetscScalar,intent(in) :: rload(0:)
real(kind=kreal),intent(in) :: rload(0:)

end subroutine petsc_set_initialguess
!=======================================================

subroutine petsc_solve(sdata,cg_iter,ireason)
implicit none
!PetscScalar sdata(:)
!PetscInt    cg_iter
real(kind=kreal),intent(in) :: sdata(:)
integer :: cg_iter,ireason

end subroutine petsc_solve
!=======================================================

subroutine petsc_load()
implicit none

end subroutine petsc_load
!=======================================================

subroutine petsc_save()
implicit none

end subroutine petsc_save
!=======================================================

subroutine petsc_destroy_vector()
implicit none

end subroutine petsc_destroy_vector
!===============================================================================

subroutine petsc_destroy_matrix()
implicit none

end subroutine petsc_destroy_matrix
!===============================================================================

subroutine petsc_destroy_solver()
implicit none

end subroutine petsc_destroy_solver
!===============================================================================

subroutine petsc_finalize()
implicit none

end subroutine petsc_finalize
!=======================================================

end module solver_petsc
!=======================================================
