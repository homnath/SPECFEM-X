! collection of solvers
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module parsolver
use set_precision
use global, only : g_num,nedof,nndof,nnode
use ksp_constants, only : KSP_MAXITER,KSP_RTOL
use math_constants, only : zero,zerotol
!use math_library !, only : dot_product_par,maxscal
use math_library_mpi
use ghost_library_mpi

contains
!-------------------------------------------------------------------------------

! conjuate-gradient solver
subroutine ksp_cg_solver(neq,nelmt,k,u_g,f,       &
gdof_elmt,ksp_iter,errcode,errtag)
!use math_library
implicit none
integer,intent(in) :: neq,nelmt ! nelmt (for intact) may not be same as global nelmt
real(kind=kreal),dimension(nedof,nedof,nelmt),intent(in) :: k ! only for intact elements
real(kind=kreal),dimension(0:neq),intent(inout) :: u_g
real(kind=kreal),dimension(0:neq),intent(in) :: f
!integer,dimension(nndof,nnode),intent(in) :: gdof
integer,dimension(nedof,nelmt),intent(in) :: gdof_elmt ! only for intact elements
integer,intent(out) :: ksp_iter
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt
integer,dimension(nedof) :: egdof
real(kind=kreal) :: alpha,beta,rz
real(kind=kreal),dimension(0:neq) :: kp,p,p_g,r,r_g
real(kind=kreal),dimension(nedof,nedof) :: km
real(kind=kreal) :: maxpg,maxug

errtag="ERROR: unknown!"
errcode=-1

!---CG solver

! check if RHS is 0
if(maxval(abs(f)).le.zerotol)then
  u_g=zero
  errcode=0
  return
endif

kp=zero
if(maxval(abs(u_g)).gt.zero)then
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,u_g(egdof))
  enddo
  kp(0)=zero
endif
r=f-kp

call assemble_ghosts(nndof,neq,r,r_g)
p=r
!----cg iteration----
cg: do ksp_iter=1,KSP_MAXITER
  call assemble_ghosts(nndof,neq,p,p_g)
  kp=zero
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,p_g(egdof))
  enddo
  kp(0)=zero

  rz=dot_product_par(r,r_g)
  alpha=rz/dot_product_par(p_g,kp)
  u_g=u_g+alpha*p_g

  maxpg=maxscal(maxval(abs(p_g)))
  maxug=maxscal(maxval(abs(u_g)))
  !if(abs(alpha)*maxscal(maxval(abs(p_g)))/maxscal(maxval(abs(u_g))).le.KSP_RTOL)then
  if(abs(alpha)*maxpg/maxug.le.KSP_RTOL)then
    errcode=0
    return
  endif
  !if(myrank==0)write(*,*)ksp_iter,abs(alpha)*maxpg/maxug
  r=r-alpha*kp
  call assemble_ghosts(nndof,neq,r,r_g)
  beta=dot_product_par(r,r_g)/rz
  p=r+beta*p
  !if(myrank==0)write(*,'(i3,f25.18,f25.18,f25.18)')ksp_iter,alpha,beta,rz
  !if(myrank==0)write(*,*)ksp_iter,alpha,beta,rz
enddo cg
write(errtag,'(a)')'ERROR: CG solver doesn''t converge!'
return
end subroutine ksp_cg_solver
!============================================

! diagonally preconditioned conjuate-gradient solver
subroutine ksp_pcg_solver(neq,nelmt,k,u_g,f,dprecon_g,       &
gdof_elmt,ksp_iter,errcode,errtag)
!use math_library
implicit none
integer,intent(in) :: neq,nelmt ! nelmt (for intact) may not be same as global nelmt
real(kind=kreal),dimension(nedof,nedof,nelmt),intent(in) :: k ! only for intact elements
real(kind=kreal),dimension(0:neq),intent(inout) :: u_g
real(kind=kreal),dimension(0:neq),intent(in) :: f,dprecon_g
!integer,dimension(nndof,nnode),intent(in) :: gdof
integer,dimension(nedof,nelmt),intent(in) :: gdof_elmt ! only for intact elements
integer,intent(out) :: ksp_iter
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt
integer,dimension(nedof) :: egdof
real(kind=kreal) :: alpha,beta,rz
real(kind=kreal),dimension(0:neq) :: kp,p,p_g,r,z,z_g
real(kind=kreal),dimension(nedof,nedof) :: km
real(kind=kreal) :: maxpg,maxug

errtag="ERROR: unknown!"
errcode=-1

!---PCG solver

! check if RHS is 0
if(maxval(abs(f)).le.zerotol)then
  u_g=zero
  errcode=0
  return
endif

kp=zero
if(maxval(abs(u_g)).gt.zero)then
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,u_g(egdof))
  enddo
  kp(0)=zero
endif
r=f-kp
z=dprecon_g*r

call assemble_ghosts(nndof,neq,z,z_g)
p=z
!----pcg iteration----
pcg: do ksp_iter=1,KSP_MAXITER
  call assemble_ghosts(nndof,neq,p,p_g)
  kp=zero
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,p_g(egdof))
  enddo
  kp(0)=zero

  rz=dot_product_par(r,z_g)
  alpha=rz/dot_product_par(p_g,kp)
  u_g=u_g+alpha*p_g

  maxpg=maxscal(maxval(abs(p_g)))
  maxug=maxscal(maxval(abs(u_g)))
  !if(abs(alpha)*maxscal(maxval(abs(p_g)))/maxscal(maxval(abs(u_g))).le.KSP_RTOL)then
  if(abs(alpha)*maxpg/maxug.le.KSP_RTOL)then
    errcode=0
    return
  endif
  !if(myrank==0)write(*,*)ksp_iter,abs(alpha)*maxpg/maxug
  r=r-alpha*kp
  z=dprecon_g*r
  call assemble_ghosts(nndof,neq,z,z_g)
  beta=dot_product_par(r,z_g)/rz
  p=z+beta*p
  !if(myrank==0)write(*,'(i3,f25.18,f25.18,f25.18)')ksp_iter,alpha,beta,rz
  !if(myrank==0)write(*,*)ksp_iter,alpha,beta,rz
enddo pcg
write(errtag,'(a)')'ERROR: PCG solver doesn''t converge!'
return
end subroutine ksp_pcg_solver
!===============================================================================

end module parsolver
!===============================================================================
