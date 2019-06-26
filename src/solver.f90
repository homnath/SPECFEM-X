! collection of solvers
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module solver
use set_precision
use global, only : g_num,nedof, GPU_pointer
use ksp_constants, only : KSP_MAXITER,KSP_RTOL
use math_constants, only : zero,zerotol

contains

! conjuate-gradient solver
subroutine ksp_cg_solver(neq,nelmt,k,u,f,gdof_elmt,ksp_iter,errcode,errtag)
implicit none
integer,intent(in) :: neq,nelmt ! nelmt (for intact) may not be same as global nelmt
real(kind=kreal),dimension(nedof,nedof,nelmt),intent(in) :: k ! only for intact elements
real(kind=kreal),dimension(0:neq),intent(inout) :: u
real(kind=kreal),dimension(0:neq),intent(in) :: f
!integer,dimension(nndof,nnode),intent(in) :: gdof
integer,dimension(nedof,nelmt),intent(in) :: gdof_elmt ! only for intact elements
integer,intent(out) :: ksp_iter
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt
integer,dimension(nedof) :: egdof
real(kind=kreal) :: alpha,beta,rz
real(kind=kreal),dimension(0:neq) :: kp,p,r
real(kind=kreal),dimension(nedof,nedof) :: km

errtag="ERROR: unknown!"
errcode=-1

!---CG solver
ksp_iter=0

! check if RHS is 0
if(maxval(abs(f)).le.zerotol)then
  u=zero
  errcode=0
  return
endif

kp=zero
if(maxval(abs(u)).gt.zero)then
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,u(egdof))
  enddo
  kp(0)=zero
endif
r=f-kp

p=r
!----cg iteration----
cg: do ksp_iter=1,KSP_MAXITER
  kp=zero
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,p(egdof))
  enddo
  kp(0)=zero
  
  rz=dot_product(r,r)
  alpha=rz/dot_product(p,kp)
  u=u+alpha*p

  if(abs(alpha)*maxval(abs(p))/maxval(abs(u)).le.KSP_RTOL)then
    errcode=0
    return
  endif

  r=r-alpha*kp
  beta=dot_product(r,r)/rz
  p=r+beta*p
  !write(*,'(i3,f25.18,f25.18,f25.18)')ksp_iter,alpha,beta,rz

enddo cg
write(errtag,'(a)')'ERROR: CG solver doesn''t converge!'
return
end subroutine ksp_cg_solver
!============================================

! diagonally preconditioned conjuate-gradient solver
subroutine ksp_pcg_solver(neq,nelmt,k,u,f,dprecon,gdof_elmt,ksp_iter,errcode,errtag)
implicit none
integer,intent(in) :: neq,nelmt ! nelmt (for intact) may not be same as global nelmt
real(kind=kreal),dimension(nedof,nedof,nelmt),intent(in) :: k ! only for intact elements
real(kind=kreal),dimension(0:neq),intent(inout) :: u
real(kind=kreal),dimension(0:neq),intent(in) :: f,dprecon
!integer,dimension(nndof,nnode),intent(in) :: gdof
integer,dimension(nedof,nelmt),intent(in) :: gdof_elmt ! only for intact elements
integer,intent(out) :: ksp_iter
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt
integer,dimension(nedof) :: egdof
real(kind=kreal) :: alpha,beta,rz,rznew,pkp
real(kind=kreal),dimension(0:neq) :: kp,p,r,z,p2,kp2
real(kind=kreal),dimension(nedof,nedof) :: km
real :: t1,t2,t3,t4

errtag="ERROR: unknown!"
errcode=-1

!---PCG solver
ksp_iter=0


! check if RHS is 0
if(maxval(abs(f)).le.zerotol)then
  u=zero
  errcode=0
  return
endif

kp=zero
if(maxval(abs(u)).gt.zero)then
  do i_elmt=1,nelmt
    egdof=gdof_elmt(:,i_elmt)
    km=k(:,:,i_elmt)
    kp(egdof)=kp(egdof)+matmul(km,u(egdof))
  enddo
  kp(0)=zero
endif
r=f-kp
z=dprecon*r

p=z

call cpu_time(t1)

call prepare_gpu(GPU_pointer,k,nedof,nelmt,gdof_elmt,neq,f,dprecon,u,r,p,KSP_RTOL)

call cpu_time(t2)

print*, "Elapsed time for GPU init : ", t2 - t1


call cpu_time(t1)

!----pcg iteration----
pcg: do ksp_iter=1,KSP_MAXITER
  kp=zero
 ! do i_elmt=1,nelmt
  !  egdof=gdof_elmt(:,i_elmt) !reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
  ! km=k(:,:,i_elmt)

 !  kp(egdof)=kp(egdof)+matmul(km,p(egdof))
 !   kp(egdof)=kp(egdof)+matmul(k(:,:,i_elmt),p(egdof))
!  enddo
call cpu_time(t3)
  call compute_matvec_prod(GPU_pointer,p,kp) 
call cpu_time(t4)

print*, 'timing loop ;', t4 - t3
  kp(0)=zero

  call gpu_dot_product(GPU_pointer,r,z,neq+1,rz)

  !print *,'GPU',rz
  !rz=0.0 
  !rz=dot_product(r,z)
  !print *,'CPU',rz
  call gpu_dot_product(GPU_pointer,p,kp,neq+1,pkp)

  ! alpha=rz/dot_product(p,kp)
  alpha=rz/pkp

   u = u+alpha*p
  ! Vector operation:  u = u + alpha*p,
  ! where u and p are vectors and alpha is scalar 
  ! call gpu_daxpy_1(GPU_pointer, u, p, alpha, neq+1) 

  if(abs(alpha)*maxval(abs(p))/maxval(abs(u)).le.KSP_RTOL)then
    errcode=0

    call cpu_time(t2)
    print*, "Elapsed time for the big loop : ", t2 - t1

    return
  endif

  r=r-alpha*kp
  z=dprecon*r

  call gpu_dot_product(GPU_pointer,r,z,neq+1,rznew)
  !beta=dot_product(r,z)/rz
  beta=rznew/rz
  p=z+beta*p
  !write(*,'(i3,f25.18,f25.18,f25.18)')ksp_iter,alpha,beta,rz

enddo pcg

write(errtag,'(a)')'ERROR: PCG solver doesn''t converge!'
return
end subroutine ksp_pcg_solver
!===============================================================================

end module solver
!===============================================================================
