! collection of solvers
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module parsolver
use set_precision
use global, only : g_num,nedof,nndof,nnode,GPU_pointer
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

integer :: i_elmt,ncuda_devices,i_gpart
integer,dimension(nedof) :: egdof
real(kind=kreal) :: alpha,beta,rz, max_p, max_u
real(kind=kreal),dimension(0:neq) :: kp,p,p_g,r,z,z_g
real(kind=kreal),dimension(nedof,nedof) :: km
real(kind=kreal) :: maxpg,maxug
real(kind=kreal),dimension(maxngnode * nndof,ngpart) :: send_buffer, recv_buffer
integer, dimension(maxngnode * nndof,ngpart) :: mpi_gdof
integer, dimension(ngpart) :: mpi_count
logical, parameter :: GPU_MODE = .true.


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

if (GPU_MODE) then

  call initialize_cuda_device(myrank,ncuda_devices)

  do i_gpart = 1,ngpart
   mpi_count(i_gpart) = gpart(i_gpart)%nnode * nndof 
   mpi_gdof(1:gpart(i_gpart)%nnode * nndof ,i_gpart) =  gpart(i_gpart)%gdof
  enddo

  call prepare_gpu(GPU_pointer,k,nedof,nelmt,gdof_elmt,neq,f,dprecon_g,u_g,r,p,z,KSP_RTOL,myrank,NPROC,maxngnode * nndof,ngpart,mpi_count,mpi_gdof)

  pcg_gpu: do ksp_iter=1,KSP_MAXITER

    call gpu_loop1(GPU_pointer, send_buffer)

    call gpu_loop2(GPU_pointer,max_p,max_u,alpha, recv_buffer)

    if (abs(alpha)*(abs(max_p))/abs(max_u) .le. KSP_RTOL) errcode=0

    if (errcode == 0 ) then
      call gpu_loop3(GPU_pointer)
      call gpu_loop4(GPU_pointer,u_g)
      return
    else
      call gpu_loop3(GPU_pointer)
    endif

  enddo pcg_gpu

else ! on CPU

!----pcg iteration----
pcg_cpu: do ksp_iter=1,KSP_MAXITER
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
enddo pcg_cpu

endif

write(errtag,'(a)')'ERROR: PCG solver doesn''t converge!'
return
end subroutine ksp_pcg_solver
!===============================================================================

end module parsolver
!===============================================================================
