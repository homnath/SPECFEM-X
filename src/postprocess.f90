! This module contains post processing library routines
module postprocess
use set_precision
use global,only:myrank
implicit none
character(len=250),private :: myfname=' => postprocess.f90'
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

! This subroutine compute the nodal tensor quantities from the elemental
! counterpart.
subroutine element_to_point_vector(velmt,vnodal)
use math_constants,only:ZERO
use global,only:nelmt,nnode,nst,ngll,g_num
implicit none
real(kind=kreal),intent(in) :: velmt(:,:,:)
real(kind=kreal),intent(out) :: vnodal(:,:)

integer :: i_elmt,num(ngll)

! check consistency of the data
if(ubound(velmt,1).ne.ubound(vnodal,1))then
  write(*,*)'ERROR: number of component differs for element and point data!'
  stop
endif
if(ubound(velmt,2).ne.ngll)then
  write(*,'(a,i0,1x,i0)')'ERROR: inconsistent ngll points!',ngll,ubound(velmt,2)
  stop
endif
if(ubound(velmt,3).ne.nelmt)then
  write(*,'(a,i0,1x,i0)')'ERROR: inconsistent nelmt!',nelmt,ubound(velmt,3)
  stop
endif
if(ubound(vnodal,2).ne.nnode)then
  write(*,'(a,i0,1x,i0)')'ERROR: inconsistent nnode!',nnode,ubound(vnodal,2)
  stop
endif
vnodal=ZERO
! sum elemental contribution on all nodes
do i_elmt=1,nelmt ! all solid elements
  num=g_num(:,i_elmt)
  vnodal(:,num)=vnodal(:,num)+velmt(:,:,i_elmt)
enddo

return

end subroutine element_to_point_vector
!===============================================================================

! This subroutine compute the nodal tensor quantities from the elemental
! counterpart.
subroutine compute_nodal_tensor(velmt,vnodal)
use math_constants,only:ZERO
use global,only:nelmt,nnode,nst,ngll,g_num
implicit none
real(kind=kreal),intent(in) :: velmt(nst,ngll,nelmt)
real(kind=kreal),intent(out) :: vnodal(nst,nnode)

integer :: i_elmt,num(ngll)

vnodal=ZERO
! sum elemental contribution on all nodes
do i_elmt=1,nelmt ! all solid elements
  num=g_num(:,i_elmt)
  vnodal(:,num)=vnodal(:,num)+velmt(:,:,i_elmt)
enddo

return

end subroutine compute_nodal_tensor
!===============================================================================

subroutine overburden_stress(nelmt,g_num,z0,p0,k0,stress_local)
use global,only:agrav,nenode,ngll,nst,g_coord,massdens_elmt
implicit none
integer,intent(in) :: nelmt
integer,intent(in) :: g_num(nenode,nelmt)
real(kind=kreal),intent(in) :: z0,p0,k0
real(kind=kreal),intent(inout) :: stress_local(nst,ngll,nelmt)

real(kind=kreal) :: z,szz
integer :: i,i_elmt,num(nenode)

do i_elmt=1,nelmt
  num=g_num(:,i_elmt)

  do i=1,ngll ! loop over integration points
    z=g_coord(3,num(i))
    szz=p0-massdens_elmt(i,i_elmt)*agrav*abs(z-z0) ! compression
    stress_local(3,i,i_elmt)=szz
    stress_local(1,i,i_elmt)=k0*szz
    stress_local(2,i,i_elmt)=k0*szz
  enddo ! i GLL
enddo ! i_elmt
return
end subroutine overburden_stress
!===============================================================================

! TODO: is it possible to compute stress_local only for intact elements just
! in case?
! it seems that the subarray cannot be a receiving array
! this subroutine computes elastic stress from the known displacement
subroutine elastic_stress(nelmt,neq,gnod,g_num,gdof_elmt,x,stress_local)
use global,only:ndim,nedof,nenode,ngnode,ngll,nst,g_coord,bulkmod_elmt,        &
shearmod_elmt
use weakform,only:compute_bmat_stress
use elastic,only:compute_cmat_elastic
use math_library,only:determinant,invert
use integration,only:dshape_hex8,lagrange_gll,dlagrange_gll,gll_weights
implicit none
integer,intent(in) :: nelmt,neq,gnod(8)
integer,intent(in) :: g_num(nenode,nelmt),gdof_elmt(nedof,nelmt)
real(kind=kreal),intent(in) :: x(0:neq)
real(kind=kreal),intent(inout) :: stress_local(nst,ngll,nelmt)

real(kind=kreal) :: detjac
real(kind=kreal) :: cmat(nst,nst),coord(ngnode,ndim),jac(ndim,ndim),           &
deriv(ndim,nenode),bmat(nst,nedof),eld(nedof),eps(nst),sigma(nst)
integer :: egdof(nedof),num(nenode)
integer :: i,i_elmt

do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))
  egdof=gdof_elmt(:,i_elmt)
  eld=x(egdof)

  do i=1,ngll ! loop over integration points
    call compute_cmat_elastic(bulkmod_elmt(i,i_elmt),shearmod_elmt(i,i_elmt),  &
    cmat)
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)
    call invert(jac)!

    deriv=matmul(jac,dlagrange_gll(:,i,:))
    call compute_bmat_stress(deriv,bmat)
    eps=matmul(bmat,eld)
    sigma=matmul(cmat,eps)
    stress_local(:,i,i_elmt)=sigma
  enddo ! i GLL
enddo ! i_elmt
return
end subroutine elastic_stress
!===============================================================================

subroutine elastic_stress_intact(nelmt_intact,neq,gnod,elmt_intact,g_num,      &
gdof_elmt,dshape_hex8,dlagrange_gll,x,stress_local)
use global,only:ndim,nedof,nelmt,nenode,ngnode,ngll,nst,g_coord,bulkmod_elmt,  &
shearmod_elmt
use weakform,only:compute_bmat_stress
use elastic,only:compute_cmat_elastic
use math_library,only:determinant,invert
implicit none
integer,intent(in) :: nelmt_intact,neq,gnod(8)
integer,intent(in) :: elmt_intact(nelmt_intact),g_num(nenode,nelmt_intact), &
gdof_elmt(nedof,nelmt_intact)
real(kind=kreal),intent(in) :: dshape_hex8(ndim,ngnode,ngll),                  &
dlagrange_gll(ndim,ngll,ngll),x(0:neq)
real(kind=kreal),intent(inout) :: stress_local(nst,ngll,nelmt)

real(kind=kreal) :: detjac
real(kind=kreal) :: cmat(nst,nst),coord(ngnode,ndim),jac(ndim,ndim),           &
deriv(ndim,nenode),bmat(nst,nedof),eld(nedof),eps(nst),sigma(nst)
integer :: egdof(nedof),num(nenode)
integer :: i,i_elmt,ielmt

do i_elmt=1,nelmt_intact
  ielmt=elmt_intact(i_elmt)
  num=g_num(:,i_elmt)
  coord=transpose(g_coord(:,num(gnod)))
  egdof=gdof_elmt(:,i_elmt)
  eld=x(egdof)
  do i=1,ngll ! loop over integration points
    call compute_cmat_elastic(bulkmod_elmt(i,i_elmt),shearmod_elmt(i,i_elmt),  &
    cmat)
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)
    call invert(jac)!

    deriv=matmul(jac,dlagrange_gll(:,i,:))
    call compute_bmat_stress(deriv,bmat)
    eps=matmul(bmat,eld)
    sigma=matmul(cmat,eps)
    stress_local(:,i,ielmt)=stress_local(:,i,ielmt)+sigma
    !if(i_elmt==1.and.i==1)then
    !      print*,'s',sigma
    !      print*,'e',eld
    !      !print*,'ev',evpt
    !     stop
    !    endif
  enddo ! i GLL
enddo ! i_elmt
!print*,size(stress_local)
return
end subroutine elastic_stress_intact
!===============================================================================

! this subroutine computes norm of the displacement error
subroutine error_norm(nelmt,gnod,g_num,gdof_elmt,diff_u,u_okada,errnorm,vol,   &
i_f_coord,fault_coords)
use global,only:ndim,nedof,nenode,ngnode,ngll,nst,g_coord,nnode
use weakform,only:compute_bmat_stress
use elastic,only:compute_cmat_elastic
use math_library,only:determinant,invert
use integration,only:dshape_hex8,gll_weights
implicit none
integer,intent(in) :: nelmt,gnod(8),i_f_coord
integer,intent(in) :: g_num(nenode,nelmt),gdof_elmt(nedof,nelmt)
real(kind=kreal),intent(in) :: diff_u(ndim,nnode),u_okada(ndim,nnode),         &
fault_coords(i_f_coord)
real(kind=kreal),intent(inout) :: errnorm,vol

real(kind=kreal) :: detjac,zero=0.0_kreal
real(kind=kreal) ::coord(ngnode,ndim),jac(ndim,ndim),du(NDIM,ngll),&
    uo(ndim,ngll)
integer :: egdof(nedof),num(nenode)
integer :: i,i_elmt

errnorm=zero
vol=zero
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)

  coord=transpose(g_coord(:,num(gnod)))
  egdof=gdof_elmt(:,i_elmt)
  du=diff_u(:,num)
  uo=u_okada(:,num)

  do i=1,ngll ! loop over integration points
    if(ANY(fault_coords==num(i)))then
        write(*,*)'singular node:',g_coord(:,num(i))
        cycle
    end if
    jac=matmul(dshape_hex8(:,:,i),coord)
    detjac=determinant(jac)

    errnorm=errnorm+dot_product(du(:,i),du(:,i))*detjac*gll_weights(i)
    vol=vol+dot_product(uo(:,i),uo(:,i))*detjac*gll_weights(i)
  enddo ! i GLL
enddo ! i_elmt
return
end subroutine error_norm
!===============================================================================

! this subroutine computes perturbed density: \nabla.(\rho s)
subroutine compute_save_density_perturbation(nodalu,errcode,errtag)
use math_constants,only:ZERO,MAG_CONS
use dimensionless,only:DIM_DENSITY
use global,only:devel_nondim,g_num,ndim,ngll,nelmt,storederiv,massdens_elmt, &
                file_head,out_path,ptail
implicit none
real(kind=kreal),intent(in) :: nodalu(:,:)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

real(kind=kreal) :: deriv(ndim,ngll),disp(ndim,ngll),rho(ngll),rhos(ndim,ngll)
real(kind=kreal) :: divrhos,divs,gradrho(ndim)
integer :: num(ngll)
integer :: i,i_dim,i_elmt

integer :: ios
character(len=250) :: fname

errcode=-1
errtag=''
! open file to write
fname=trim(out_path)//trim(file_head)//'_pdensity_gll_'//trim(ptail)
open(unit=11,file=trim(fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  disp=nodalu(:,num)
  rho=massdens_elmt(:,i_elmt)
  do i=1,ngll ! loop over integration points
    deriv=storederiv(:,:,i,i_elmt)
    ! \nabla.s
    divs=dot_product(deriv(1,:),disp(1,:)) + &
         dot_product(deriv(2,:),disp(2,:)) + &
         dot_product(deriv(3,:),disp(3,:))
    ! \nabla\rho
    gradrho=matmul(deriv,rho)
    ! \nabla.(\rho s)
    divrhos=rho(i)*divs+dot_product(gradrho,disp(:,i))
    if(devel_nondim)then
      divrhos=DIM_DENSITY*divrhos
    endif
    write(11)divrhos
  enddo ! i GLL
enddo ! i_elmt
close(11)

errcode=0
return
end subroutine compute_save_density_perturbation
!===============================================================================

! This subroutine computes acceleration due to gravity from potential
! nodals: nodal scalar.
! nodalg: nodal gradient vector.
subroutine compute_gradient_of_scalar(nodals,nodalg)
use math_constants,only:ZERO,MAG_CONS
use global,only:g_num,ndim,ngll,nelmt,storederiv,POT_TYPE,PMAGNETIC
implicit none
real(kind=kreal),intent(in) :: nodals(:)
real(kind=kreal),intent(inout) :: nodalg(:,:)

real(kind=kreal) :: deriv(ndim,ngll),ephi(ngll),gradphi(ndim)
integer :: num(ngll)
integer :: i,i_elmt

nodalg=ZERO
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  ephi=nodals(num)

  do i=1,ngll ! loop over integration points
    deriv=storederiv(:,:,i,i_elmt)
    gradphi=matmul(deriv,ephi)
    nodalg(:,num(i))=nodalg(:,num(i))+gradphi
  enddo ! i GLL
enddo ! i_elmt

if(POT_TYPE==PMAGNETIC)then
  nodalg=MAG_CONS*nodalg
endif
return
end subroutine compute_gradient_of_scalar
!===============================================================================

! This subroutine computes magnetic field from magnetic potential and
! magnetization
! B=\mu_0 (H+M)
! But this routine only computed (H+M)
! nodals: nodal scalar.
! nodalB: magnetic field without \mu_0.
subroutine compute_premagnetic_field(nodals,nodalB)
use math_constants,only:ZERO,MAG_CONS
use global,only:g_num,ndim,ngll,nelmt,storederiv, &
                magnetization_elmt
implicit none
real(kind=kreal),intent(in) :: nodals(:)
real(kind=kreal),intent(inout) :: nodalB(:,:)

real(kind=kreal) :: deriv(ndim,ngll),ephi(ngll),gradphi(ndim)
integer :: num(ngll)
integer :: i,i_elmt

nodalB=ZERO
do i_elmt=1,nelmt
  num=g_num(:,i_elmt)
  ephi=nodals(num)

  do i=1,ngll ! loop over integration points
    deriv=storederiv(:,:,i,i_elmt)
    gradphi=matmul(deriv,ephi)
    nodalB(:,num(i))=nodalB(:,num(i))+gradphi+magnetization_elmt(:,i,i_elmt)
  enddo ! i GLL
enddo ! i_elmt
return
end subroutine compute_premagnetic_field
!===============================================================================

! This routine writes scalar data to files in Ensight Gold format
subroutine write_scalar_to_file(nnode,datav,ext,stag,istep)
use global,only:infbc,nst,out_path,file_head,ptail, &
tstep_sformat,nnode_finite,node_finite,nnode_infinite,node_infinite,savedata
use visual
implicit none
integer,intent(in) :: nnode
real(kind=kreal),intent(in) :: datav(:)
character(len=*) :: ext
character(len=*),optional :: stag
integer,intent(in),optional :: istep

integer :: i,n,npart
integer,allocatable :: ipart(:)
! The variable "spart" must be 80 characters long
character(len=80),allocatable :: spart(:)
character(len=20) :: format_str
character(len=20) :: ttag,vtag
character(len=250) :: out_fname,out_fname_inf

! check for time step tag
format_str='('//trim(tstep_sformat)//')'
ttag=''
if(present(istep))then
  write(ttag,format_str)istep
  ttag='_step'//trim(ttag)
endif
! check for varible name tag
vtag=''
if(present(stag))then
  if(len(trim(stag)).ne.0)then
    vtag='_'//trim(stag)
  endif
endif
! check for extension tag
! following statement cause segmentation error!
!ext=trim(adjustl(ext))
if(len(trim(ext)).le.0)then
  write(*,'(a)')'ERROR: ext is empty!'
  stop
endif

! file name
out_fname=trim(out_path)//trim(file_head)//trim(ttag)//trim(vtag)//&
          trim(ptail)//'.'//trim(ext)
if(savedata%infinite)then
  out_fname_inf=trim(out_path)//trim(file_head)//'_inf'//trim(ttag)//trim(vtag)//&
          trim(ptail)//'.'//trim(ext)
endif

n=ubound(datav,1)
if(n.ne.nnode)then
  write(*,'(a)')'ERROR: wrong bounds of "datav"!'
  stop
endif

if(infbc)then
  npart=2
  allocate(ipart(npart),spart(npart))
  ! Note: we will NOT write all parts in a singel file. Therefore, ipart must be
  ! 1 for all.
  ipart=(/ (1,i=1,npart) /)
  spart(1)='finite region'
  spart(2)='infinite region'
else
  npart=1
  allocate(ipart(npart),spart(npart))
  ipart=(/ (i,i=1,npart) /)
  spart(1)='finite region'
endif

! separately plot finite and infinite regions
if(infbc)then
  ! write only the finite region
  call write_ensight_pernodeSCALAS_part1(out_fname,ipart,spart,1, &
  nnode_finite,node_finite,nnode,real(datav))
  ! write only the infinite region
  if(savedata%infinite)then
    call write_ensight_pernodeSCALAS_part1(out_fname_inf,ipart,spart,2, &
    nnode_infinite,node_infinite,nnode,real(datav))
  endif
else
  ! write only the finite region
  call write_ensight_pernodeSCALAS(out_fname,ipart,spart, &
  nnode,real(datav))
endif
deallocate(ipart,spart)
end subroutine write_scalar_to_file
!===============================================================================

! This routine writes scalar data to files in Ensight Gold format for the free
! surface file.
subroutine write_scalar_to_file_freesurf(nnode,datav,ext,stag,istep,plane)
use global,only:infbc,nst,out_path,file_head,ptail, &
tstep_sformat,nnode_finite,node_finite,nnode_infinite,node_infinite,savedata
use visual
implicit none
integer,intent(in) :: nnode
real(kind=kreal),intent(in) :: datav(:)
character(len=*) :: ext
character(len=*),optional :: stag
integer,intent(in),optional :: istep
logical,intent(in),optional :: plane

integer :: i,n,npart
integer,allocatable :: ipart(:)
logical :: isplane
! The variable "spart" must be 80 characters long
character(len=80),allocatable :: spart(:)
character(len=20) :: format_str
character(len=20) :: ttag,vtag
character(len=250) :: out_fname

! Determine if it is for a plane
isplane=.false.
if(present(plane))then
  if(plane)isplane=.true.
endif
! check for time step tag
format_str='('//trim(tstep_sformat)//')'
ttag=''
if(present(istep))then
  write(ttag,format_str)istep
  ttag='_step'//trim(ttag)
endif
! check for varible name tag
vtag=''
if(present(stag))then
  if(len(trim(stag)).ne.0)then
    vtag='_'//trim(stag)
  endif
endif
! check for extension tag
! following statement cause segmentation error!
!ext=trim(adjustl(ext))
if(len(trim(ext)).le.0)then
  write(*,'(a)')'ERROR: ext is empty!'
  stop
endif

! file name
if(isplane)then
  out_fname=trim(out_path)//trim(file_head)//'_free_surface_plane'//trim(ttag)//trim(vtag)//&
        trim(ptail)//'.'//trim(ext)
else
  out_fname=trim(out_path)//trim(file_head)//'_free_surface'//trim(ttag)//trim(vtag)//&
        trim(ptail)//'.'//trim(ext)
endif

n=ubound(datav,1)
if(n.ne.nnode)then
  write(*,'(a)')'ERROR: wrong bounds of "datav"!'
  stop
endif

npart=1
allocate(ipart(npart),spart(npart))
ipart=(/ (i,i=1,npart) /)
spart(1)='free surface'

! Write only on the free surface
call write_ensight_pernodeSCALAS(out_fname,ipart,spart, &
nnode,real(datav))
deallocate(ipart,spart)
end subroutine write_scalar_to_file_freesurf
!===============================================================================

! This routine writes vector data to files in Ensight Gold format.
subroutine write_vector_to_file(nnode,datav,ext,stag,istep)
use global,only:infbc,nst,out_path,file_head,benchmark_okada,ptail, &
tstep_sformat,nnode_finite,node_finite,nnode_infinite,node_infinite,savedata
use visual
implicit none
integer,intent(in) :: nnode
real(kind=kreal),intent(in) :: datav(:,:)
character(len=*) :: ext
character(len=*),optional :: stag
integer,intent(in),optional :: istep

integer :: i,n,ncomp,npart
integer,allocatable :: ipart(:)
character(len=80),allocatable :: spart(:) ! this must be 80 characters long
character(len=20) :: format_str
character(len=20) :: ttag,vtag
character(len=250) :: out_fname,out_fname_inf

! check for time step tag
format_str='('//trim(tstep_sformat)//')'
ttag=''
if(present(istep))then
  write(ttag,format_str)istep
  ttag='_step'//trim(ttag)
endif
! check for varible name tag
vtag=''
if(present(stag))then
  if(len(trim(stag)).ne.0)then
    vtag='_'//trim(stag)
  endif
endif
! check for extension tag
! following statement cause segmentation error!
!ext=trim(adjustl(ext))
if(len(trim(ext)).le.0)then
  write(*,'(a)')'ERROR: ext is empty!'
  stop
!else
!  if(ext(1:1).ne.'.')ext='.'//trim(ext)
endif
! file name
out_fname=trim(out_path)//trim(file_head)//trim(ttag)//trim(vtag)//&
          trim(ptail)//'.'//trim(ext)
if(savedata%infinite)then
  out_fname_inf=trim(out_path)//trim(file_head)//'_inf'//trim(ttag)//trim(vtag)//&
          trim(ptail)//'.'//trim(ext)
endif

! determine data type
ncomp=ubound(datav,1)
if(ncomp/=1 .and. ncomp/=3 .and. ncomp/=6)then
  write(*,'(a)')'ERROR: invalid number of components for "datav"!'
  print*,ncomp
  stop
endif

n=ubound(datav,2)
if(n.ne.nnode)then
  write(*,'(a)')'ERROR: wrong bounds of "datav"!'
  stop
endif

if(infbc)then
  npart=2
  allocate(ipart(npart),spart(npart))
  ! Note: we will NOT write all parts in a single file. Therefore, ipart must be
  ! 1 for all.
  ipart=(/ (1,i=1,npart) /)
  spart(1)='finite region'
  spart(2)='infinite region'
else
  npart=1
  allocate(ipart(npart),spart(npart))
  ipart=(/ (i,i=1,npart) /)
  spart(1)='finite region'
endif

! separately plot finite and infinite regions
if(infbc)then
  ! write only the finite region
  call write_ensight_pernodeVECAS_part1(out_fname,ipart,spart,1, &
  nnode_finite,node_finite,ncomp,nnode,real(datav))
  ! write only the infinite region
  if(savedata%infinite)then
    call write_ensight_pernodeVECAS_part1(out_fname_inf,ipart,spart,2, &
    nnode_infinite,node_infinite,ncomp,nnode,real(datav))
  endif
else
  ! write only the finite region
  call write_ensight_pernodeVECAS(out_fname,ipart,spart, &
  ncomp,nnode,real(datav))
endif
deallocate(ipart,spart)
end subroutine write_vector_to_file
!===============================================================================

! This routine writes vector data to files in Ensight Gold format
subroutine write_vector_to_file_freesurf(nnode,datav,ext,stag,istep,plane)
use global,only:infbc,nst,out_path,file_head,benchmark_okada,ptail, &
tstep_sformat,nnode_finite,node_finite,nnode_infinite,node_infinite,savedata
use visual
implicit none
integer,intent(in) :: nnode
real(kind=kreal),intent(in) :: datav(:,:)
character(len=*) :: ext
character(len=*),optional :: stag
integer,intent(in),optional :: istep
logical,intent(in),optional :: plane

integer :: i,n,ncomp,npart
integer,allocatable :: ipart(:)
logical :: isplane
character(len=80),allocatable :: spart(:) ! this must be 80 characters long
character(len=20) :: format_str
character(len=20) :: ttag,vtag
character(len=250) :: out_fname

! Determine if it is for a plane
isplane=.false.
if(present(plane))then
  if(plane)isplane=.true.
endif
! check for time step tag
format_str='('//trim(tstep_sformat)//')'
ttag=''
if(present(istep))then
  write(ttag,format_str)istep
  ttag='_step'//trim(ttag)
endif
! check for varible name tag
vtag=''
if(present(stag))then
  if(len(trim(stag)).ne.0)then
    vtag='_'//trim(stag)
  endif
endif
! check for extension tag
! following statement cause segmentation error!
!ext=trim(adjustl(ext))
if(len(trim(ext)).le.0)then
  write(*,'(a)')'ERROR: ext is empty!'
  stop
!else
!  if(ext(1:1).ne.'.')ext='.'//trim(ext)
endif
! file name
if(isplane)then
  out_fname=trim(out_path)//trim(file_head)//'_free_surface_plane'//trim(ttag)//trim(vtag)//&
          trim(ptail)//'.'//trim(ext)
else
  out_fname=trim(out_path)//trim(file_head)//'_free_surface'//trim(ttag)//trim(vtag)//&
          trim(ptail)//'.'//trim(ext)
endif

! determine data type
ncomp=ubound(datav,1)
if(ncomp/=1 .and. ncomp/=3 .and. ncomp/=6)then
  write(*,'(a)')'ERROR: invalid number of components for "datav"!'
  print*,ncomp
  stop
endif

n=ubound(datav,2)
if(n.ne.nnode)then
  write(*,'(a)')'ERROR: wrong bounds of "datav"!'
  stop
endif

npart=1
allocate(ipart(npart),spart(npart))
ipart=(/ (i,i=1,npart) /)
spart(1)='free surface'

! write only the finite region
call write_ensight_pernodeVECAS(out_fname,ipart,spart, &
ncomp,nnode,real(datav))
deallocate(ipart,spart)
end subroutine write_vector_to_file_freesurf
!===============================================================================

!! This routine save data to files Ensight Gold format
!! TODO: make it optional
!subroutine save_data(istep,nnode,nodalu,nodalu_okada,  &
!stress_global,nodalphi,nodalg)
!use global,only:infbc,nst,out_path,file_head,savedata,benchmark_okada,ptail, &
!tstep_sformat,nnode_finite,node_finite
!use dimensionless,only:DIM_L,DIM_GPOT,DIM_G,DIM_MOD
!use math_constants
!use visual
!implicit none
!integer,intent(in) :: istep,nnode
!real(kind=kreal),intent(in) :: nodalu(:,:),nodalu_okada(:,:),          &
!stress_global(:,:),nodalphi(:),nodalg(:,:)
!
!integer :: i,npart
!integer,allocatable :: ipart(:)
!character(len=80),allocatable :: spart(:) ! this must be 80 characters long
!character(len=20) :: format_str
!character(len=250) :: out_fname
!
!format_str='(a,'//trim(tstep_sformat)//',a)'
!if(infbc)then
!  npart=2
!  allocate(ipart(npart),spart(npart))
!  ipart=(/ (i,i=1,npart) /)
!  spart(1)='finite region'
!  spart(2)='infinite region'
!else
!  npart=1
!  allocate(ipart(npart),spart(npart))
!  ipart=(/ (i,i=1,npart) /)
!  spart(1)='finite region'
!endif
!
!if(infbc)then
!  ! write only the finite region
!  ! write displacement vector
!  if(savedata%disp)then
!    write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',   &
!    istep,trim(ptail)//'.dis'
!    call write_ensight_pernodeVECAS_part1(out_fname,ipart,spart, &
!    nnode_finite,node_finite,3,nnode,real(DIM_L*nodalu))
!    if(benchmark_okada.and.istep==0)then
!      ! write okada displacement vector
!      write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_okada'//&
!      trim(ptail)//'.dis'
!      call write_ensight_pernodeVECAS_part1(out_fname,ipart,spart, &
!      nnode_finite,node_finite,3,nnode,real(DIM_L*nodalu_okada))
!    endif
!  endif
!
!  ! write stress tensor
!  if(savedata%stress)then
!    write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',   &
!    istep,trim(ptail)//'.sig'
!    call write_ensight_pernodeVECAS_part1(out_fname,ipart,spart, &
!    nnode_finite,node_finite,6,nnode,real(DIM_MOD*stress_global))
!  endif
!
!  ! write gravity potential
!  if(savedata%gpot)then
!    write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',   &
!    istep,trim(ptail)//'.gpot'
!    call write_ensight_pernodeSCALAS_part1(out_fname,ipart,spart, &
!    nnode_finite,node_finite,nnode,real(DIM_GPOT*nodalphi))
!  endif
!
!  if(savedata%agrav)then
!    write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',   &
!    istep,trim(ptail)//'.grav'
!    call write_ensight_pernodeVECAS_part1(out_fname,ipart,spart, &
!    nnode_finite,node_finite,3,nnode,real(DIM_G*nodalg))
!  endif
!else
!  ! write displacement vector
!  if(savedata%disp)then
!    write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',   &
!    istep,trim(ptail)//'.dis'
!    call write_ensight_pernodeVECAS(out_fname,ipart,spart, &
!    3,nnode,real(DIM_L*nodalu))
!    if(benchmark_okada.and.istep==0)then
!      ! write okada displacement vector
!      write(out_fname,'(a)')trim(out_path)//trim(file_head)//'_okada_step'//&
!      trim(ptail)//'.dis'
!      call write_ensight_pernodeVECAS(out_fname,ipart,spart, &
!      3,nnode,real(DIM_L*nodalu_okada))
!    endif
!  endif
!
!  ! write stress tensor
!  if(savedata%stress)then
!    write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',   &
!    istep,trim(ptail)//'.sig'
!    call write_ensight_pernodeVECAS(out_fname,ipart,spart, &
!    6,nnode,real(DIM_MOD*stress_global))
!  endif
!
!  ! write gravity potential
!  if(savedata%gpot)then
!    write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',   &
!    istep,trim(ptail)//'.gpot'
!    call write_ensight_pernodeSCALAS(out_fname,ipart,spart, &
!    nnode,real(DIM_GPOT*nodalphi))
!  endif
!
!  if(savedata%agrav)then
!    write(out_fname,fmt=format_str)trim(out_path)//trim(file_head)//'_step',   &
!    istep,trim(ptail)//'.grav'
!    call write_ensight_pernodeVECAS(out_fname,ipart,spart, &
!    3,nnode,real(DIM_G*nodalg))
!  endif
!endif
!end subroutine save_data
!!===============================================================================

end module postprocess
!===============================================================================
