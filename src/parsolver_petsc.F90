!AUTHORS:
!Hom Nath Gharti
!Stefano Zhampini
!REFERENCE:
!PETSC documentation
!-------------------------------------------------------------------------------
module parsolver_petsc
!-------------------------------------------------------------------------------
!                    Include files
!-------------------------------------------------------------------------------
!#include <finclude/petsckspdef.h>
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscvec.h90"
!#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
!#include "petsc/finclude/petscpc.h"
!#include "petsc/finclude/petscviewer.h"
!#include "petsc/finclude/petscviewer.h90"
use ksp_constants                                                                
use global                                                                       
use math_library_mpi,only:maxscal,minscal                                          
use mpi_library,only:check_allocate,sync_process                                 
use ghost_library_mpi,only:ngpart,gpart                                          
                                                                                 
!use petscsys                                                                    
!use petscvec                                                                    
!use petscmat                                                                    
use petscksp  
implicit none

PetscBool      flg,flg_ch,flg_lu,flg_ilu,mat_symmetry
PetscInt       petsc_solver_type
integer,parameter :: SUPERLU=2,MUMPS=3
PetscInt       ival,icntl
PetscReal      val

Vec              xvec,bvec,local_vec!,gxvec
Mat              Amat,Fmat!,AmatT
KSP              ksp
PC               pc
PetscReal        mnorm,atol,dtol,rtol
PetscInt         iter,maxiter
! For communications from local to global   
VecScatter             vscat!,pscat,,vscat_all
! Stores l2g map info 
ISLocalToGlobalMapping l2gmap                    
!PetscBool        flg

PetscInt :: nzeros_max,nzeros_min,nzerosoff_max
PetscInt :: ngdof_part
PetscInt :: ig0,ig1
PetscInt,allocatable :: nnzero_diag(:),nnzero_offdiag(:)
PetscErrorCode   ierr
!integer :: ierr
character(len=250),private :: myfname=" => parsolver_petsc.f90"
character(len=500),private :: errsrc
contains
!-------------------------------------------------------------------------------

subroutine petsc_initialize()
implicit none
errsrc=trim(myfname)//' => petsc_initialize'

! initialize petsc
call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

end subroutine petsc_initialize
!===============================================================================

subroutine petsc_create_vector()
implicit none
IS global_is,local_is

errsrc=trim(myfname)//' => petsc_create_vector'

! create vector objects
call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,ngdof,xvec,ierr)
CHKERRA(ierr)
call VecDuplicate(xvec,bvec,ierr)
CHKERRA(ierr)

! local vector
call VecCreateSeq(PETSC_COMM_SELF,neq,local_vec,ierr)
CHKERRA(ierr)

! objects needed for global vector scattering to local vector
! create local and global IS (index set) objects from the array of local and
! global indices
call ISCreateGeneral(PETSC_COMM_WORLD,neq,l2gdof(1:),PETSC_COPY_VALUES,global_is,ierr)
CHKERRA(ierr)
call ISCreateStride(PETSC_COMM_SELF,neq,0,1,local_is,ierr);
CHKERRA(ierr)
! create VecScatter object which is needed to scatter PETSc parallel vectors
call VecScatterCreate(bvec,global_is,local_vec,local_is,vscat,ierr)
CHKERRA(ierr)
call ISDestroy(global_is,ierr) ! no longer necessary
call ISDestroy(local_is,ierr)  ! no longer necessary

end subroutine petsc_create_vector
!===============================================================================

subroutine petsc_matrix_preallocate_size()
implicit none
Vec         nzeror_gvec,nzeror_dvec,nzeror_ovec,iproc_gvec,               &
            interface_gvec,ninterface_dvec,ninterface_ovec,nself_gvec
PetscInt :: i !,ncol_part,nrow_part
PetscInt,allocatable :: nzeros(:),ig_array(:) 
PetscScalar rval

PetscInt :: igdof,maxrank0,n,ng,ng0,ng1,np0
PetscInt,allocatable :: iproc_array(:)
PetscScalar,pointer :: nzeror_array(:),rproc_array(:)
PetscScalar,pointer :: nzeror_darray(:),nzeror_oarray(:),rnself_array(:)
PetscReal :: fac_ni,max_ni,pmax,pmin,rnid,rnioffd,rnd,rnoffd

PetscInt :: ir,ic,igr,igc,ir0,ic0,igr0,igc0
PetscInt :: nd,noffd
PetscInt :: i_bool,i_ndof

PetscInt :: nibool
PetscInt,allocatable :: ibool_interface(:),isg_interface(:),   &
                        nself_array(:)
PetscInt, allocatable :: ninterface_darray(:),ninterface_oarray(:)
PetscScalar,allocatable :: rg_interface(:),rnself_lgarray(:)
PetscScalar,pointer :: rninterface_darray(:),rninterface_oarray(:)

character(len=10) :: ptail
!character(len=60) :: outf_name
PetscInt :: count_diag,count_nsparse

errsrc=trim(myfname)//' => petsc_matrix_preallocate_size'

write(ptail,'(i4)')myrank                                                        
ptail=adjustl(ptail)

! predetermine or read compressed size of the sparse matrix

! count number of nonzeros per row
allocate(nzeros(neq),stat=ierr)
call check_allocate(ierr,errsrc)

nzeros=0;
do i=1,nsparse
  nzeros(krow_sparse(i))=nzeros(krow_sparse(i))+1
enddo
nzeros_max=maxscal(maxval(nzeros))
nzeros_min=minscal(minval(nzeros))
nzerosoff_max=nzeros_max
if(myrank==0)then
  write(logunit,'(a,i0,a,i0,a,i0,1x,i0)')' ngdof: ',ngdof, &
  ' nzeros_max: ',nzeros_max,' nzeros_min: ',nzeros_min,count(krow_sparse==1)
  flush(logunit)
endif
deallocate(nzeros)

! precompute ownership range OR partion layout
ng1=ngdof/nproc
ng0=ceiling(real(ngdof)/real(nproc))

np0=ngdof-nproc*ng1

if(np0.eq.0)then
! ng0=ng1
! all processors have equal gdofs
  ng=ng0
  ig0=myrank*ng0 ! 0-based index
  ig1=ig0+ng0-1
elseif(np0.gt.0)then
! first np0 processors have ng0 gdofs each and remainging processors have ng1
! gdofs each
  maxrank0=np0-1 ! myrank is 0-based
  if(myrank.le.maxrank0)then
    ng=ng0
    ig0=myrank*ng0 ! 0-based index
    ig1=ig0+ng0-1
  else !myrank.gt.maxrank0
    ng=ng1
    ig0=np0*ng0+(myrank-np0)*ng1 ! 0-based index
    ig1=ig0+ng1-1
  endif
else
! Error
  write(logunit,*)'ERROR: illegal value of "np0"!'
  flush(logunit)
  stop
endif
call VecDuplicate(xvec,nzeror_gvec,ierr)
CHKERRA(ierr)
call VecDuplicate(xvec,nzeror_dvec,ierr)
CHKERRA(ierr)
call VecDuplicate(xvec,nzeror_ovec,ierr)
CHKERRA(ierr)
call VecDuplicate(xvec,iproc_gvec,ierr)
CHKERRA(ierr)
call VecDuplicate(xvec,interface_gvec,ierr)
CHKERRA(ierr)
call VecDuplicate(xvec,nself_gvec,ierr)
CHKERRA(ierr)
call VecDuplicate(xvec,ninterface_dvec,ierr)
CHKERRA(ierr)
call VecDuplicate(xvec,ninterface_ovec,ierr)
CHKERRA(ierr)

! assign owner processor ID to each gdof (or row)
allocate(ig_array(ng),rproc_array(ng),stat=ierr)
call check_allocate(ierr,errsrc)
ig_array=(/ (i,i=ig0,ig1) /)
rproc_array=real(myrank)
call VecSetValues(iproc_gvec,ng,ig_array,rproc_array,INSERT_VALUES,ierr);
CHKERRA(ierr)
deallocate(ig_array,rproc_array)
call VecAssemblyBegin(iproc_gvec,ierr)
CHKERRA(ierr)
call VecAssemblyEnd(iproc_gvec,ierr)
CHKERRA(ierr)
call VecMin(iproc_gvec,PETSC_NULL_INTEGER,pmin,ierr)
call VecMax(iproc_gvec,PETSC_NULL_INTEGER,pmax,ierr)
! copy solution to local array
allocate(iproc_array(neq),rproc_array(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call scatter_globalvec(iproc_gvec,rproc_array)
iproc_array=int(rproc_array)
deallocate(rproc_array)
call VecDestroy(iproc_gvec,ierr)
! assign interface ID to each gdofs
rval=1.0
! all DOFs
do i=1,ngpart                                                                    
   nibool=gpart(i)%nnode                                                         
    allocate(ibool_interface(nibool),stat=ierr)                                  
    call check_allocate(ierr,errsrc)                                             
    ibool_interface=gpart(i)%node                                                
    do i_bool=1,nibool                                                           
      do i_ndof=1,NNDOF                                                          
        igdof=ggdof(i_ndof,ibool_interface(i_bool))-1                            
        if(igdof.ge.0)call VecSetValues(interface_gvec,1,igdof,rval,           & 
        INSERT_VALUES,ierr);                                                     
      enddo                                                                      
    enddo                                                                        
    deallocate(ibool_interface)                                                  
enddo     
call VecAssemblyBegin(interface_gvec,ierr)
CHKERRA(ierr)
call VecAssemblyEnd(interface_gvec,ierr)
CHKERRA(ierr)

! copy solution to local array
allocate(isg_interface(neq),rg_interface(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call scatter_globalvec(interface_gvec,rg_interface)
isg_interface=int(rg_interface)
deallocate(rg_interface) 
call VecDestroy(interface_gvec,ierr)

! estimate correction for the number of nonzero entries in the diagonal and
! nondiagonal portion
! self interface
rval=1.0
do i=1,neq
  if(isg_interface(i).eq.1)then    
    call VecSetValues(nself_gvec,1,l2gdof(i),rval,ADD_VALUES,ierr);
  endif
enddo
call VecAssemblyBegin(nself_gvec,ierr)
CHKERRA(ierr)
call VecAssemblyEnd(nself_gvec,ierr)
CHKERRA(ierr)
call VecGetLocalSize(nself_gvec,n,ierr)

allocate(rnself_lgarray(neq),stat=ierr)
call check_allocate(ierr,errsrc)
call scatter_globalvec(nself_gvec,rnself_lgarray)
call VecGetArrayF90(nself_gvec,rnself_array,ierr)
allocate(nself_array(n),stat=ierr)
call check_allocate(ierr,errsrc)
nself_array=int(rnself_array(1:n))
where(nself_array.gt.0)nself_array=nself_array-1 ! subtract self
call VecRestoreArrayF90(nself_gvec,rnself_array,ierr)
call VecDestroy(nself_gvec,ierr)

!outf_name='isg_interface'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')isg_interface
!close(1)
deallocate(isg_interface)

! factor for maximum number of interfaces for each nondiagonal entry of the
! stiffness matrix
! the factor below is valid ONLY for rectagular partitioning of the global model
max_ni=8.0
fac_ni=0.0

! first element
igr0=kgrow_sparse(1)-1
igc0=kgcol_sparse(1)-1
ir0=krow_sparse(1)
ic0=kcol_sparse(1)
nd=0; noffd=0
rnid=0.; rnioffd=0.
if(iproc_array(ir0).eq.iproc_array(ic0))then
  nd=1;
  if(igr0.ne.igc0 .and. rnself_lgarray(ir0).gt.1.0 .and. &
  rnself_lgarray(ic0).gt.1.0)then
    fac_ni=min(max_ni,min(rnself_lgarray(ir0),rnself_lgarray(ic0)))
    rnid=1.0/fac_ni
  endif
else
  noffd=1
  if(igr0.ne.igc0 .and. rnself_lgarray(ir0).gt.1.0 .and. &
  rnself_lgarray(ic0).gt.1.0)then
    fac_ni=min(max_ni,min(rnself_lgarray(ir0),rnself_lgarray(ic0)))
    rnioffd=1.0/fac_ni
  endif
endif
count_diag=0
count_nsparse=nd+noffd
do i=2,nsparse
  igr=kgrow_sparse(i)-1
  igc=kgcol_sparse(i)-1
  ir=krow_sparse(i)
  ic=kcol_sparse(i)
  if(l2gdof(ir).ne.igr.or.l2gdof(ic).ne.igc)then
    write(logunit,'(a,i0,1x,i0,1x,i0,1x,i0)')'strange error in l2gdof: ',l2gdof(ir),igr,l2gdof(ic),igc
    flush(logunit)
    stop
  endif
  if(igr.ne.igr0)then
    ! new row starts
    ! set values computed so far
    count_nsparse=count_nsparse+nd+noffd
    rnd=real(nd)
    rnoffd=real(noffd)
    call VecSetValues(nzeror_dvec,1,igr0,rnd,ADD_VALUES,ierr)
    CHKERRA(ierr)
    call VecSetValues(nzeror_ovec,1,igr0,rnoffd,ADD_VALUES,ierr)
    CHKERRA(ierr)
    
    call VecSetValues(ninterface_dvec,1,igr0,rnid,ADD_VALUES,ierr)
    CHKERRA(ierr)
    call VecSetValues(ninterface_ovec,1,igr0,rnioffd,ADD_VALUES,ierr)
    CHKERRA(ierr)

    ! reset
    nd=0; noffd=0
    rnid=0.; rnioffd=0.
    igr0=igr !kgrow_sparse(i)-1
    igc0=igc !kgcol_sparse(i)-1
    ir0=ir !krow_sparse(i)
    ic0=ic !kcol_sparse(i)

    if(iproc_array(ir0).eq.iproc_array(ic0))then
      nd=1;
      if(igr0.ne.igc0 .and. rnself_lgarray(ir0).gt.0.0 .and. &
      rnself_lgarray(ic0).gt.0.0)then
        fac_ni=min(max_ni,min(rnself_lgarray(ir0),rnself_lgarray(ic0)))
        rnid=1.0/fac_ni
      endif
    else
      noffd=1
      if(igr0.ne.igc0 .and. rnself_lgarray(ir0).gt.0.0 .and. &
      rnself_lgarray(ic0).gt.0.0)then
        fac_ni=min(max_ni,min(rnself_lgarray(ir0),rnself_lgarray(ic0)))
        rnioffd=1.0/fac_ni
      endif
    endif
  else
    ! count
    if(iproc_array(ir).eq.iproc_array(ic))then
      nd=nd+1;
      if(igr.ne.igc .and. rnself_lgarray(ir).gt.0.0 .and. &
      rnself_lgarray(ic).gt.0.0)then
        fac_ni=min(max_ni,min(rnself_lgarray(ir),rnself_lgarray(ic)))
        rnid=rnid+(1.0/fac_ni)
      endif
    else
      noffd=noffd+1
      if(igr.ne.igc .and. rnself_lgarray(ir).gt.0.0 .and. &
      rnself_lgarray(ic).gt.0.0)then
        fac_ni=min(max_ni,min(rnself_lgarray(ir),rnself_lgarray(ic)))
        rnioffd=rnioffd+(1.0/fac_ni)
      endif
    endif
  endif
  if(i.eq.nsparse)then
    ! for last
    count_nsparse=count_nsparse+nd+noffd
    rnd=real(nd)
    rnoffd=real(noffd)
    call VecSetValues(nzeror_dvec,1,igr0,rnd,ADD_VALUES,ierr)
    CHKERRA(ierr)
    call VecSetValues(nzeror_ovec,1,igr0,rnoffd,ADD_VALUES,ierr)
    CHKERRA(ierr)

    call VecSetValues(ninterface_dvec,1,igr0,rnid,ADD_VALUES,ierr)
    CHKERRA(ierr)
    call VecSetValues(ninterface_ovec,1,igr0,rnioffd,ADD_VALUES,ierr)
    CHKERRA(ierr)
  endif
enddo
deallocate(rnself_lgarray)
deallocate(iproc_array)
!if(DEBUG)then
!  print*,'counted nsparse:',myrank,nsparse,count_nsparse,count_nsparse-nsparse
!  print*,'counted ndiag:',myrank,count_diag
!endif
deallocate(krow_sparse,kcol_sparse)
call sync_process
call VecAssemblyBegin(nzeror_dvec,ierr)
call VecAssemblyEnd(nzeror_dvec,ierr)
call VecAssemblyBegin(nzeror_ovec,ierr)
call VecAssemblyEnd(nzeror_ovec,ierr)

call VecAssemblyBegin(ninterface_dvec,ierr)
call VecAssemblyEnd(ninterface_dvec,ierr)
call VecAssemblyBegin(ninterface_ovec,ierr)
call VecAssemblyEnd(ninterface_ovec,ierr)

call VecGetLocalSize(nzeror_dvec,n,ierr)

call VecGetArrayF90(nzeror_dvec,nzeror_darray,ierr)
allocate(nnzero_diag(n),stat=ierr)
call check_allocate(ierr,errsrc)
nnzero_diag=int(nzeror_darray(1:n))

!outf_name='ninterface_self'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')nself_array
!close(1)
!! apply correction for repeatition due to interfaces
where(nnzero_diag.gt.0)nnzero_diag=nnzero_diag-nself_array
deallocate(nself_array)
!if(myrank==0)print*,n,minval(nzeror_darray),maxval(nzeror_darray),&
!minval(nnzero_diag),maxval(nnzero_diag)
!call sync_process

!outf_name='nzeror_diagonal'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')nnzero_diag
!close(1)
call VecRestoreArrayF90(nzeror_dvec,nzeror_darray,ierr)
call VecDestroy(nzeror_dvec,ierr)

call VecGetArrayF90(nzeror_ovec,nzeror_oarray,ierr)
allocate(nnzero_offdiag(n))
nnzero_offdiag=int(nzeror_oarray(1:n))
!outf_name='nzeror_offdiagonal'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')nnzero_offdiag
!close(1)
call VecRestoreArrayF90(nzeror_ovec,nzeror_oarray,ierr)
call VecDestroy(nzeror_ovec,ierr)

! correction
! I do not know why but there are some DOFs where the correction exceeds by 4 or
! 8 therefore to be safe we need to subtract this from all
call VecGetArrayF90(ninterface_dvec,rninterface_darray,ierr)
!where(rninterface_darray.gt.0.0 .and. rninterface_darray.lt.1.0)rninterface_darray=1.0
allocate(ninterface_darray(n),stat=ierr)
call check_allocate(ierr,errsrc)
ninterface_darray=int(rninterface_darray(1:n))
call VecRestoreArrayF90(ninterface_dvec,rninterface_darray,ierr)
call VecDestroy(ninterface_dvec,ierr)
where(ninterface_darray.gt.0)ninterface_darray=ninterface_darray-4
where(ninterface_darray.lt.0)ninterface_darray=0
!outf_name='ninterface_diagonal'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')ninterface_darray
!close(1)
deallocate(ninterface_darray)

call VecGetArrayF90(ninterface_ovec,rninterface_oarray,ierr)
!where(rninterface_oarray.gt.0.0 .and. rninterface_oarray.lt.1.0)rninterface_oarray=1.0
allocate(ninterface_oarray(n),stat=ierr)
call check_allocate(ierr,errsrc)
ninterface_oarray=int(rninterface_oarray(1:n))
call VecRestoreArrayF90(ninterface_ovec,rninterface_oarray,ierr)
call VecDestroy(ninterface_ovec,ierr)
where(ninterface_oarray.gt.0)ninterface_oarray=ninterface_oarray-8
where(ninterface_oarray.lt.0)ninterface_oarray=0
!outf_name='ninterface_offdiagonal'//trim(ptail)
!open(1,file=outf_name,action='write',status='replace')
!write(1,'(i4)')ninterface_oarray
!close(1)
deallocate(ninterface_oarray)
call sync_process

do i=1,nsparse
  rval=1.
  igdof=kgrow_sparse(i)-1 ! fortran index
  call VecSetValues(nzeror_gvec,1,igdof,rval,ADD_VALUES,ierr);
  CHKERRA(ierr)
enddo
call VecAssemblyBegin(nzeror_gvec,ierr)
CHKERRA(ierr)
call VecAssemblyEnd(nzeror_gvec,ierr)
CHKERRA(ierr)
call VecGetLocalSize(nzeror_gvec,n,ierr)
CHKERRA(ierr)
if(myrank==0)then
  write(logunit,'(a,4(i0,1x))')' size of vector: ',ng,n,minval(kgrow_sparse),ig0
  flush(logunit)
endif
deallocate(kgrow_sparse,kgcol_sparse)
call VecGetArrayF90(nzeror_gvec,nzeror_array,ierr)
CHKERRA(ierr)

call VecRestoreArrayF90(nzeror_gvec,nzeror_array,ierr)
CHKERRA(ierr)
call VecDestroy(nzeror_gvec,ierr)
CHKERRA(ierr)
where(nnzero_diag.lt.0)nnzero_diag=0
where(nnzero_offdiag.lt.0)nnzero_offdiag=0
where(nnzero_diag.gt.ng)nnzero_diag=ng
call sync_process
if(myrank==0)then
  write(logunit,'(a,i0)')'success!',nzeros_max
  flush(logunit)
endif

end subroutine petsc_matrix_preallocate_size
!===============================================================================

subroutine petsc_create_matrix()
implicit none
PetscInt :: istart,iend

errsrc=trim(myfname)//' => petsc_create_matrix'

! create the matrix and preallocate
call MatCreate(PETSC_COMM_WORLD,Amat,ierr)
call MatSetType(Amat,MATMPIAIJ,ierr)
CHKERRA(ierr)
call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,ngdof,ngdof,ierr)
CHKERRA(ierr)

!does not work before preallocation
!call MatGetLocalSize(Amat,nrow_part,ncol_part,ierr)

! preallocation
!call MatMPIAIJSetPreallocation(Amat,nzeros_max,PETSC_NULL_INTEGER,nzeros_max, &
!PETSC_NULL_INTEGER,ierr)
call MatMPIAIJSetPreallocation(Amat,0,nnzero_diag,         &
0,nnzero_offdiag,ierr)
CHKERRA(ierr)
deallocate(nnzero_diag,nnzero_offdiag)
!if(myrank==0)print*,'ngdof:',ngdof
call MatSetFromOptions(Amat,ierr)
CHKERRA(ierr)

call MatGetOwnershipRange(Amat,istart,iend,ierr)
CHKERRA(ierr)
call sync_process

end subroutine petsc_create_matrix
!===============================================================================

subroutine petsc_create_solver()
implicit none

errsrc=trim(myfname)//' => petsc_create_solver'

! Create the linear solver and set various options

! Create linear solver context
call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
CHKERRA(ierr)
call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)

! diagonally scale the matrix
! since the euqutions are nondimensionalized, the scaling is not necessary?
call KSPSetDiagonalScale(ksp,PETSC_TRUE,ierr)
CHKERRA(ierr)

! define solver type
! this will be overwritten by the command line arguments if provided
! default solver type
petsc_solver_type=1

if(petsc_solver_type==0)then
  if(myrank==0)then
    write(logunit,'(a)')'Solver type: GMRES'
    flush(logunit)
  endif
  !call KSPSetType(ksp,KSPMINRES,ierr);
  !call KSPSetType(ksp,KSPBCGSL,ierr);
  call KSPSetType(ksp,KSPGMRES,ierr);
  !call KSPSetType(ksp,KSPFGMRES,ierr);
  !call KSPSetType(ksp,KSPLGMRES,ierr);
  !call KSPSetType(ksp,KSPGCR,ierr);
  !call KSPSetType(ksp,KSPCG,ierr);
  !call KSPSetType(ksp,KSPPREONLY,ierr);
  CHKERRA(ierr)
  call KSPGetPC(ksp,pc,ierr)
  CHKERRA(ierr)
  !call PCSetType(pc,PCGAMG,ierr)
  call PCSetType(pc,PCPBJACOBI,ierr)
  CHKERRA(ierr)
elseif(petsc_solver_type==1)then
  if(myrank==0)then
    write(logunit,'(a)')'Solver type: CG'
    flush(logunit)
  endif
  call KSPSetType(ksp,KSPCG,ierr);
  CHKERRA(ierr)
  call KSPGetPC(ksp,pc,ierr)
  CHKERRA(ierr)
  !call PCSetType(pc,PCNONE,ierr)
  !call PCSetType(pc,PCHYPRE,ierr)
  call PCSetType(pc,PCBJACOBI,ierr)
  !call PCSetType(pc,PCGAMG,ierr)
  CHKERRA(ierr)
  call PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
  CHKERRA(ierr)
elseif(petsc_solver_type.eq.SUPERLU)then
  if(myrank==0)then
    write(logunit,'(a)')'Solver type: SUPERLU'
    flush(logunit)
  endif
  flg_ilu = PETSC_FALSE;
  flg_lu     = PETSC_FALSE;
  ! version < 3.8.0
  ! call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_superlu_lu",flg_lu,flg,ierr);
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
  "-use_superlu_lu",flg_lu,flg,ierr);
  CHKERRA(ierr)
  if(flg_lu .or. flg_ilu)then
    call KSPSetType(ksp,KSPPREONLY,ierr);
    CHKERRA(ierr)
    call KSPGetPC(ksp,pc,ierr);
    CHKERRA(ierr)
    if(flg_lu)then
      call PCSetType(pc,PCLU,ierr);
      CHKERRA(ierr)
    elseif(flg_ilu)then
      call PCSetType(pc,PCILU,ierr);
      CHKERRA(ierr)
    endif
    call PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
    CHKERRA(ierr)
    ! version < 3.9
    !call PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU,ierr);
    call PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU,ierr);
    CHKERRA(ierr)
    ! version < 3.9
    !call PCFactorSetUpMatSolverPackage(pc,ierr); ! call MatGetFactor() to create F
    call PCFactorSetUpMatSolverType(pc,ierr); ! call MatGetFactor() to create F
    CHKERRA(ierr)
 
    call PCFactorGetMatrix(pc,Fmat,ierr);
    CHKERRA(ierr)
    !call MatSuperluSetILUDropTol(Fmat,1.e-8,ierr);
    !CHKERRA(ierr)
  endif
elseif(petsc_solver_type.eq.MUMPS)then
  if(myrank==0)then
    write(logunit,'(a)')'Solver type: MUMPS'
    flush(logunit)
  endif
  flg_lu    = PETSC_FALSE;
  flg_ch = PETSC_FALSE;
  ! version < 3.8.0
  !call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_mumps_ch",flg_ch,flg,ierr);
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
  "-use_mumps_ch",flg_ch,flg,ierr);
  if(flg_lu .or. flg_ch)then
    call KSPSetType(ksp,KSPPREONLY,ierr);
    call KSPGetPC(ksp,pc,ierr);
    if(flg_lu)then
      call PCSetType(pc,PCLU,ierr);
    elseif(flg_ch)then
      call MatSetOption(Amat,MAT_SPD,PETSC_TRUE,ierr); ! set MUMPS id%SYM=1
      call PCSetType(pc,PCCHOLESKY,ierr);
    endif
    call PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
    CHKERRA(ierr)
    ! version < 3.9
    !call PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS,ierr);
    !call PCFactorSetUpMatSolverPackage(pc,ierr); ! call MatGetFactor() to create F
    call PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr);
    call PCFactorSetUpMatSolverType(pc,ierr); ! call MatGetFactor() to create F
  
    call PCFactorGetMatrix(pc,Fmat,ierr);
    icntl = 7; ival = 2;
    call MatMumpsSetIcntl(Fmat,icntl,ival,ierr);
    icntl = 1; val = 0.0;
    call MatMumpsSetCntl(Fmat,icntl,val,ierr);
  endif
endif

call KSPSetTolerances(ksp,KSP_RTOL,KSP_ATOL,KSP_DTOL,KSP_MAXITER,ierr)
CHKERRA(ierr)

!  Set runtime options, e.g.,
!    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
!  These options will override those specified above as long as
!  KSPSetFromOptions() is called _after_ any other customization
!  routines.
call KSPSetFromOptions(ksp,ierr)

end subroutine petsc_create_solver
!===============================================================================

subroutine petsc_set_ksp_operator(reuse_pc)
implicit none
logical,intent(in) :: reuse_pc
errsrc=trim(myfname)//' => petsc_set_ksp_operator'
if(reuse_pc)then
  ! reuse preconditioner
  call KSPSetReusePreconditioner(ksp,PETSC_TRUE,ierr)
else
  ! do not reuse preconditioner
  call KSPSetReusePreconditioner(ksp,PETSC_FALSE,ierr)
endif
CHKERRA(ierr)

! set ksp operators
! version < 3.5
!call KSPSetOperators(ksp,Amat,Amat,SAME_PRECONDITIONER,ierr)
call KSPSetOperators(ksp,Amat,Amat,ierr) !version >= 3.5.0
CHKERRA(ierr)
end subroutine petsc_set_ksp_operator
!===============================================================================

subroutine petsc_set_stiffness_matrix(storekmat)
use math_library_mpi,only:sumscal
use ieee_arithmetic
implicit none

real(kind=kreal),intent(in) :: storekmat(:,:,:)                                  
integer :: i,i_elmt,ielmt,j,n,ndzero                                             
integer :: ggdof_elmt(NEDOF)                                                     
                                                                                 
PetscInt irow,jcol                                                               
Vec   vdiag                                                                      
PetscScalar rval                                                                 
PetscScalar,pointer :: diag_array(:)                                             
                                                                                 
real(kind=8) :: xval

! Set and assemble matrix.
!  - Note that MatSetValues() uses 0-based row and column numbers
!  in Fortran as well as in C (as set here in the array "col").

call MatZeroEntries(Amat,ierr)
CHKERRA(ierr)
call sync_process
rval=1.0

! entirely in solid                                                              
do i_elmt=1,nelmt                                                                
  ielmt=i_elmt                                                                   
  ggdof_elmt=reshape(ggdof(:,g_num(:,ielmt)),(/NEDOF/))                          
  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0                            
  do i=1,NEDOF                                                                   
    do j=1,NEDOF                                                                 
    irow=i; jcol=j                                                               
    if(ggdof_elmt(irow).ge.0.and.ggdof_elmt(jcol).ge.0)then                      
    !.and.storekmat_intact_ic(i,j,i_elmt).ne.0.0_kreal)then                      
      xval=storekmat(i,j,ielmt)                                                  
      if(ieee_is_nan(xval).or. .not.ieee_is_finite(xval))then                    
        write(logunit,*)'ERROR: stiffness matrix has nonfinite value/s!',myrank,ielmt,&
        mat_id(ielmt),xval,minval(abs(storekmat)),maxval(abs(storekmat))         
        flush(logunit)
        stop                                                                     
      endif                                                                     
      call MatSetValues(Amat,1,ggdof_elmt(irow),1,ggdof_elmt(jcol),           &  
      storekmat(i,j,ielmt),ADD_VALUES,ierr)                                      
      CHKERRA(ierr)                                                              
    endif                                                                        
    enddo                                                                        
  enddo                                                                          
enddo    

call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
CHKERRA(ierr)
call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
CHKERRA(ierr)

if(symmetric_solver)then
  call MatSetOption(Amat,MAT_SYMMETRIC,PETSC_TRUE,ierr)                            
  CHKERRA(ierr)  
else
  call MatSetOption(Amat,MAT_SYMMETRIC,PETSC_FALSE,ierr)                            
  CHKERRA(ierr)  
endif

!! check symmetry                                                                
!call MatDuplicate(Amat,MAT_DO_NOT_COPY_VALUES,AmatT,ierr)                        
!CHKERRA(ierr)                                                                   
!call MatTranspose(Amat,MAT_INITIAL_MATRIX,AmatT,ierr)                            
!CHKERRA(ierr)                                                                   
!rval=-1.0                                                                       
!call MatAXPY(AmatT,rval,Amat,SAME_NONZERO_PATTERN,ierr)                          
!call MatNorm(Amat,NORM_FROBENIUS,mnorm,ierr)                                    
!if(myrank==0)print*,'Matrix norm:',mnorm                                        
!call MatNorm(AmatT,NORM_FROBENIUS,mnorm,ierr)                                    
!if(myrank==0)print*,'Symmetry norm:',mnorm                                      
!call MatDestroy(AmatT,ierr)                                                      
!if(myrank==0)print*,'matrix setting & assembly complete11!'                     
!call sync_process                                                               
                                                                                 
!call MatCreateVecs(Amat,vdiag,PETSC_NULL_OBJECT,ierr) <3.8.0 version            
call MatCreateVecs(Amat,vdiag,PETSC_NULL_VEC,ierr)                               
call MatGetDiagonal(Amat,vdiag,ierr)                                             
call VecGetLocalSize(vdiag,n,ierr)                                               
CHKERRA(ierr)                                                                    
call VecGetArrayF90(vdiag,diag_array,ierr)                                       
CHKERRA(ierr)                                                                    
ndzero=count(diag_array==0.)                                                     
if(ndzero.gt.0)then                                                              
  write(logunit,*)'WARNING: NZEROs in diagonal:',myrank,n, &
  count(diag_array==0.),minval(abs(diag_array)),maxval(abs(diag_array))                                
  flush(logunit)
endif                                                                            
call VecRestoreArrayF90(vdiag,diag_array,ierr)                                   
call sync_process                                                                
call VecDestroy(vdiag,ierr)
                                                      
end subroutine petsc_set_stiffness_matrix
!===============================================================================

subroutine petsc_set_vector(rload)
use ieee_arithmetic
implicit none
PetscScalar,intent(in) :: rload(0:)
PetscScalar zero
PetscInt    istart,iend

call VecGetOwnershipRange(bvec,istart,iend,ierr)
CHKERRA(ierr)

zero=0.0
call VecSet(bvec,zero,ierr)
call VecSetValues(bvec,neq,l2gdof(1:),rload(1:),ADD_VALUES,ierr);

! assemble vector
call VecAssemblyBegin(bvec,ierr)
call VecAssemblyEnd(bvec,ierr)

end subroutine petsc_set_vector
!===============================================================================

subroutine petsc_set_initialguess(rload)
implicit none
PetscScalar,intent(in) :: rload(0:)
PetscScalar zero
PetscInt    istart,iend

call VecGetOwnershipRange(bvec,istart,iend,ierr)
CHKERRA(ierr)

zero=0.0
call VecSet(bvec,zero,ierr)
call VecSetValues(bvec,neq,l2gdof(1:),rload(1:),ADD_VALUES,ierr);

! assemble vector
call VecAssemblyBegin(bvec,ierr)
call VecAssemblyEnd(bvec,ierr)
if(myrank==0)then
  write(logunit,'(a)')'initial guess setting & assembly complete!'
  flush(logunit)
endif

end subroutine petsc_set_initialguess
!===============================================================================

subroutine petsc_solve(sdata,cg_iter,ireason)
implicit none
PetscScalar sdata(:)
PetscInt    cg_iter
PetscInt    ireason

!! null space
!call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,nullspace,ierr)
!!call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_OBJECT, nullspace,ierr);
!call MatSetNullSpace(Amat,nullspace,ierr)
!call MatNullSpaceRemove(nullspace,bvec,PETSC_NULL_OBJECT,ierr)
!call MatNullSpaceDestroy(nullspace,ierr)
!TMP 
!TMP !call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL_OBJECT, nullspace,ierr);
!TMP !call KSPSetNullSpace(ksp, nullspace,ierr);
!TMP !call MatNullSpaceDestroy(nullspace,ierr);

! Solve the linear system
call KSPSolve(ksp,bvec,xvec,ierr)

! View solver info; we could instead use the option -ksp_view
!call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)

! Check solution and clean up
call KSPGetConvergedReason(ksp,ireason,ierr)
call KSPGetIterationNumber(ksp,cg_iter,ierr)

! copy solution to local array
call scatter_globalvec(xvec,sdata)

end subroutine petsc_solve
!===============================================================================

subroutine scatter_globalvec(global_vec,larray)
implicit none

Vec         global_vec
PetscScalar larray(:)
PetscInt    n

PetscScalar,pointer :: array_data(:)
!call VecScatterBegin(vscat,global_vec,local_vec,INSERT_VALUES,SCATTER_REVERSE,ierr)
!CHKERRA(ierr)
!call VecScatterEnd(vscat,global_vec,local_vec,INSERT_VALUES,SCATTER_REVERSE,ierr)
!CHKERRA(ierr)
!call VecGetSize(local_vec,n,ierr)
!CHKERRA(ierr)
!call VecGetArray(local_vec,a_v,a_i,ierr)
!CHKERRA(ierr)
!larray(1:n)=a_v(a_i+1:a_i+n) ! TODO: use FORTRAN POINTER
!if(myrank==0)print*,'larray:',minval(larray(1:n)),maxval(larray(1:n))
!call VecRestoreArray(local_vec,a_v,a_i,ierr)
!CHKERRA(ierr)

call VecScatterBegin(vscat,global_vec,local_vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
CHKERRA(ierr)
call VecScatterEnd(vscat,global_vec,local_vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
CHKERRA(ierr)
call VecGetSize(local_vec,n,ierr)
call VecGetArrayF90(local_vec,array_data,ierr)
CHKERRA(ierr)
larray(1:n)=array_data(1:n)
call VecRestoreArrayF90(local_vec,array_data,ierr)
CHKERRA(ierr)

end subroutine scatter_globalvec
!===============================================================================

subroutine petsc_load()
PetscViewer viewer

call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"PetscMatVecStore", &
FILE_MODE_READ,viewer,ierr)
call MatLoad(Amat,viewer,ierr)
call VecLoad(bvec,viewer,ierr)
call PetscViewerDestroy(viewer,ierr)
end subroutine petsc_load
!===============================================================================

subroutine petsc_save()
PetscViewer viewer

call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"PetscMatVecStore", &
FILE_MODE_WRITE,viewer,ierr)
call MatView(Amat,viewer,ierr)
call VecView(bvec,viewer,ierr)
call PetscViewerDestroy(viewer,ierr)
end subroutine petsc_save
!===============================================================================

subroutine petsc_destroy_vector()
implicit none

! Free work space.  All PETSc objects should be destroyed when they
! are no longer needed.
call VecDestroy(local_vec,ierr)
call VecDestroy(xvec,ierr)
call VecDestroy(bvec,ierr)

end subroutine petsc_destroy_vector
!===============================================================================

subroutine petsc_destroy_matrix()
implicit none

! Free work space.  All PETSc objects should be destroyed when they
! are no longer needed.
call MatDestroy(Amat,ierr)

end subroutine petsc_destroy_matrix
!===============================================================================

subroutine petsc_destroy_solver()
implicit none

! Free work space.  All PETSc objects should be destroyed when they
! are no longer needed.
call KSPDestroy(ksp,ierr)

end subroutine petsc_destroy_solver
!===============================================================================

subroutine petsc_finalize()
implicit none

! Free work space.  All PETSc objects should be destroyed when they
! are no longer needed.
call VecScatterDestroy(vscat,ierr)
call PetscFinalize(ierr)

end subroutine petsc_finalize
!===============================================================================

end module parsolver_petsc
!===============================================================================
