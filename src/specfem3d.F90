! include 'license.txt'
! REVISION:
!   HNG, Aug 25,2011; HNG, Jul 14,2011; HNG, Jul 11,2011; Apr 09,2010
subroutine specfem3d()
! import necessary libraries
use dimensionless
use global
use string_library, only : parse_file
use math_constants
use conversion_constants
use gll_library
use shape_library
use infinite_element
use math_library
use element,only:hex8_gnode
use dof
use fault
use weakform
use gravity
use preprocess
use elastic
use traction
use mtraction
use viscoelastic
use global_dof
use plastic_library
use map_location
use earthquake
#if (USE_MPI)
use mpi_library
use ghost_library_mpi
use math_library_mpi
use sparse
use parsolver
use parsolver_petsc
#else
use serial_library
use math_library_serial
use sparse_serial
use solver
use solver_petsc
#endif
use bc
use free_surface
use benchmark
use visual
use postprocess
implicit none

character(len=250) :: myfname=' => specfem3d.f90'
character(len=500) :: errsrc

!i,j are dummy vars for interation
integer :: i,j
!istat: status indicator for allocation (can be used in other contexts)
integer :: istat

!do-loop indices
integer :: i_elmt,i_nliter,i_node,i_tstep,i_comp
integer :: ielmt,imat,idof,iedof!element ID for gdof, node, etc.

real(kind=kreal),dimension(nst),parameter :: unit_voigt=(/one,one,one,ZERO,    &
ZERO,ZERO /)!Voigt representation for vector
real(kind=kreal) :: jacw !determinant of Jacobian*gll_weight
real(kind=kreal) :: dt ! time step

real(kind=kreal) :: sfac ! slip factor

! KSP convergence reason
integer :: ksp_convreason
!ksp_iter: conjugate gradient iteration,
!ksp_tot: total cg interation, nl_iter: nonlinear iteration, nl_tot: total nl
!iteration
integer :: ksp_iter,ksp_tot,nl_iter,nl_tot
logical :: nl_isconv ! logical variable to check convergence of
! nonlinear (NL) iterations

real(kind=kreal) :: G,K
real(kind=kreal) :: cmat(nst,nst),estrain(nst),dev_strain(nst),  &
esigma(nst),sigma(nst),vsigma(nst)
!cmat: elastic matric (Cijkl) in Voigt notation
!estrain: elastic strain
!esigma: elastic stress
!sigma: stress (used in different context than esigma?)
!vsigma: viscose stress

! dynamic arrays
!num: g_num for particular element.
!node_valency: number of elements that share each node.
integer,allocatable::num(:),node_valency(:)

! factored parameters. Only for elasto-plastic implementation
real(kind=kreal),allocatable :: cohf(:),nuf(:),phif(:),psif(:),ymf(:)

integer :: nzero_dprecon
real(kind=kreal),allocatable :: dprecon(:),ndscale(:)

real(kind=kreal),allocatable :: bmat(:,:),coord(:,:),deriv(:,:),    &
jac(:,:)
!km: stiffness matrix for each element
!storekmat: stiffness matrix for all elements
real(kind=kreal),allocatable :: kmat(:,:),storekmat(:,:,:)

!uerr: used to check convergence
!umax: max of displacement magnitude, uxmax: max of displacement components
real(kind=kreal) :: uerr,maxu,maxdu
!du: incremental soln
!u: solution (summed over du)
real(kind=kreal),allocatable :: du(:),u(:)

real(kind=kreal),allocatable :: eld(:),eload(:),bload(:),   &
vload(:),rhoload(:),resload(:)
!eld: elastic displacement on all nodes of the element
!eload: elastic load on element
!extload: external load
!ubcload: load contributed by displacement BC
!load: like resload. Not currently used
real(kind=kreal),allocatable :: slipload(:),extload(:),bodyload(:), &
viscoload(:),ubcload(:),load(:)

!strain_elmt: strain for all elements
!stress_elmt: stress for each element
!stress_nodal: nodal stress for all elements in processor
real(kind=kreal),allocatable :: strain_elmt(:,:,:),strain_nodal(:,:),      &
stress_elmt(:,:,:),stress_nodal(:,:)
!bcnodalv: prescribed BC nodal variables
!nodalu: nodal displacement for all nodes (not just BC)
real(kind=kreal),allocatable :: bcnodalv(:,:),nodalu(:,:)
real(kind=kreal),allocatable :: nodalphi(:),nodalg(:,:)
! magnetization
real(kind=kreal),allocatable :: nodalB(:,:)
!,psigma(:,:),psigma0(:,:),taumax(:),nsigma(:)
!bodyload: load computed on all nodes of the element
!viscoload: load contributed (on all nodes) by visco elements
!bmat: strain displacement matrix (per node). bmat*displacement vector = strain
!bload: load computed on each dof on a particular element
!vload: like bload, for viscose load
!resload: residual load
!coord: coordinates of geometrical nodes
!deriv: derivative of interpolation functions
!dprecon: diagonal preconditioner. Not used for petsc solver
!ndscale: nondimensionlize scale
!ubcload: load contributed by displacement BC
!jac: Jacobian
integer,allocatable :: egdof(:),egdofu(:)
! placeholder array. holds values of gdof_elmt for a given element.

logical :: isgravity,ispseudoeq ! gravity load and pseudostatic load

! viscoelastic parameters
integer :: mdomain
integer :: ielmt_elas,ielmt_viscoelas,imatve,iviscoelas,nelmt_elas,            &
nelmt_viscoelas
integer :: tot_nelmt_elas,max_nelmt_elas,min_nelmt_elas
integer :: tot_nelmt_viscoelas,max_nelmt_viscoelas,min_nelmt_viscoelas
integer :: nmatblk_elas
real(kind=kreal) :: min_relaxtime,max_relaxtime
! time at current time step
real(kind=kreal) :: t
!factor for time unit conversion
real(kind=kreal) :: tunitfac

integer :: i_maxwell
integer,allocatable :: eid_elas(:),eid_viscoelas(:)
real(kind=kreal),allocatable :: relaxtime(:,:),muratio(:),tratio(:)
! q0: initial state variable for viscoelastic rheology
real(kind=kreal),allocatable :: visco_q0(:,:,:,:),q0(:,:),elas_e0(:,:,:)
real(kind=kreal) :: e0(nst) !e0: initial strain

real(kind=kreal) :: vesigma(nst)
real(kind=kreal) :: trace_vsigma0,trace_strain
real(kind=kreal) :: esigma0_dev(nst),esigma_dev(nst)
real(kind=kreal) :: maxres

integer :: geq,inum,nequ
logical,allocatable :: iseq(:)
integer,allocatable :: gdofu(:) 
integer :: tot_neq,max_neq,min_neq
! number of active ghost partitions for a node
integer,allocatable :: ngpart_node(:)
character(len=250) :: errtag ! error message
integer :: errcode
errtag=""; errcode=-1

! Relaxation time
! Convert relaxation time unit to the time step time unit for consistency
! TODO: check for nondimensionalizing the time step and relaxation time
! Note: for viscosity it doesn't matter because the only relevant parameter is 
! 'tratio' which is dimensionless in itself.
allocate(relaxtime(nmaxwell,nmatblk_viscoelas),muratio(nmaxwell),tratio(nmaxwell))
relaxtime=inftol
tunitfac=ONE
if(index(tunit,'sec').gt.0 .or. index(tunit,'second').gt.0)then
  ! Do not convert, because origin unit of relaxtime is second.
  tunitfac=ONE
elseif(index(tunit,'hr').gt.0 .or. index(tunit,'hour').gt.0)then
  ! Convert second to hour
  tunitfac=SEC2HOUR
elseif(index(tunit,'day').gt.0)then
  ! Convert second to day
  tunitfac=SEC2DAY
elseif(index(tunit,'month').gt.0)then
  ! Convert second to month
  tunitfac=SEC2MONTH
elseif(index(tunit,'yr').gt.0 .or. index(tunit,'year').gt.0)then
  ! Convert second to year
  tunitfac=SEC2YEAR
else
  write(logunit,*)'ERROR: unrecognized time unit: ',trim(tunit)
  flush(logunit)
endif
! compute relax time
iviscoelas=0
do i=1,nmatblk_viscoelas
!  if(mat_domain(i)==VISCOELASTIC_DOMAIN .or. &
!    mat_domain(i)==VISCOELASTIC_TRINFDOMAIN .or. &
!    mat_domain(i)==VISCOELASTIC_INFDOMAIN)then
!    iviscoelas=iviscoelas+1
    imat=imatve_to_imat(i)
    ! Variables viscosity_blk and shearmod_blk are NOT nondimensionalized.
    ! Therefore the relaxtime will be in seconds.
    ! If nondimensionalized it may be better to dimensionalize again to compute
    ! the relax time without any confusion.
    ! WRONG
    !relaxtime(:,iviscoelas)=tunitfac*viscosity_blk(:,iviscoelas)/shearmod_blk(imat)
    ! APPROXIMATE
    ! CORRECT
    if(trim(devel_example).eq.'axial_rod')then
      relaxtime(:,i)=tunitfac*viscosity_blk(:,i)/ym_blk(imat)
    else  
      !relaxtime(:,i)=tunitfac*TWO*viscosity_blk(:,i)/shearmod_blk(imat)
      relaxtime(:,i)=devel_rtfac*tunitfac*viscosity_blk(:,i)/shearmod_blk(imat)
    endif
!  endif
enddo
min_relaxtime=minscal(minval(relaxtime))
max_relaxtime=maxscal(maxval(relaxtime))
if(myrank.eq.0)then
  write(logunit,'(a,g0.6,1x,a,g0.6)')'Relax time => min: ',min_relaxtime,' max: ',max_relaxtime
  write(logunit,'(a)')'Time unit: '//trim(tunit)
  flush(logunit)
endif

! count elastic and viscoelastic elements
nelmt_elas=0; nelmt_viscoelas=0
do i_elmt=1,nelmt
  mdomain=mat_domain(mat_id(i_elmt))
  ! infinite elements included in elastic domain
  if(mdomain==ELASTIC_DOMAIN .or.  &
    mdomain==ELASTIC_TRINFDOMAIN .or.  &
    mdomain==ELASTIC_INFDOMAIN)then
    nelmt_elas=nelmt_elas+1
  elseif(mdomain==VISCOELASTIC_DOMAIN .or. &
    mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
    mdomain==VISCOELASTIC_INFDOMAIN)then
    nelmt_viscoelas=nelmt_viscoelas+1
  else
    write(logunit,*)'ERROR: unrecognized material domain',mdomain,'!'
    flush(logunit)
    stop
  endif
enddo
if(nelmt/=nelmt_elas+nelmt_viscoelas)then
  write(logunit,*)'ERROR: total number of elements mismatch!'
  flush(logunit)
  stop
endif
tot_nelmt_elas=sumscal(nelmt_elas)
max_nelmt_elas=maxscal(nelmt_elas)
min_nelmt_elas=minscal(nelmt_elas)

tot_nelmt_viscoelas=sumscal(nelmt_viscoelas)
max_nelmt_viscoelas=maxscal(nelmt_viscoelas)
min_nelmt_viscoelas=minscal(nelmt_viscoelas)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0,1x,a,i0)')'elements elastic => total: ', &
  tot_nelmt_elas,' max: ',max_nelmt_elas,' min: ',min_nelmt_elas
  write(logunit,'(a,i0,1x,a,i0,1x,a,i0)')'elements viscoelastic => total: ', &
  tot_nelmt_viscoelas,' max: ',max_nelmt_viscoelas,' min: ',min_nelmt_viscoelas
  flush(logunit)
endif
allocate(eid_elas(nelmt_elas),eid_viscoelas(nelmt_viscoelas))
allocate(elas_e0(nst,ngll,nelmt_viscoelas), &
visco_q0(nst,nmaxwell,ngll,nelmt_viscoelas),q0(nst,nmaxwell))

! save element ID separately for elastic and viscoelastic elements
ielmt_elas=0; ielmt_viscoelas=0
do i_elmt=1,nelmt
  mdomain=mat_domain(mat_id(i_elmt))
  ! infinite elements included in elastic domain
  if(mdomain==ELASTIC_DOMAIN .or.  &
    mdomain==ELASTIC_TRINFDOMAIN .or.  &
    mdomain==ELASTIC_INFDOMAIN)then
    ielmt_elas=ielmt_elas+1
    eid_elas(ielmt_elas)=i_elmt
  elseif(mdomain==VISCOELASTIC_DOMAIN .or. &
    mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
    mdomain==VISCOELASTIC_INFDOMAIN)then
    ielmt_viscoelas=ielmt_viscoelas+1
    eid_viscoelas(ielmt_viscoelas)=i_elmt
  endif
enddo
! check element counts

! prepare ghost partitions for the communication
if(nproc.gt.1)then
  call prepare_ghost()
endif
! Old connectivity no longer necessary
deallocate(g_num0)
call sync_process

! prepare fault
if( iseqsource .and. (eqsource_type.eq.3 .or. eqsource_type.eq.4) )then
  if(myrank==0)then
    write(logunit,'(a)',advance='no')'preparing split fault...'
    flush(logunit)
  endif
  call prepare_fault(errcode,errtag)
  call sync_process
  call control_error(errcode,errtag,stdout,myrank)
  if(myrank==0)then
    write(logunit,*)'complete!'
    flush(logunit)
  endif
endif

! Apply displacement boundary conditions
allocate(bcnodalv(nndof,nnode))
bcnodalv=ZERO
if(myrank==0)then
  write(logunit,'(a)',advance='no')'applying BC...'
  flush(logunit)
endif
allocate(gdof(nndof,nnode),stat=istat)
if (istat/=0)then
  write(logunit,*)'ERROR: cannot allocate memory!'
  flush(logunit)
  stop
endif

allocate(infinite_iface(6,nelmt),infinite_face_idir(6,nelmt))
infinite_iface=.false.
infinite_face_idir=-9999
! activate degrees of freedom
!gdof=0
!gdof(idofu,:)=1
!
!gdof=1
call activate_dof(errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)

! This will ensure that the gdof IDs are same in the finite/infinite interface
! nodes if they lie across different processors.
! This can also be done if we explicitly define the gdof ON/OFF state on those
! inteface nodes across all the processors.
call assemble_ghosts_gdof(nndof,gdof,gdof)
where(gdof>0)gdof=1
call sync_process

! At this point, all godf IDs are consistent across the parallel interfaces
! having the values either 0 or 1.

! Apply Dirichlet boundary conditions
call apply_bc(bcnodalv,errcode,errtag)
call sync_process
call control_error(errcode,errtag,stdout,myrank)

!! Undo the unmatching dipalcement BCs. This may occur in fault implementation
!call sync_process
!call undo_unmatching_displacementBC(bcnodalv)
!call sync_process

call finalize_gdof(errcode,errtag)
call control_error(errcode,errtag,stdout,myrank)
if(myrank==0)then
  write(logunit,*)'complete!'
  flush(logunit)
endif
!-------------------------------------
! modify ghost GDOFs
if(nproc.gt.1)then
  call prepare_ghost_gdof()
endif
allocate(num(nenode),coord(ngnode,ndim),jac(ndim,ndim),deriv(ndim,nenode),     &
bmat(nst,nedofu),eld(nedofu),bload(nedofu),vload(nedofu),eload(nedofu),        &
nodalu(nndofu,nnode),egdof(nedof),egdofu(nedofu),stat=istat)
if (istat/=0)then
  write(logunit,*)'ERROR: cannot allocate memory!'
  flush(logunit)
  stop
endif
if(ISPOT_DOF)then
  allocate(nodalphi(nnode),nodalg(ndim,nnode),nodalB(ndim,nnode),stat=istat)
  if(istat/=0)then
    write(logunit,*)'ERROR: cannot allocate memory!'
    flush(logunit)
    stop
  endif
endif
!--------------------------------


!!-------------------------------------------------------------------------------
!! TEST for map_point2naturalhex8 function
!print*,hex8_gnode
!num=g_num(:,1)
!xelm=g_coord(:,num(hex8_gnode))
!xp=xelm(:,1)
!print*,'testing:'
!print*,xp
!print*,ngnode
!print*,transpose(xelm)
!print*,shape(xelm)
!call  map_point2naturalhex8(xelm,xp,xip,errcode,errtag)
!call control_error(errcode,errtag,stdout,myrank)
!print*,'map location:',xp,xip
!print*,'-----------------------------------------------------------------------'
!xp=g_coord(:,num(27))
!xelm=g_coord(:,num(hex8_gnode))
!call  map_point2naturalhex8(xelm,xp,xip,errcode,errtag)
!call control_error(errcode,errtag,stdout,myrank)
!print*,'map location:',xp,xip
!print*,'-----------------------------------------------------------------------'
!xp=g_coord(:,num(2))
!xelm=g_coord(:,num(hex8_gnode))
!call  map_point2naturalhex8(xelm,xp,xip,errcode,errtag)
!call control_error(errcode,errtag,stdout,myrank)
!print*,'map location:',xp,xip
!print*,'-----------------------------------------------------------------------'
!stop
!!-------------------------------------------------------------------------------

! store elemental global degrees of freedoms from nodal gdof
! this removes the repeated use of reshape later but it has larger size than
! gdof!!!
allocate(gdof_elmt(nedof,nelmt))
gdof_elmt=0
do i_elmt=1,nelmt
  gdof_elmt(:,i_elmt)=reshape(gdof(:,g_num(:,i_elmt)),(/nedof/))
enddo

!---------------------------------------
! global indexing
! this process is necessary for petsc implementation and is done in a single
! processor
call sync_process

if(ismpi .and. nproc.gt.1 .and. myrank.eq.0)call gindex()
call sync_process
!----------------------------------------------

tot_neq=sumscal(neq); max_neq=maxscal(neq); min_neq=minscal(neq)
if(myrank==0)then
  write(logunit,'(a,i0,a,i0,a,i0)')'degrees of freedoms => total:',tot_neq,&
                                  ' max:',max_neq,' min:',min_neq
  flush(logunit)
endif

if(myrank==0)then
  write(logunit,'(a)',advance='no')'preprocessing...'
  flush(logunit)
endif

! compute initial stress assuming elastic domain
if(savedata%stress)then
  allocate(stress_elmt(nst,ngll,nelmt),stress_nodal(nst,nnode))
  stress_elmt=ZERO
endif
if(savedata%strain)then
  allocate(strain_elmt(nst,ngll,nelmt),strain_nodal(nst,nnode))
  strain_elmt=ZERO
endif
if(isstress0)then
if(s0_type==0)then
  ! compute initial stress using SEM itself

  allocate(extload(0:neq),du(0:neq),dprecon(0:neq), &
  storekmat(nedof,nedof,nelmt),stat=istat)
  if (istat/=0)then
    write(logunit,*)'ERROR: cannot allocate memory!'
    flush(logunit)
    stop
  endif
  extload=ZERO; isgravity=.true.; ispseudoeq=.false.
  call stiffness_bodyload(nelmt,neq,hex8_gnode,g_num,gdof_elmt,mat_id,gam_blk, &
  storekmat,dprecon,extload,isgravity,ispseudoeq)

  if(myrank==0)then
    write(logunit,*)'complete!'
    flush(logunit)
  endif
  !-------------------------------

  if(myrank==0)then
    write(logunit,'(a)')'--------------------------------------------'
    flush(logunit)
  endif

  ! assemble from ghost partitions
  if(nproc.gt.1)then
    call assemble_ghosts(nndof,neq,dprecon,dprecon)
  endif

  dprecon(1:)=one/dprecon(1:); dprecon(0)=ZERO

  ! compute displacement due to graviy loading to compute initial stress
  du=ZERO
  call ksp_pcg_solver(neq,nelmt,storekmat,du,extload,   &
  dprecon,gdof_elmt,ksp_iter,errcode,errtag)
  call control_error(errcode,errtag,stdout,myrank)
  
  du(0)=ZERO

  call elastic_stress(nelmt,neq,hex8_gnode,g_num,gdof_elmt,du,stress_elmt)
  deallocate(extload,dprecon,du,storekmat)
elseif(s0_type==1)then
  ! compute initial stress using simple relation for overburden pressure
  call overburden_stress(nelmt,g_num,z_datum,s0_datum,epk0,stress_elmt)
else
  write(logunit,*)'ERROR: s0_type:',s0_type,' not supported!'
  flush(logunit)
  stop
endif
endif
!-------------------------------

allocate(slipload(0:neq),extload(0:neq),rhoload(0:neq),ubcload(0:neq))

extload=ZERO
rhoload=ZERO

allocate(node_valency(nnode))

! compute node valency only once
node_valency=0
do i_elmt=1,nelmt
  ielmt=i_elmt
  num=g_num(:,ielmt)
  node_valency(num)=node_valency(num)+1
enddo
! assemble all node_valceny across the processors
call assemble_ghosts_nodal_iscalar(node_valency,node_valency)

! open summary file
if(myrank==0)then
  write(logunit,'(a)')'KSP_MAXITER, KSP_TOL, NL_MAXITER, NL_TOL'
  write(logunit,'(i0,1x,g0.6,1x,i0,1x,g0.6)')KSP_MAXITER,KSP_RTOL,NL_MAXITER,NL_TOL
  !write(logunit,'(a)')'Number of SRFs'
  !write(logunit,'(i0)')nsrf
  write(logunit,'(a,i0)')'Number of time steps:',ntstep
  !write(logunit,'(a)')'STEP, CGITER, NLITER, UXMAX, UMAX'
  flush(logunit)
endif

allocate(cohf(nmatblk),nuf(nmatblk),phif(nmatblk),psif(nmatblk),ymf(nmatblk))

allocate(load(0:neq),bodyload(0:neq),viscoload(0:neq),             &
resload(0:neq),du(0:neq),u(0:neq),kmat(nedof,nedof),            &
storekmat(nedof,nedof,nelmt),stat=istat)
if(istat/=0)then
  write(logunit,*)'ERROR: cannot allocate memory!'
  flush(logunit)
  stop
endif

allocate(ngpart_node(nnode))

! compute stable time step for implicit integration
dt=dtstep

nl_tot=0

elas_e0=ZERO
visco_q0=ZERO

nodalu=ZERO
bodyload=ZERO
viscoload=ZERO
slipload=ZERO ! slip load
extload=ZERO ! incremental external load
ubcload=ZERO
load=ZERO
u=ZERO

!! WARNING: this is temporary
!! WARNING: must remove this because it is already withing the time loop
!! earthquake source
!call earthquake_load(neq,extload,errcode,errtag)
!call sync_process
!call control_error(errcode,errtag,stdout,myrank)

if(solver_type.eq.petsc_solver)then
  ! prepare sparsity of the stiffness matrix
  call prepare_sparse()

  ! petsc solver
  call petsc_initialize()
  if(myrank==0)then
    write(logunit,'(a)')'petsc_initialize: SUCCESS!'
    flush(logunit)
  endif
  ! create sparse vector, matrix, and preallocate                                                         
  ! TODO: following call is not necessary for RECYCLE                            
  call petsc_create_vector()                                                     
  call petsc_matrix_preallocate_size()                                           
  call petsc_create_matrix()                                                     
  call petsc_create_solver()                                                     
  if(myrank==0)then
    write(logunit,'(a)')'petsc_preallocate_matrix_size: SUCCESS!'
    flush(logunit)
  endif
endif

! WARNING: TODO
! slip gdof for split PC
if(ISDISP_DOF)then
  allocate(iseq(neq))
  iseq=.false.
  do i_node=1,nnode
    do i=1,nndofu
      geq=gdof(idofu(i),i_node)
      if(geq.gt.0)iseq(geq)=.true.
    enddo
  enddo
  nequ=count(iseq)
  allocate(gdofu(nequ))
  gdofu=-9999
  inum=0
  do i=1,neq
    if(iseq(i))then
      inum=inum+1
      gdofu(inum)=i
    endif
  enddo
  if(nequ.ne.inum)then
    if(myrank==0)then
      write(logunit,'(a)')'ERROR: nequ & inum mismatch!'
      flush(logunit)
    endif
    call close_process
  endif
  deallocate(iseq)
endif
! prepare background gravity data
call prepare_gravity()

allocate(storederiv(ndim,ngll,ngll,nelmt),storejw(ngll,nelmt))
if(solver_type.eq.builtin_solver .or. solver_diagscale)then
  allocate(dprecon(0:neq))
endif
if(solver_diagscale)then
  allocate(ndscale(0:neq))
endif
if(trim(devel_example).eq.'axial_rod')then
  open(77,file=trim(file_head)//"_strain.dat",action="write",status="replace")
endif
time_step: do i_tstep=1,ntstep
  t=dt*real(i_tstep,kreal)
  nodalu=ZERO
  if(ISPOT_DOF)then
    nodalphi=ZERO
  endif
  ubcload=ZERO
  !extload=ZERO
  rhoload=ZERO
  if(myrank==0)then
    write(logunit,'(a,i0,a,g0.6)')'step: ',i_tstep,' t: ',t
    flush(logunit)
  endif

  if(i_tstep==1)then
    ! compute elastic stiffness matrix for time = 0
    call compute_stiffness_elastic(storekmat,rhoload,errcode,errtag)
    !if(istraction.and.trcase.eq.TRACTION_INTERNAL)then
    !  call compute_surface_stiffness(storekmat,errcode,errtag)
    !  symmetric_solver=.false.
    !  call control_error(errcode,errtag,stdout,myrank)
    !endif
    if(solver_type.eq.petsc_solver)then
      call petsc_set_stiffness_matrix(storekmat)
      if(myrank==0)then
        write(logunit,'(a)')' petsc_set_stiffness_matrix: SUCCESS!'
        flush(logunit)
      endif
      call petsc_set_ksp_operator(reuse_pc=.false.)
    endif
  elseif(i_tstep==2)then
    ! Since we use a uniform dt, following routine has to be called only once 
    ! for a linear viscoelastic model. For nonlinear or nonuniform time steps
    ! it has to be called for every time steps or every changing time step.
    ! This will simply overwrite the storekmat for viscoelastic elements.
    call compute_stiffness_viscoelastic(nelmt_viscoelas,eid_viscoelas,         &
         dt,relaxtime,storekmat,errcode,errtag)
    if(solver_type.eq.petsc_solver)then
      call petsc_set_stiffness_matrix(storekmat)
      if(myrank==0)then
        write(logunit,'(a)')' petsc_set_stiffness_matrix: SUCCESS!'
        flush(logunit)
      endif
      call petsc_set_ksp_operator(reuse_pc=.true.)
    endif
  endif

  ! apply traction boundary conditions
  ! WARNING: i_tstep==1 is ONLY for rod example
  if((istraction.or.isfstraction).and.i_tstep==1)then
    if(myrank==0)then
      write(logunit,'(a)',advance='no')'applying traction...'
      flush(logunit)
    endif
    call apply_traction(extload,errcode,errtag)
    call control_error(errcode,errtag,stdout,myrank)
    if(myrank==0)then
      write(logunit,*)'complete!',maxval(abs(extload))
      flush(logunit)
    endif
  endif
  if(trim(devel_example).eq.'axial_rod')then
    if(i_tstep>600)extload=ZERO 
  endif
  ! apply magnetic traction
  if(ismtraction)then
    if(myrank==0)then
      write(logunit,'(a)',advance='no')'applying magnetic traction...'
      flush(logunit)
    endif
    call apply_mtraction(extload,errcode,errtag)
    call sync_process
    call control_error(errcode,errtag,stdout,myrank)
    if(myrank==0)then
      write(logunit,*)'complete!',maxval(abs(extload))
      flush(logunit)
    endif
  endif
  ! compute load contributed by the earthquake slip
  ! split-node apparoch: prescribe the slip on the fault explicitly
  if(iseqsource.and.eqsource_type.eq.3)then
    if(i_tstep==1)then
      if(myrank==0)then
        write(logunit,'(a)')'Earthquake source type: slip with split node'
        write(logunit,'(a,1x,i2)')'Slip taper option: ',itaper_slip
        flush(logunit)
      endif
      if(divide_slip)then
        sfac=HALF
        ! Plus side
        call compute_fault_slip_load(-1,sfac,storekmat,slipload,errcode,errtag)
        ! Minus side
        call compute_fault_slip_load(1,sfac,storekmat,slipload,errcode,errtag)
      else
        sfac=ONE
        ! Plus side
        call compute_fault_slip_load(1,sfac,storekmat,slipload,errcode,errtag)
      endif
      ! The "sync" here is very important because some processors arrive this 
      ! stage faster than other. This may hang going to control_error routine!
      call sync_process
      call control_error(errcode,errtag,stdout,myrank)
    endif
    if(srate)then
      extload=t*slipload
    else
      extload=slipload
    endif
  endif
  ! moment-density tensor apparoch: compute equivalent moment-density tensor
  ! from the prescribe slip on the fault
  if(iseqsource.and.eqsource_type.lt.3.and.i_tstep==1)then
    if(myrank==0)then
      write(logunit,'(a)')'Earthquake source type: moment-density tensor'
      flush(logunit)
    endif
    call earthquake_load(neq,extload,errcode,errtag)
    call sync_process
    call control_error(errcode,errtag,stdout,myrank)
  endif
  
  ! Modify RHS vector for prescribed displacements
  ! WARNING: need to check for nedofu
  do i_elmt=1,nelmt
    ielmt=i_elmt ! all elements
    num=g_num(:,ielmt)
    egdof=gdof_elmt(:,ielmt)
   
    kmat=storekmat(:,:,ielmt)
    iedof=0
    do j=1,nenode
      do i=1,nndof !nndofu
        iedof=iedof+1
        if(bcnodalv(i,num(j))/=ZERO)then
          ubcload(egdof)=ubcload(egdof)-kmat(:,iedof)*bcnodalv(i,num(j))
        endif
      enddo
    enddo
  enddo ! i_elmt
  if(solver_type.eq.builtin_solver)then
    ! Compute diagonal precoditioner
    dprecon=ZERO
    do i_elmt=1,nelmt
      ielmt=i_elmt ! all elements
      egdof=gdof_elmt(:,ielmt)
      do j=1,nedof
        dprecon(egdof(j))=dprecon(egdof(j))+storekmat(j,j,ielmt)
      enddo
    enddo ! i_elmt
    dprecon(0)=ZERO

    ! assemble diagonal preconditioner
    if(nproc.gt.1)then
      call assemble_ghosts(nndof,neq,dprecon,dprecon)
    endif
    call sync_process()
    dprecon(0)=ZERO
    nzero_dprecon=count(dprecon(1:).eq.ZERO)
    if(nzero_dprecon.gt.0)then
      write(logunit,*)'WARNING: nzero in dprecon:', &
      nzero_dprecon,minval(abs(dprecon(1:))),maxval(abs(dprecon(1:)))
      flush(logunit)
    endif

    if(solver_diagscale)then
      ! regularize linear equations,
      ndscale=ONE
      do i=1,neq
        ndscale(i)=one/sqrt(abs(dprecon(i)))
      enddo
      ! nondimensionalize stiffness matrix
      ! elastic region
      do i_elmt=1,nelmt_elas
        ielmt=eid_elas(i_elmt)
        egdof=gdof_elmt(:,ielmt)
        do i=1,nedof
          do j=1,nedof
            storekmat(i,j,ielmt)=ndscale(egdof(i))*storekmat(i,j,ielmt)*       &
            ndscale(egdof(j))
          enddo
        enddo
      enddo
      ! viscoelastic region
      do i_elmt=1,nelmt_viscoelas
        ielmt=eid_viscoelas(i_elmt)
        egdof=gdof_elmt(:,ielmt)
        do i=1,nedof
          do j=1,nedof
            storekmat(i,j,ielmt)=ndscale(egdof(i))*storekmat(i,j,ielmt)*       &
            ndscale(egdof(j))
          enddo
        enddo
      enddo
    else
      dprecon(1:)=one/dprecon(1:)
    endif
  endif !(solver_type.eq.builtin_solver)

  extload(0)=ZERO
  ! set BC nodal displacements to nodalu array
  if(ISDISP_DOF)then
    do i=1,nndofu
      idof=idofu(i)
      do j=1,nnode
        if(bcnodalv(idof,j)/=ZERO)nodalu(i,j)=bcnodalv(idof,j)
      enddo
    enddo
  endif
  ! set BC nodal potential to nodalphi array
  if(ISPOT_DOF)then
    do i=1,nndofphi
      idof=idofphi(i)
      do j=1,nnode
        if(bcnodalv(idof,j)/=ZERO)nodalphi(j)=bcnodalv(idof,j)
      enddo
    enddo
  endif
  ksp_tot=0; nl_iter=0

  ! only for fault
  load=extload+ubcload+rhoload

  du=ZERO; u=ZERO
  
  load(0)=ZERO

  !bodyload=ZERO; bodyload(0)=ZERO
  ! nonlinear iteration loop
  nonlinear: do i_nliter=1,NL_MAXITER
    nl_iter=nl_iter+1
    resload=load-bodyload
    resload(0)=ZERO
    maxres=maxscal(maxval(abs(resload)))
    if(myrank==0)then
      write(logunit,'(a,i0,1x,g0.6,1x,g0.6)')' Residual NL: ',i_nliter, &
      maxres,maxval(abs(bodyload))
      flush(logunit)
    endif

    if(solver_type.eq.builtin_solver)then
      ! builtin solver
      if(solver_diagscale)then
        ! nondimensionlize load
        resload=ndscale*resload
        call ksp_cg_solver(neq,nelmt,storekmat,du,resload,     &
        gdof_elmt,ksp_iter,errcode,errtag)
        du=ndscale*du
        call control_error(errcode,errtag,stdout,myrank)
      else
        ! pcg solver
        call ksp_pcg_solver(neq,nelmt,storekmat,du,resload,    &
        dprecon,gdof_elmt,ksp_iter,errcode,errtag)
        call control_error(errcode,errtag,stdout,myrank)
      endif
    else
      !petsc solver
      !call petsc_set_stiffness_matrix(storekmat)
      !if(myrank==0)print*,'petsc_set_stiffness_matrix: SUCCESS!'
      call petsc_set_vector(resload)
      if(myrank==0)then
        write(logunit,'(a)')' petsc_set_vector: SUCCESS!'
        flush(logunit)
      endif
      !call petsc_set_ksp_operator()

      call petsc_solve(du(1:),ksp_iter,ksp_convreason)
      if(myrank==0)then
        write(logunit,'(a)')' petsc_solve: SUCCESS!'
        flush(logunit)
      endif
    endif

    ksp_tot=ksp_tot+ksp_iter
    du(0)=ZERO
    maxdu=maxscal(maxval(abs(du)))
    if(myrank==0)then
      write(logunit,'(a,i0,1x,a,g0.6)')' KSP iters: ',ksp_iter, &
      'max du: ',maxdu
      if(solver_type.eq.petsc_solver)then
        write(logunit,'(a,i0)')' convergence reason: ',ksp_convreason
      endif
      flush(logunit)
    endif
    u=u+du

    maxu=maxscal(maxval(abs(u)))
    ! check convergence
    uerr=ZERO
    if(maxu.eq.ZERO)then
      uerr=one
    else
      uerr=maxdu/maxu
    endif
    nl_isconv=uerr.le.NL_TOL
    if(i_nliter>1.and.maxscal(maxval(abs(resload))).le.ZEROtol)nl_isconv=.true.
    if(myrank==0)then
      write(logunit,'(a,g0.6,1x,a,g0.6)')' UErr:',maxdu/maxu,'maxu:',maxu
      flush(logunit)
    endif

    call sync_process()
    ! update total nodal solution vector
    ! time steps are not incremental!!
    ! therefore, NO u(t+1)=u(t)+du
    ! u contains both diaplacement and/or gravity
    ! displacement
    if(ISDISP_DOF)then
      do i=1,nndofu
        idof=idofu(i)
        do i_node=1,nnode
          if(gdof(idof,i_node)/=0)then
            nodalu(i,i_node)=u(gdof(idof,i_node))
          endif
        enddo
      enddo
    endif
    ! gravity
    if(ISPOT_DOF)then
      do i=1,nndofphi
        idof=idofphi(i)
        do i_node=1,nnode
          if(gdof(idof,i_node)/=0)then
            ! \phi is a scalar
            nodalphi(i_node)=u(gdof(idof,i_node))
          endif
        enddo
      enddo
    endif
    bodyload=ZERO; !viscoload=ZERO

    if(ISDISP_DOF)then
      if(myrank==0)then
        write(logunit,*)'computing elemental stress'
        flush(logunit)
      endif
      ! compute stress
      ! elastic elements
      ! This part is repeated for the first step. We should change this for
      ! efficiency.
      do i_elmt=1,nelmt_elas
        ielmt=eid_elas(i_elmt)
        num=g_num(:,ielmt)
        egdofu=gdof_elmt(edofu,ielmt)
        eld=reshape(nodalu(:,g_num(:,ielmt)),(/nedofu/))
        bload=ZERO
        do i=1,ngll ! integration loop
          call compute_cmat_elastic(bulkmod_elmt(i,ielmt), &
          shearmod_elmt(i,ielmt),cmat)
          
          deriv=storederiv(:,:,i,ielmt)
          jacw=storejw(i,ielmt)
        
          call compute_bmat_stress(deriv,bmat)
          estrain=matmul(bmat,eld)
          sigma=matmul(cmat,estrain)
          if(savedata%strain)strain_elmt(:,i,ielmt)=estrain
          if(savedata%stress)stress_elmt(:,i,ielmt)=sigma

          eload=matmul(sigma,bmat)
          bload=bload+eload*jacw
        enddo ! i
        bodyload(egdofu)=bodyload(egdofu)+bload
      enddo
      !---------------------------------------------------------------------------
      
      if(allelastic)exit nonlinear

      ! viscoelastic elemenets
      do i_elmt=1,nelmt_viscoelas
        ielmt=eid_viscoelas(i_elmt)
        imat=mat_id(ielmt)
        imatve=imat_to_imatve(imat)

        muratio=muratio_blk(:,imatve)
        tratio=dt/relaxtime(:,imatve)
        num=g_num(:,ielmt)
        egdofu=gdof_elmt(edofu,ielmt)
        eld=reshape(nodalu(:,g_num(:,ielmt)),(/nedofu/))
        bload=ZERO; vload=ZERO
        do i=1,ngll ! integration loop
          K=bulkmod_elmt(i,ielmt) 
          G=shearmod_elmt(i,ielmt)
          !call compute_cmat_maxwell(K,G,tratio,cmat)
          !call compute_cmat_elastic(K,G,cmat)
          deriv=storederiv(:,:,i,ielmt)
          jacw=storejw(i,ielmt)
        
          call compute_bmat_stress(deriv,bmat)
          estrain=matmul(bmat,eld) ! strain at current time step
          if(savedata%strain.and.i_tstep==1.and.i_nliter==1)then
            ! store elastic strain
            strain_elmt(:,i,ielmt)=estrain
          endif
          trace_strain=estrain(1)+estrain(2)+estrain(3)
          dev_strain(1:3)=(estrain(1:3)-ONE_THIRD*trace_strain)
          dev_strain(4:6)=estrain(4:6)*HALF
          if(savedata%stress.and.i_tstep==1.and.i_nliter==1)then
            ! store elastic stress
            esigma=TWO*G*dev_strain
            esigma(1:3)=esigma(1:3)+K*trace_strain
            stress_elmt(:,i,ielmt)=esigma
          endif
          !esigma=matmul(cmat,estrain)
          !esigma=TWO*G*dev_strain
          !esigma(1:3)=esigma(1:3)+K*trace_strain

          !----------------------------ZIENCKIEWICZ---------------------------
          e0=elas_e0(:,i,i_elmt)
          q0=visco_q0(:,:,i,i_elmt)
          if(i_tstep==1)then !.and.i_nliter==1)then
            ! initialize
            e0=dev_strain
            do i_maxwell=1,nmaxwell                                                
              q0(:,i_maxwell)=e0                                                     
            enddo
          endif
          ! compute stress and update state variables
          call compute_stress_genmaxwell(K,G,tratio,muratio,estrain,e0,q0,&
          vesigma,vsigma)
          ! compute load
          !eload=matmul(vsigma,bmat)
          !vload=vload+eload*jacw
          eload=matmul(vesigma,bmat)
          bload=bload+eload*jacw
          ! update
          if(nl_isconv.or.nl_iter==NL_MAXITER)then
            elas_e0(:,i,i_elmt)=e0
            visco_q0(:,:,i,i_elmt)=q0
            if(savedata%strain)strain_elmt(:,i,ielmt)=estrain
            if(savedata%stress)stress_elmt(:,i,ielmt)=vesigma
          endif
          !----------------------------ZIENCKIEWICZ---------------------------
        enddo ! i
        !----compute the total bodyload vector----
        !viscoload(egdofu)=viscoload(egdofu)+vload
        bodyload(egdofu)=bodyload(egdofu)+bload
      enddo ! i_elmt

      bodyload(0)=ZERO
      !viscoload(0)=ZERO
    endif !(ISDISP_DOF)

    ! time step 0  and i_nliter 0 is entirely elastic
    ! write data for tiem step 0
    if(i_tstep==1.and.i_nliter==1)then
      if(ISDISP_DOF)then
        ! write displacement field
        if(savedata%disp)then
          call write_vector_to_file(nnode,DIM_L*nodalu,&
          ext='dis',istep=0)
          ! On the free surface
          if(savedata%fsplot)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_L*nodalu(:,gnode_fs),&
            ext='dis',istep=0)
          endif
          if(savedata%fsplot_plane)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_L*nodalu(:,gnode_fs),&
            ext='dis',istep=0,plane=.true.)
          endif
        endif
        ! plot stress
        if(savedata%stress)then
          call compute_nodal_tensor(stress_elmt,stress_nodal)  
          if(nproc.gt.1)then
            call assemble_ghosts_nodal_vectorn(NST,stress_nodal,stress_nodal)
          endif
          ! compute average on the sharing nodes
          do i_comp=1,NST
            stress_nodal(i_comp,:)=stress_nodal(i_comp,:)/real(node_valency,kreal)
          enddo
          call write_vector_to_file(nnode,DIM_MOD*stress_nodal,&
          ext='sig',istep=0)
          ! On the free surface
          if(savedata%fsplot)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
            ext='sig',istep=0)
          endif
          if(savedata%fsplot_plane)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
            ext='sig',istep=0,plane=.true.)
          endif
        endif
        !if(devel_mgll)then
        !  call compute_save_density_perturbation(nodalu,errcode,errtag)
        !endif
        ! save strain at a point
        if(savedata%strain)then
          call compute_nodal_tensor(strain_elmt,strain_nodal)  
          if(nproc.gt.1)then
            call assemble_ghosts_nodal_vectorn(NST,strain_nodal,strain_nodal)
          endif
          ! compute average on the sharing nodes
          do i_comp=1,NST
            strain_nodal(i_comp,:)=strain_nodal(i_comp,:)/real(node_valency,kreal)
          enddo
          call write_vector_to_file(nnode,strain_nodal,&
          ext='eps',istep=0)
          ! On the free surface
          if(savedata%fsplot)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*strain_nodal(:,gnode_fs),&
            ext='eps',istep=0)
          endif
          if(savedata%fsplot_plane)then
            call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*strain_nodal(:,gnode_fs),&
            ext='eps',istep=0,plane=.true.)
          endif
          write(77,*)0.0,strain_nodal(1,2099)                               
          flush(77)
        endif
        ! Benchmark calculation for elastic result
        if(benchmark_okada .and. ISDISP_DOF)then
          call compute_okada_solution()
        endif
      endif
      if(ISPOT_DOF)then
        ! plot gravity potential
        if(savedata%gpot)then
          call write_scalar_to_file(nnode,DIM_GPOT*nodalphi,ext='gpot',istep=0) 
          ! On the free surface
          if(savedata%fsplot)then
            call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs), &
            ext='gpot',istep=0) 
          endif
          if(savedata%fsplot_plane)then
            call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs), &
            ext='gpot',istep=0,plane=.true.) 
          endif
        endif
        if(savedata%mpot)then
          call write_scalar_to_file(nnode,DIM_MPOT*nodalphi,ext='mpot',istep=0)
          ! On the free surface
          if(savedata%fsplot)then
            call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs), &
            ext='mpot',istep=0) 
          endif
          if(savedata%fsplot_plane)then
            call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs), &
            ext='mpot',istep=0,plane=.true.) 
          endif
        endif
        ! gravitational
        if(savedata%agrav)then
          ! compute acceleration due to gravity
          call compute_gradient_of_scalar(nodalphi,nodalg)
          if(nproc.gt.1)then
            call assemble_ghosts_nodal_vector(nodalg,nodalg)
          endif
          ! compute average on the sharing nodes
          do i_comp=1,ndim
            nodalg(i_comp,:)=nodalg(i_comp,:)/real(node_valency,kreal)
          enddo
          ! plot gravity accelration
          if(savedata%agrav)then
            call write_vector_to_file(nnode,DIM_G*nodalg,ext='grav',istep=0)
            ! On the free surface
            if(savedata%fsplot)then
              call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs), &
              ext='grav',istep=0)
            endif
            if(savedata%fsplot_plane)then
              call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs), &
              ext='grav',istep=0,plane=.true.)
            endif
          endif
        endif
        ! magnetic
        if(savedata%magb)then
          ! compute magnetic field
          call compute_premagnetic_field(nodalphi,nodalB)
          if(nproc.gt.1)then
            call assemble_ghosts_nodal_vector(nodalB,nodalB)
          endif
          ! compute average on the sharing nodes
          do i_comp=1,ndim
            nodalB(i_comp,:)=nodalB(i_comp,:)/real(node_valency,kreal)
          enddo
          ! multiply by \mu_0
          nodalB=MAG_CONS*nodalB
          ! plot magnetic field
          if(savedata%magb)then
            call write_vector_to_file(nnode,DIM_B*nodalB,ext='magb',istep=0)
            ! plot magnetic field on the free surface
            if(savedata%fsplot)then
              call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs), &
              ext='magb',istep=0)
            endif
            if(savedata%fsplot_plane)then
              call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs), &
              ext='magb',istep=0,plane=.true.)
            endif
          endif
        endif
      endif
      
      if(ntstep.le.1.and.NL_MAXITER.le.1)then
        exit time_step 
      endif
    endif ! i_tstep==1.and.i_nliter==NL_MAXITER
    ! Exit nonlinear loop if converged
    if(nl_isconv)exit nonlinear
  enddo nonlinear ! i_nliter=1,NL_MAXITER

  !!compute nodal strain--------------------------------------------------
  !if(ISDISP_DOF)then
  !  ! strain 
  !  if(savedata%strain)then
  !    if(myrank==0)then
  !      write(logunit,*)'computing nodal strain'
  !      flush(logunit)
  !    endif
  !    call compute_nodal_tensor(strain_elmt,strain_nodal)  
  !    if(nproc.gt.1)then
  !      call assemble_ghosts_nodal_vectorn(NST,strain_nodal,strain_nodal)
  !    endif
  !    ! compute average on the sharing nodes
  !    do i_comp=1,NST
  !      strain_nodal(i_comp,:)=strain_nodal(i_comp,:)/real(node_valency,kreal)
  !    enddo
  !    write(77,*)dt*real(i_tstep),strain_nodal(1,2099)                               
  !    flush(77)
  !  endif
  !endif
  !-----------------------------------------------------------------------------

  if(nl_iter>=NL_MAXITER .and. .not.nl_isconv)then
    if(myrank==0)then
      write(logunit,*)'WARNING: nonconvergence in nonlinear iterations!'
      write(logunit,*)'desired tolerance:',NL_TOL,' achieved tolerance:',uerr
      flush(logunit)
    endif
  endif
  nl_tot=nl_tot+nl_iter

  if(ISDISP_DOF)then
    ! plot displacement
    if(savedata%disp)then
      call write_vector_to_file(nnode,DIM_L*nodalu,ext='dis',istep=i_tstep) 
      ! On the free surface
      if(savedata%fsplot)then
        call write_vector_to_file_freesurf(nnode_fs,DIM_L*nodalu(:,gnode_fs),&
        ext='dis',istep=i_tstep)
      endif
      if(savedata%fsplot_plane)then
        call write_vector_to_file_freesurf(nnode_fs,DIM_L*nodalu(:,gnode_fs),&
        ext='dis',istep=i_tstep,plane=.true.)
      endif
    endif
    ! plot stress
    if(savedata%stress)then
      call compute_nodal_tensor(stress_elmt,stress_nodal)  
      if(nproc.gt.1)then
        call assemble_ghosts_nodal_vectorn(NST,stress_nodal,stress_nodal)
      endif
      ! compute average on the sharing nodes
      do i_comp=1,NST
        stress_nodal(i_comp,:)=stress_nodal(i_comp,:)/real(node_valency,kreal)
      enddo
      call write_vector_to_file(nnode,DIM_MOD*stress_nodal,&
      ext='sig',istep=i_tstep)
      ! On the free surface
      if(savedata%fsplot)then
        call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
        ext='sig',istep=i_tstep)
      endif
      if(savedata%fsplot_plane)then
        call write_vector_to_file_freesurf(nnode_fs,DIM_MOD*stress_nodal(:,gnode_fs),&
        ext='sig',istep=i_tstep,plane=.true.)
      endif
    endif
    ! plot strain
    if(savedata%strain)then
      call compute_nodal_tensor(strain_elmt,strain_nodal)  
      if(nproc.gt.1)then
        call assemble_ghosts_nodal_vectorn(NST,strain_nodal,strain_nodal)
      endif
      ! compute average on the sharing nodes
      do i_comp=1,NST
        strain_nodal(i_comp,:)=strain_nodal(i_comp,:)/real(node_valency,kreal)
      enddo
      call write_vector_to_file(nnode,strain_nodal,&
      ext='eps',istep=i_tstep)
      ! On the free surface
      if(savedata%fsplot)then
        call write_vector_to_file_freesurf(nnode_fs,strain_nodal(:,gnode_fs),&
        ext='eps',istep=i_tstep)
      endif
      if(savedata%fsplot_plane)then
        call write_vector_to_file_freesurf(nnode_fs,strain_nodal(:,gnode_fs),&
        ext='eps',istep=i_tstep,plane=.true.)
      endif
    endif
    !if(devel_mgll)then
    !  call compute_save_density_perturbation(nodalu,errcode,errtag)
    !endif
  endif
  if(ISPOT_DOF)then
    ! plot gravity potential
    if(savedata%gpot)then
      call write_scalar_to_file(nnode,DIM_GPOT*nodalphi,ext='gpot',istep=i_tstep) 
      ! On the free surface
      if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs), &
        ext='gpot',istep=i_tstep) 
      endif
      if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs,DIM_GPOT*nodalphi(gnode_fs), &
        ext='gpot',istep=i_tstep,plane=.true.) 
      endif
    endif
    if(savedata%mpot)then
      call write_scalar_to_file(nnode,DIM_MPOT*nodalphi,ext='mpot',istep=i_tstep) 
      ! On the free surface
      if(savedata%fsplot)then
        call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs), &
        ext='mpot',istep=i_tstep) 
      endif
      if(savedata%fsplot_plane)then
        call write_scalar_to_file_freesurf(nnode_fs,DIM_MPOT*nodalphi(gnode_fs), &
        ext='mpot',istep=i_tstep,plane=.true.) 
      endif
    endif
    
    ! gravitational
    if(savedata%agrav)then
      ! compute acceleration due to gravity
      call compute_gradient_of_scalar(nodalphi,nodalg)
      if(nproc.gt.1)then
        call assemble_ghosts_nodal_vector(nodalg,nodalg)
      endif
      ! compute average on the sharing nodes
      do i_comp=1,ndim
        nodalg(i_comp,:)=nodalg(i_comp,:)/real(node_valency,kreal)
      enddo
      ! plot gravity accelration
      if(savedata%agrav)then
        call write_vector_to_file(nnode,DIM_G*nodalg,ext='grav',istep=i_tstep)
        ! On the free surface
        if(savedata%fsplot)then
          call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs), &
          ext='grav',istep=i_tstep)
        endif
        if(savedata%fsplot_plane)then
          call write_vector_to_file_freesurf(nnode_fs,DIM_G*nodalg(:,gnode_fs), &
          ext='grav',istep=i_tstep,plane=.true.)
        endif
      endif
    endif
    ! magnetic
    if(savedata%magb)then
      ! compute magnetic field
      call compute_premagnetic_field(nodalphi,nodalB)
      if(nproc.gt.1)then
        call assemble_ghosts_nodal_vector(nodalB,nodalB)
      endif
      ! compute average on the sharing nodes
      do i_comp=1,ndim
        nodalB(i_comp,:)=nodalB(i_comp,:)/real(node_valency,kreal)
      enddo
      ! multiply by \mu_0
      nodalB=MAG_CONS*nodalB
      ! plot magnetic field
      if(savedata%magb)then
        call write_vector_to_file(nnode,DIM_B*nodalB,ext='magb',istep=i_tstep)
        ! On the free surface
        if(savedata%fsplot)then
          call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs), &
          ext='magb',istep=i_tstep)
        endif
        if(savedata%fsplot_plane)then
          call write_vector_to_file_freesurf(nnode_fs,DIM_B*nodalB(:,gnode_fs), &
          ext='magb',istep=i_tstep,plane=.true.)
        endif
      endif
    endif
  endif
  if(myrank==0)then
    write(logunit,*)' ' 
    flush(logunit)
  endif
enddo time_step ! i_tstep time stepping loop
if(savedata%strain)then
  close(77)
  deallocate(strain_elmt,strain_nodal)
endif
! cleanup solver
if(solver_type.eq.petsc_solver)then
  call petsc_destroy_vector()                                                      
  call petsc_destroy_matrix()                                                      
  call petsc_destroy_solver()                                                      
  call petsc_finalize()
endif

call cleanup_fault()
deallocate(egdof,egdofu)
if(allocated(gdofu))deallocate(gdofu)
deallocate(extload,load,resload,rhoload,ubcload)
deallocate(du,u)
deallocate(nodalu,bcnodalv)
if(ISPOT_DOF)then
  deallocate(nodalphi)
endif

if(solver_type.eq.builtin_solver .and.solver_diagscale)then
  deallocate(dprecon,ndscale)
endif
deallocate(mat_id,mat_domain,gam_blk,ym_blk,coh_blk,nu_blk,phi_blk,psi_blk,srf)
if(allocated(imat_to_imatve))deallocate(imat_to_imatve)
if(allocated(imatve_to_imat))deallocate(imatve_to_imat)
deallocate(g_coord,g_num)
deallocate(node_valency)
deallocate(bmat,deriv,eld,num)
deallocate(kmat)
if(allocated(infinite_iface))deallocate(infinite_iface)
if(allocated(infinite_face_idir))deallocate(infinite_face_idir)
call cleanup_ghost()

return
end subroutine specfem3d
!===============================================================================
