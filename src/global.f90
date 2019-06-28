! DESCRIPTION
!  These modules contain global variables, constants, and parameters. 
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Jul 07,2011; HNG, Apr 09,2010
! TODO

!  Precision parameters
module set_precision
implicit none
integer,parameter :: kreal=8 !selected_real_kind(15)
end module set_precision
!===============================================================================

! This module contains math constants and parameters
module math_constants
use set_precision
implicit none
real(kind=kreal),parameter :: ZERO=0.0_kreal,HALF=0.5_kreal,ONE=1.0_kreal,     &
TWO=2.0_kreal,THREE=3.0_kreal,FOUR=4.0_kreal,ONE_THIRD=ONE/THREE,              &
TWO_THIRD=TWO/THREE,FOUR_THIRD=FOUR/THREE
real(kind=kreal),parameter :: PI=3.141592653589793_kreal,HALFPI=HALF*PI, &
TWOPI=TWO*PI

! tolerance value for zero
real(kind=kreal),parameter :: INFTOL=1.0e32_kreal,ZEROTOL = 1.0e-12_kreal

! Gravitational constant: G ( m^3 kg^{-1} s^{-2} )
! source: 2014 CODATA recommended values
! http://www.physics.nist.gov/cgi-bin/cuu/Value?bg
real(kind=kreal),parameter :: GRAV_CONS=6.67408e-11_kreal
! vacuum permeability, permeability of free space, permeability of vacuum,
! or magnetic constant
real(kind=kreal),parameter :: MAG_CONS=FOUR*PI*1.0e-7_kreal
end module math_constants
!===============================================================================

! This module contains conversion parameters for some physical quantities
module conversion_constants
use set_precision
use math_constants,only:ONE,PI
implicit none
! CGS to SI units
real(kind=kreal),parameter :: M2KM=1e-3_kreal
real(kind=kreal),parameter :: CGS2SI_MOMENT=1e-7_kreal
! SI to CGS units
real(kind=kreal),parameter :: KM2M=1e+3_kreal
real(kind=kreal),parameter :: SI2CGS_MOMENT=1e+7_kreal
real(kind=kreal),parameter :: DEG2RAD=PI/180.0_kreal,RAD2DEG=180.0_kreal/PI
! Time unit conversion
! 1 min   = 60 secs
! 1 hour  = 60 mins
! 1 day   = 24 hours
! 1 month = 30 days
! 1 year  = 365 days
real(kind=kreal),parameter :: SEC2MIN   = ONE/60.0_kreal,        &
                              SEC2HOUR  = SEC2MIN/60.0_kreal,    &
                              SEC2DAY   = SEC2HOUR/24.0_kreal,   &
                              SEC2MONTH = SEC2DAY/30.0_kreal,    &
                              SEC2YEAR  = SEC2DAY/365.0_kreal

end module conversion_constants
!===============================================================================

! This module contains KSP solver parameters
module ksp_constants
use set_precision
implicit none

! Following two parameters can be replaced via input file
integer :: KSP_MAXITER=20000
real(kind=kreal) :: KSP_RTOL=1.0e-8_kreal
real(kind=kreal),parameter :: KSP_ATOL=1.0e-30_kreal
real(kind=kreal),parameter :: KSP_DTOL=1.0e30_kreal
end module ksp_constants
!===============================================================================

! (Non)dimensionalize
module dimensionless
use set_precision
implicit none
! density
real(kind=kreal) :: DIM_DENSITY
real(kind=kreal) :: NONDIM_DENSITY
! mass
real(kind=kreal) :: DIM_M
! time
real(kind=kreal) :: DIM_T
real(kind=kreal) :: NONDIM_T
! length,displacement
real(kind=kreal) :: DIM_L
real(kind=kreal) :: NONDIM_L
! velocity
real(kind=kreal) :: DIM_VEL
! acceleration
real(kind=kreal) :: DIM_ACCEL
real(kind=kreal) :: NONDIM_ACCEL
! elastic modulus
real(kind=kreal) :: DIM_MOD
real(kind=kreal) :: NONDIM_MOD
! moment tensor
real(kind=kreal) :: DIM_MTENS
real(kind=kreal) :: NONDIM_MTENS
! gravity potential
real(kind=kreal) :: DIM_GPOT
! gravity acceleration
real(kind=kreal) :: DIM_G
! magnetic potential
real(kind=kreal) :: DIM_MPOT
! magnetic field or induction
real(kind=kreal) :: DIM_B
end module dimensionless
!===============================================================================

! This model contains Earth related constants
module earth_constants
use set_precision
use math_constants,only: ONE,PI
! R_EARTH is the radius of the bottom of the oceans (radius of Earth in m)
double precision, parameter :: R_EARTH = 6371000.d0
real(kind=kreal),parameter :: R_UNIT_SPHERE = ONE
end module earth_constants
!===============================================================================

! This module global parameters/variables
module global
use set_precision
use math_constants,only: ONE
implicit none

integer(kind=8) :: GPU_pointer

! UTM projection zone (Optional input)
integer :: UTM_ZONE
character(len=20) :: proc_str,ptail,ptail_inp
logical,parameter :: off=.false., on=.true.
character(len=3) :: method
integer,parameter :: NDIM=3 !number of dimensions. this is a 3D version
integer,parameter :: NDIM2=NDIM*NDIM

! degrees of freedoms
logical :: ISPOT_DOF,ISDISP_DOF
integer,parameter :: PGRAVITY=1,PMAGNETIC=2
! potential type: 1: 'gravity' and 2: 'magnetic'
integer :: POT_TYPE
character(len=20) :: POT_STRING
! background gravity
! if ISPOT_DOF is .FALSE., this serves as the Cowling approximation
logical :: ISGRAV0

! number of displacement degrees of freedom per node - ux, uy, uz
integer,parameter :: nndofu=3
! number of gravity potential degrees of freedom per node - \phi
integer,parameter :: nndofphi=1
! number of displacement degrees of freedom per node - ux, uy, uz
integer :: nndof
!IDs for displacement degrees of freedom per node
integer,dimension(NDIM) :: idofu
!IDs for gravity potential degrees of freedom per node
integer,dimension(nndofphi) :: idofphi
!IDs for displacement degrees of freedom per element
integer,allocatable :: edofu(:)
!IDs for gravity potential degrees of freedom per element
integer,allocatable :: edofphi(:)
! number of gauss-lobatto-legendre points along the x, y, and z directions.
! we always take ngllx=nglly=ngllz for the simplicity
! ngll=ngllx*nglly*ngllz
integer :: ngllx,nglly,ngllz,ngll
! total number of GLL points on XY, YZ, and ZX planes
integer :: ngllxy,ngllyz,ngllzx,maxngll2d
integer,parameter :: ng=8 ! number of gauss points for FEM

logical :: ismpi !.true. : MPI, .false. : serial
integer :: myrank,nproc !myrank is indexed from 0
integer :: ngdof !Number of nodal degrees of freedom per processor = nndof*nnode
integer :: neq !number of equations per processor = ngdof - degrees of freedom
! lost due to constraints
integer,allocatable :: l2gdof(:)!map from local dof (in processor) to global dof
! (in entire system). l2gdof contains the global (system-wide) indices for the 
! dof in that processor.

integer :: nsparse
integer,allocatable :: kcol_sparse(:),krow_sparse(:),kgcol_sparse(:),          &
kgrow_sparse(:)

integer :: nenode,nedof ! number of elemental nodes (ie. number of nodes per
!element = ngll), number of elemental degrees of freedom -> nedof=nndof*nenode

integer :: ngnode ! number of geometrical nodes. usually, for FEM ngnode=nenode
! total # of nodes, ie. gll pts, per processor, and number of elements per
! processor
integer :: nnode,nelmt
integer,allocatable :: mat_id(:)

integer,allocatable :: g_num0(:,:),g_num(:,:),gdof(:,:),ggdof(:,:)
!g_num: global node IDs for each element (per processor).
!gdof: matrix of nodal dof (per processor).
!eg. gdof(:,1) gives all dof in 1st node of that processor.
!ggdof is the same as gdof, but for all processors (global global dof matrix).
integer,allocatable :: gdof_elmt(:,:)
!gdof_elmt: matrix of elemental degrees of freedom (per processor), e.g.,
!gdof_elmt(:,1) gives all dof in 1st element in the processsor.

! number of elemental degrees of freedoms for displacement
integer :: nedofu
! number of elemental degrees of freedoms for gravity
integer :: nedofphi

! acceleration due to gravity
real(kind=kreal),parameter :: agrav=9.82_kreal
real(kind=kreal),allocatable :: g_coord(:,:) ! global coordinates
! coordinate extents of partitioned model of the finite region
real(kind=kreal) :: pmodel_minx,pmodel_maxx,pmodel_miny,pmodel_maxy, &
pmodel_minz,pmodel_maxz
! coordinate extents of model
real(kind=kreal) :: model_minx,model_maxx,model_miny,model_maxy, &
model_minz,model_maxz
real(kind=kreal) :: mincoord,maxcoord,absmaxcoord ! maximum value of coordinates
real(kind=kreal) :: maxsize_elmt,sqmaxsize_elmt ! maximum size of the element across the digonal

! model properties
! bulk modulus, shear modulus, mass density, magnetization
logical :: isbulkmod,isshearmod,ismassdens,ismagnetization
! minimum, maximum value of density
real(kind=kreal) :: mindensity,maxdensity
! minimum, maximum value of bulk modulus
real(kind=kreal) :: minbulkmod,maxbulkmod
! minimum, maximum value of shear modulus
real(kind=kreal) :: minshearmod,maxshearmod
integer :: nmatblk !number of material domains
logical :: isdensity
! Flag to check if the domain is empty
logical,allocatable :: isempty_blk(:)
integer,allocatable :: mat_domain(:),type_blk(:)
real(kind=kreal),allocatable :: gam_blk(:),rho_blk(:),ym_blk(:),nu_blk(:),     &
coh_blk(:),phi_blk(:),psi_blk(:)
character(len=60),allocatable :: mfile_blk(:)
real(kind=kreal),allocatable :: bulkmod_blk(:),shearmod_blk(:)
real(kind=kreal),allocatable :: massdens_elmt(:,:),bulkmod_elmt(:,:),          &
shearmod_elmt(:,:)
real(kind=kreal),allocatable :: grav0_nodal(:,:),dgrav0_elmt(:,:,:)
! magnetization
real(kind=kreal),allocatable :: magnetization_elmt(:,:,:)
integer :: nwmat
integer,allocatable :: waterid(:)
logical,allocatable :: water(:)

integer,parameter :: ELASTIC_DOMAIN=1,VISCOELASTIC_DOMAIN=11
integer,parameter :: ACOUSTIC_DOMAIN=2
! ELASTIC_TRINFDOMAIN=100
! ELASTIC_INFDOMAIN=1000
! VISCOELASTIC_TRINFDOMAIN=1100
! VISCOELASTIC_INFDOMAIN=11000
!   Transition/Infinite domain element can have either only one DOF type,
!   i.e., gravity perturbation or both DOF types, i.e., gravity perturbation 
!   and displacement
integer,parameter :: ELASTIC_TRINFDOMAIN=100*ELASTIC_DOMAIN, &
                     ELASTIC_INFDOMAIN=1000*ELASTIC_DOMAIN, &
                     VISCOELASTIC_TRINFDOMAIN=100*VISCOELASTIC_DOMAIN, &
                     VISCOELASTIC_INFDOMAIN=1000*VISCOELASTIC_DOMAIN
integer,parameter :: nmaxwell=1
integer :: nmatblk_viscoelas
integer,allocatable :: imat_to_imatve(:),imatve_to_imat(:)
real(kind=kreal),allocatable :: muratio_blk(:,:),viscosity_blk(:,:)

integer :: nmatblk_magnet
integer,allocatable :: imat_to_imatmag(:),imatmag_to_imat(:)
real(kind=kreal),allocatable :: magnetization_blk(:,:),Mmag_blk(:)
logical,allocatable :: ismagnet_blk(:)
! model types
character(len=20) :: model_type

! custome model
character(len=20) :: cmodel
! reference z coordinate from which depth is measured down positive
real(kind=kreal) :: cmodel_zref
! alpha parameter, density anomaly on the top
real(kind=kreal) :: cmodel_alpha,cmodel_drho0

real(kind=kreal),allocatable :: storederiv(:,:,:,:)
real(kind=kreal),allocatable :: storejw(:,:)

logical :: allelastic,iseqload,iseqsource,iswater,phinu
! pseudostatic coefficients for earthquake loading eqkh=ah/g, eqkv=av/g
real(kind=kreal) :: eqkx,eqky,eqkz
! where ah and av are horizontal and vertical pseudostatic accelerations

! Number of unique stress components. Stress tensor is symmetric and its
! components are in order:
! 1: sxx
! 2: syy
! 3: szz
! 4: sxy
! 5: syz
! 6: szx
! First three components are axial components. Last three components are
! shear components.
integer,parameter :: NST=6
character(len=250) :: file_head,inp_path,out_path,part_path
! displacement BC, ghost, traction, and water surface files
character(len=250) :: confile,idfile
character(len=250),dimension(3) :: coordfile
! Essential boundary conditions (Primary)
logical :: isubc ! truncated BC. 0: No, 1: yes [Default]
! Displacement BC defined on the all unique SEM nodes on the free surface.
! 0: No [Default], 1: yes
logical :: isfsubc
character(len=250) :: ufspath
character(len=250) :: matfile,uxfile,uyfile,uzfile,ufsfile,gfile, &
fsfile,wsfile
character(len=250) :: matfile_inherit

! Traction (Natural or Secondary) boundary conditions
logical :: istraction,ismtraction
logical :: isfstraction,isfstraction0
character(len=250) :: trfspath,trfspath0
character(len=250) :: trfsfile,trfsfile0
character(len=250) :: trfile,mtrfile
! Main traction cases
integer,parameter :: TRACTION_EXTERNAL=0,TRACTION_INTERNAL=1,TRACTION_BOTH=2
integer :: trcase
! infinite BC file
logical :: infbc ! Infinite BC. 0: No [Default], 1: yes
logical :: zmaxinf ! Infinite BC on the zmax surface. 0: No [Default] 1: Yes
logical,allocatable :: infinite_iface(:,:)
integer,allocatable :: infinite_face_idir(:,:)
character(len=250) :: inffile
integer :: nelmt_finite,nelmt_trinfinite,nelmt_infinite
integer,allocatable :: g_num_finite(:,:),g_num_trinfinite(:,:),g_num_infinite(:,:)
integer,allocatable :: elmt_finite(:),elmt_trinfinite(:),elmt_infinite(:) 
integer :: nnode_finite,nnode_trinfinite,nnode_infinite
integer,allocatable :: node_finite(:),node_trinfinite(:),node_infinite(:) 
integer :: nl_maxiter
real(kind=kreal) :: nl_tol
integer :: nexcav,ninc,nsrf,ntstep
real(kind=kreal) :: dtstep
! time unit: 'second', 'minute', 'hour', 'day', 'month', 'year'
character(len=10) :: tunit
! Excavation ID (regions), number
! of excavation IDs (regions) in each stage
! Excavation ID (regions), number of excavation IDs (regions) in each stage
integer,allocatable :: excavid(:),nexcavid(:)
! strength reduction factors
real(kind=kreal),allocatable :: srf(:)

! dervative = 0 contraints
logical :: isdxval,isdyval,isdzval

integer :: ielmtINF1,ielmtINF2,ifaceINF
! material identification type for the infinite-elements
! 'define': material block ID/s are specified by imat_inf (and) imat_trinf 
! 'inherit': material block ID/s will be inherited from the parent element 
character(len=10) :: matinf_type
integer :: imat_trinf,imat_inf
character(len=250) :: infrfile ! infinite-element layer reference file
! radius for the infinite-element mesh layer
real(kind=kreal) :: rinf
! BC value for the infinite-element layer surface, 
real(kind=kreal) :: valINF
! switch to check if traction is zero on the infinite surfaces 
logical :: isdxINF,isdyINF,isdzINF
! Pole coordinates for the decay function
! which is generally zero
real(kind=kreal) :: pole_coord0(NDIM),pole_coord1(NDIM)
! add infinite element layer internally
logical :: add_infmesh
! switch to check if any of the axis has to be fixed for infinite elements
logical :: isyfixINF,isxfixINF,iszfixINF

! type of initial pole for infinite element
character(len=10) :: pole0
! type of pole for infinite element
! plane, axis, point
character(len=10) :: pole_type
! pole axis
! 1: X, 2: Y, 3: Z
integer :: pole_axis
real(kind=kreal) :: axis_range(2)
! Type of integration along the infinite direction 
integer,parameter :: INFINITE_RADAU=0, INFINITE_GAUSS=1
! infinite quadrature. DEFAULT: INFINITE_RADAU
integer :: infquad

! initial stress
logical :: isstress0,usek0 ! use k0 to compute horizontal stress also for
! s0_type=0
real(kind=kreal) :: s0_type ! initial stress type
! 0: by default compute by SEM,
! 1: simple overburden pressure use s0+gamma*z
!    only for horizontal and homogeneous, and
! 2: read from file
real(kind=kreal) :: z_datum,s0_datum,epk0
! z-coordinate at free surface, stress at free surface, and
! earth-pressure coefficient

! earthquake fault variables
character(len=250) :: slipfile,cmtfile,faultfile,faultmetafile,                &
faultslipfile_plus,faultslipfile_minus
integer :: nsource,fault_npatch
!0: Fault slip, 1: CMT solution, 2: Finite fault, 3: Node split
integer :: eqsource_type
integer :: itaper_slip
!0: No taper, 1: Set slip on boundary nodes to ZERO
! Currently implemented ONLY for node split
logical :: divide_slip
! Devide slip in equal and opposite on + and - sides of the fault surface
! .false.: No (Default), .true.: Yes
logical :: srate
! uniform slip rate
! isgslip: Switch for geometrical splitting.
! isnoslip: Switch for noslip imposition. This is useful for adjoint simulation.
logical :: isgsplit,isnoslip
! station
character(len=250) :: stationfile
integer :: nstation
logical :: isstation
!Benchmarking
! .TRUE. : Okada benchmark, .FALSE. : no benchmark (Default)
logical :: benchmark_okada
logical :: benchmark_error
!parameters for okada
real(kind=kreal) :: okada_origin(NDIM) !okada_origin(3): Depth + down
real(kind=kreal) :: okada_aw1,okada_aw2,okada_al1,okada_al2,okada_depth,       &
okada_dip,okada_disl1,okada_disl2,okada_disl3

! order of viscoelastic algorithm
! 1: First order, 2: Second order Simo and Hughes, 21: Second order Zienckiewicz
integer,parameter :: VISCO_ORDER=21
integer,parameter :: VISCO_MAXWELL=0,VISCO_ZENER=1,VISCO_GENMAXWELL=2
integer :: visco_model
character(len=20) :: visco_rheology
! solver types
! DEVELOPER OPTIONS
! diagonally scale equations and use CG solver without preconditioning
logical,parameter ::  solver_diagscale=.false.
!for builtin solver
integer,parameter :: smart_solver=0  ! select appropriate solver automatically
integer,parameter :: builtin_solver=1! select builtin conjugate gradient solver
integer,parameter :: petsc_solver=2  ! select PETSC solver
integer :: solver_type=builtin_solver !smart_solver !builtin_solver !petsc_solver !smart_solver 
! By default solver is symmetric but it may be changed later depending on the
! conditions.
logical :: symmetric_solver=.true.
! save options
type savedata_options
  logical :: model,disp,stress,porep,psigma,maxtau,nsigma,scf,vmeps
  logical :: model_cell
  logical :: strain
  logical :: gpot,agrav
  logical :: mpot,magb
  logical :: infinite
  ! Free surface plot. If this option is .TRUE., and the free surface file is
  ! given, the result will be plotted on the free surface.
  logical :: fsplot,fsplot_plane 
end type savedata_options
type(savedata_options) :: savedata

! others
character(len=1),parameter :: CR=achar(13) ! carriage return to overwrite
!previous line

! format string for time step
character(len=20) :: tstep_sformat
! Log file all information
character(len=250) :: log_file
! file unit ID for log file
integer :: logunit=7
integer :: stdout=6

! developement variables
! By default model is nondimensionalized unless the "devel_nondim" is .false. 
logical :: devel_nondim
! save model properties on GLL points
logical :: devel_mgll
real(kind=kreal) :: devel_gaminf
! relaxation time factor: ONE or TWO only
real(kind=kreal) :: devel_rtfac
! example: axial_rod
character(len=20) :: devel_example
end module global
!===============================================================================
