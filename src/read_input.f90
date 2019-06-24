! DESCRIPTION
!  This module conains the routine to read main input file.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!   HNG, Jul 07,2011; HNG, Apr 09,2010
! TODO
!  - 
module input
implicit none
character(len=250),private :: myfname=' => read_input.f90'
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

! This subroutine sets the default values of input parameters.
subroutine set_default_values()
use global
use math_constants,only:inftol,zero,zerotol,ONE,HALF,TWO,THREE
implicit none

end subroutine set_default_values
!===============================================================================

! This subroutine reads the input information from a structured ASCII text file
! REVISION:
!   HNG, Jul 07,2011; HNG, Apr 09,2010
! TODO:
!   - prompt warning or error for unknown argument/s
subroutine read_input(inp_fname,errcode,errtag,ispartmesh)
use global
use ksp_constants
use math_constants,only:inftol,zero,zerotol,ONE,HALF,TWO,THREE
use conversion_constants,only:DEG2RAD
use math_library,only:magnetic_unitvec
use string_library
use infinite_layer
use recreate_fault
implicit none
integer :: i,i_line
character(len=*),intent(in) :: inp_fname
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
logical,optional,intent(in) :: ispartmesh

logical :: isfrom_partmesh
character(len=250) :: line
character(len=60) :: lineword(9)
character(len=800) ::tag
character(len=80) :: strval,token
character(len=1) :: tmp_char
character(len=80),dimension(50) :: args
integer :: id,ind,ios,narg,slen

integer :: bc_stat,preinfo_stat,mesh_stat,material_stat,control_stat,          &
eqload_stat,stress0_stat,traction_stat,mtraction_stat,water_stat,save_stat,    &
eqsource_stat,benchmark_stat,station_stat,devel_stat,mag_stat
integer :: mat_count
integer :: ielmt,i_node,inode,imat,tmp_nelmt,tmp_nnode !,mat_domain

character(len=250) :: fname
character(len=150) :: data_path,mat_path

integer,parameter :: NMAXLINE=1000
integer :: nproc_inp ! partition ID
integer :: iselastic,istat,ival,issave,nexcavid_all
integer,allocatable :: ivect(:)
real(kind=kreal) :: rval
real(kind=kreal),allocatable :: rvect(:)

logical,allocatable :: ismat(:)

! Magnetization
! Inclination (0) or latitude (1). If latitude, the inclination is obtained
! using the relation: inclination = ATAN(2*TAN(latitude))
integer :: incORlat
! magnitude of the magnetization
real(kind=kreal) :: M0
real(kind=kreal) :: inc,dec,azim

errtag="ERROR: unknown!"
errcode=-1
errsrc=trim(myfname)//' => read_input'

! determine if the routine is called from the partmesh
isfrom_partmesh=.false.
if(present(ispartmesh))then
  if(ispartmesh)then
    isfrom_partmesh=.true.
    ! for partmesh write log in the terminal
    logunit=stdout
  endif
endif

! reading main input information
if(myrank==0)then
  write(logunit,'(a)')'reading main input file...'
  flush(logunit)
endif
preinfo_stat=-1
bc_stat=-1
traction_stat=0
mtraction_stat=0
eqload_stat=0
water_stat=0
mesh_stat=-1
material_stat=-1
stress0_stat=-1
control_stat=-1
eqsource_stat=0
station_stat=0
save_stat=0
benchmark_stat=0
devel_stat=0
mag_stat=0

ISDISP_DOF=.true.
ISPOT_DOF=.false.
POT_TYPE=PGRAVITY
POT_STRING='gravity'
ISGRAV0=.false.

iselastic=0
allelastic=.false.
istraction=.false.
isfstraction=.false.
isfstraction0=.false.
trcase=TRACTION_EXTERNAL
ismtraction=.false.
iseqload=.false.
iseqsource=.false.
iswater=.false.
isstress0=.false.
usek0=.false.
phinu=.false.
add_infmesh=.false.
pole0='origin'
pole_type='point'
pole_axis=-9999
benchmark_okada=.false.
benchmark_error=.false.
isstation=.false.

! Default savedata options
savedata%model=.false.
savedata%model_cell=.false.
savedata%disp=.false.
savedata%stress=.false.
savedata%strain=.false.
savedata%porep=.false.
savedata%psigma=.false.
savedata%maxtau=.false.
savedata%nsigma=.false.
savedata%scf=.false.
! gravity potential
savedata%gpot=.false.
savedata%agrav=.false.
! magnetic potential
savedata%mpot=.false.
savedata%magb=.false.
savedata%infinite=.false.
savedata%fsplot=.false.
savedata%fsplot_plane=.false.

! Default development variables
devel_nondim=.true.
devel_gaminf=-INFTOL
devel_mgll=.false.
devel_example=''
devel_rtfac=ONE

! by default 4th column of the material list is assumed to be unit weight
! if .true., density is assumed
isdensity=.false.

s0_type=0 ! by default compute initial stress using SEM

mat_count=0

! default value
method='sem'
! input path
if(ismpi.and.nproc.gt.1)then
  inp_path='./partition/'
else
  inp_path='./input/'
endif
! output path
out_path='./output/'
! partititon path
part_path='./partition/'

eqkx=0.0_kreal
eqky=0.0_kreal
eqkz=0.0_kreal

! default KSP parameters are set in global.f90
! default NL parameters
nl_tol=zerotol; nl_maxiter=1

nexcav=0
nsrf=1 ! number of strength reduction factors
ninc=1 ! number of load increments
ntstep=1 ! number of time steps
dtstep=zero ! time step interval
tunit='sec'
isdxval=.false.
isdyval=.false.
isdzval=.false.

isubc=.true.
isfsubc=.false.

! BC value for the infinite-element layer surface, 
! which is generally zero
infbc=.false.
zmaxinf=.false.
imat_trinf=-9999
imat_inf=-9999
valINF=zero
! default pole coordinates for infinite elements
pole_coord0=ZERO
pole_coord1=ZERO
isxfixINF=.false.
isyfixINF=.false.
iszfixINF=.false.
infquad=0
! custome model
cmodel=''
model_type='block'
matinf_type=''
! open file to read
open(unit=11,file=trim(inp_fname),status='old', action='read',iostat=ios)
if (ios /= 0)then
  write(errtag,'(a)')'ERROR: input file "'//trim(inp_fname)//&
  '" cannot be opened!'
  return
endif
do
  read(11,'(a)',iostat=ios)line
  ! This will read a line and proceed to next line
  if (ios/=0)exit
  ! check for blank and comment line
  if (isblank(line) .or. iscomment(line,'#'))cycle

  ! check for line continuation
  tag=''
  do
    call last_char(line,tmp_char,ind) 
    if(tmp_char.eq.'&')line(ind:ind)=''
    tag=trim(tag)//trim(line)    
    if(tmp_char.ne.'&')exit
    read(11,'(a)',iostat=ios)line
    ! This will read a line and proceed to next line
    if(ios /= 0)then
      write(errtag,'(a)')'ERROR: line continuation incomplete!'
      return
    endif    
  enddo
  call first_token(tag,token)

  ! read pre information
  if (trim(token)=='preinfo:')then
    if(preinfo_stat==1)then
      write(errtag,*)'ERROR: copy of line type preinfo: not permitted!'
      return
    endif
    preinfo_stat=-1;
    call split_string(tag,',',args,narg)
    call seek_integer('utm_zone',ival,args,narg,istat)
    if(istat==0)then
      if(ival.lt.1 .or. ival.gt.60)then
        write(errtag,'(a,i0)')'ERROR: invalid UTM ZONE (valid range [1,60]! ',ival
        return
      endif
      UTM_ZONE=ival
    endif
    if(ismpi.or. isfrom_partmesh)then
      nproc_inp=get_integer('nproc',args,narg);
      if(isfrom_partmesh)nproc=nproc_inp
      if(nproc_inp/=nproc)then
        write(errtag,*)'ERROR: number of processors and images must be equal!',nproc,nproc_inp
        return
      endif
    endif
    call seek_string('method',strval,args,narg)
    if (.not. isblank(strval))method=trim(strval)
    
    call seek_integer('pot_dof',ival,args,narg,istat)
    if(istat==0)then
      if(ival.eq.0)then
        ISPOT_DOF=.false.
      elseif(ival.eq.1)then
        ISPOT_DOF=.true.
      else
        write(errtag,*)'ERROR: pot_dof must be 0 or 1!',ival
        return
      endif
    endif
        
    call seek_integer('disp_dof',ival,args,narg,istat)
    if(istat==0)then
      if(ival.eq.0)then
        ISDISP_DOF=.false.
      elseif(ival.eq.1)then
        ISDISP_DOF=.true.
      else
        write(errtag,*)'ERROR: disp_dof must be 0 or 1!',ival
        return
      endif
    endif
   
    if(ISPOT_DOF)then
      strval=get_string('pot_type',args,narg)
      if(trim(strval).eq.'gravity' .or. &
         trim(strval).eq.'Gravity' .or. &
         trim(strval).eq.'GRAVITY')then
        POT_TYPE=PGRAVITY
        POT_STRING='gravity'
      elseif(trim(strval).eq.'magnetic' .or. &
             trim(strval).eq.'Magnetic' .or. &
             trim(strval).eq.'MAGNETIC')then
        POT_TYPE=PMAGNETIC
        POT_STRING='magnetic'
      else
        write(errtag,*)'ERROR: unsupported "pot_type": ',trim(strval)
        return
      endif
    endif
    
    call seek_integer('grav0',ival,args,narg,istat)
    if(istat==0)then
      if(ival.eq.0)then
        ISGRAV0=.false.
      elseif(ival.eq.1)then
        ISGRAV0=.true.
      else
        write(errtag,*)'ERROR: grav0 must be 0 or 1!',ival
        return
      endif
    endif
        
    !method=get_string('method',args,narg)
    if (trim(method)=='sem')then
      ngllx=get_integer('ngllx',args,narg);
      nglly=get_integer('nglly',args,narg);
      ngllz=get_integer('ngllz',args,narg);

      ngllxy=ngllx*nglly ! total GLL points on XY plane
      ngllyz=nglly*ngllz ! total GLL points on YZ plane
      ngllzx=ngllz*ngllx ! total GLL points on ZX plane
      ngll=ngllx*nglly*ngllz ! total GLL points
    endif
    nenode=get_integer('nenode',args,narg);
    ! number of elemental degrees of reedom
    ! nedof=nndof*nenode
    if (method=='sem')then
      ! number of geometrical nodes
      ngnode=get_integer('ngnode',args,narg);
    elseif(method=='fem')then
      ngnode=nenode ! default number of geometrical nodes
    else
      write(errtag,*)'ERROR: wrong value for method!'
      return
    endif
    call seek_string('inp_path',strval,args,narg)
    if (.not. isblank(strval))inp_path=trim(strval)
    slen=len_trim(inp_path)
    if(inp_path(slen:slen)/='/')inp_path=trim(inp_path)//'/'
    if(ismpi .or. isfrom_partmesh)then
      call seek_string('part_path',strval,args,narg)
      if (.not. isblank(strval))part_path=trim(strval)
      slen=len_trim(part_path)
      if(part_path(slen:slen)/='/')part_path=trim(part_path)//'/'
    endif

    call seek_string('out_path',strval,args,narg)
    if (.not. isblank(strval))out_path=trim(strval)
    slen=len_trim(out_path)
    if(out_path(slen:slen)/='/')out_path=trim(out_path)//'/'

    preinfo_stat=1
    cycle
  endif
  ! read mesh information
  if (trim(token)=='mesh:')then
    if(mesh_stat==1)then
      write(errtag,*)'ERROR: copy of line type mesh: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    coordfile(1)=get_string('xfile',args,narg)
    coordfile(2)=get_string('yfile',args,narg)
    coordfile(3)=get_string('zfile',args,narg)
    confile=get_string('confile',args,narg)
    fsfile=get_string('fsfile',args,narg)
    idfile=get_string('idfile',args,narg)
    if(ismpi.and.nproc.gt.1)gfile=get_string('gfile',args,narg)
    !nmatblk=get_integer('nmatblk',args,narg)

    mesh_stat=1
    cycle
  endif

  ! read bc information
  if (trim(token)=='bc:')then
    if(bc_stat==1)then
      write(errtag,*)'ERROR: copy of line type bc: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    call seek_integer('ubc',ival,args,narg,istat)
    if(istat==0)then
      if(ival.eq.0)then
        isubc=.false.
      elseif(ival.eq.1)then
        isubc=.true.
      else
        write(errtag,*)'ERROR: ubc must be 0 or 1!',ival
        return
      endif
    endif
    uxfile=get_string('uxfile',args,narg)
    uyfile=get_string('uyfile',args,narg)
    uzfile=get_string('uzfile',args,narg)
    call seek_string('ufsfile',strval,args,narg)
    if(.not. isblank(strval))then
      isfsubc=.true.
      ufsfile=trim(strval)
      ufspath=get_string('ufspath',args,narg)
    endif
    call seek_real('dxval',rval,args,narg,istat)
    if(istat==0 .and. rval==zero)isdxval=.true.
    call seek_real('dyval',rval,args,narg,istat)
    if(istat==0 .and. rval==zero)isdyval=.true.
    call seek_real('dzval',rval,args,narg,istat)
    if(istat==0 .and. rval==zero)isdzval=.true.
    ! information for infinite BC-----------------------------------------------
    call seek_integer('infbc',ival,args,narg,istat)
    if(istat==0.and.ival.gt.0)infbc=.true.

    if(myrank==0)then
      if(infbc)then
        write(logunit,'(a)')' infinite elements: YES'
      else
        write(logunit,'(a)')' infinite elements: NO'
      endif
      flush(logunit)
    endif

    if(infbc)then
      ! By default no traction BC on infinite surface
      isdxINF=.false.
      isdyINF=.false.
      isdzINF=.false.

      inffile=trim(file_head)//'_infinite_recreated'
      call seek_integer('zmaxinf',ival,args,narg,istat)
      if(istat==0 .and. ival>0)zmaxinf=.true.
      call seek_integer('add_infmesh',ival,args,narg,istat)
      if(istat==0 .and. ival>0)add_infmesh=.true.
      call seek_string('pole0',strval,args,narg)
      if (.not. isblank(strval))pole0=trim(strval)
      if(trim(pole0)=='origin')then
        pole_coord0(1)=ZERO
        pole_coord0(2)=ZERO
        pole_coord0(3)=ZERO
      elseif(trim(pole0)=='center')then
        ! the pole coordinates will be reassigned to model center in
        ! semgeotech.F90
        !pole_coord0(1)=HALF*(model_minx+model_maxx)
        !pole_coord0(2)=HALF*(model_miny+model_maxy)
        !pole_coord0(3)=HALF*(model_minz+model_maxz)
      elseif(trim(pole0)=='user')then
        pole_coord0=get_real_vect('coord0',NDIM,args,narg)
        if(myrank==0)then
          write(logunit,'(a,3(1x,f0.6))')' pole_coord0: ',pole_coord0
          flush(logunit)
        endif
      else
        write(errtag,'(a)')'ERROR: pole0 "'//trim(pole0)//'" is not supported!'
        return
      endif
      if(myrank==0)then
        write(logunit,'(a)')' initial pole: '//trim(pole0)
        flush(logunit)
      endif
      call seek_string('pole_type',strval,args,narg)
      if (.not. isblank(strval))pole_type=trim(strval)
      if(trim(pole_type)=='axis'.or.trim(pole_type)=='pointaxis')then
        pole_axis=get_integer('pole_axis',args,narg)
        if(pole_axis.lt.1 .and. pole_axis.gt.3)then
          write(errtag,'(a,i0)')'ERROR: invalid "pole_axis" value!',pole_axis
          return
        endif
      endif
      if(trim(pole_type)=='pointaxis')then
        axis_range=get_real_vect('axis_range',2,args,narg)
        if(axis_range(1).gt.axis_range(2))then
          write(errtag,'(a,1x,f0.6,1x,f0.6)')'ERROR: axis_range must be in the ascending order!',axis_range
          return
        endif
        pole_coord1=get_real_vect('coord1',NDIM,args,narg)
        if(myrank==0)then
          write(logunit,'(a,3(1x,f0.6))')' pole_coord1: ',pole_coord1
          write(logunit,'(a,2(1x,f0.6))')' axis_range: ',axis_range
          flush(logunit)
        endif
      endif
      if(myrank==0)then
        write(logunit,'(a)')' pole type: '//trim(pole_type)
        write(logunit,'(a,i0)')' pole axis: ',pole_axis
        if(add_infmesh)then
          write(logunit,'(a)')' infinite mesh: INTERNAL'
        else
          write(logunit,'(a)')' infinite mesh: EXTERNAL'
        endif
        flush(logunit)
      endif
      if(add_infmesh)then
        
        ! for internal infinite meshing only the pole types point and axis are
        ! implemented
        if(trim(pole_type)=='plane')then
          write(errtag,'(a,i0)')'ERROR: for internal infinite meshing, "plane" &
          & is not yet supported!'
          return
        endif
        matinf_type=get_string('mat_type',args,narg)
        if(trim(matinf_type)=='inherit' .and. POT_TYPE==PMAGNETIC)then
          write(errtag,'(a,i0)')'ERROR: "inherit" option for infinite elements &
          & is not supported for magnetic DOF.'
          return
        endif
        if(trim(matinf_type).eq.'define')then
          call seek_integer('imat_trinf',ival,args,narg,istat)
          if(istat==0)imat_trinf=ival
          imat_inf=get_integer('imat_inf',args,narg)
        endif
        infrfile=get_string('infrfile',args,narg)
        rinf=get_real('rinf',args,narg)
        call seek_real('valinf',rval,args,narg,istat)
        if(istat==0)valINF=rval
      endif
      call seek_string('infquad',strval,args,narg)
      if (.not. isblank(strval))then
        if(trim(strval).eq.'gauss')then
          infquad=INFINITE_GAUSS
        elseif(trim(strval).eq.'radau')then
          infquad=INFINITE_RADAU
        else
          write(errtag,'(a)')'ERROR: invalid infinite quadrature:'//trim(strval)
          return
        endif
      endif
      call seek_real('dxmin',rval,args,narg,istat)
      if(istat==0 .and. rval==zero)isdxINF=.true.
      call seek_real('dymin',rval,args,narg,istat)
      if(istat==0 .and. rval==zero)isdyINF=.true.
      call seek_real('dzmin',rval,args,narg,istat)
      if(istat==0 .and. rval==zero)isdzINF=.true.
    endif
    !---------------------------------------------------------------------------
    bc_stat=1
    cycle
  endif

  ! read initial stress information
  if (trim(token)=='stress0:')then
    if(stress0_stat==1)then
      write(errtag,*)'ERROR: copy of line type stress0: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    isstress0=.true.
    call seek_integer('type',ival,args,narg,istat)
    if(istat==0)s0_type=ival
    if(s0_type==1)then
      z_datum=get_real('z0',args,narg)
      s0_datum=get_real('s0',args,narg)
    endif
    call seek_real('k0',rval,args,narg,istat)
    if(istat==0)epk0=rval
    call seek_integer('usek0',ival,args,narg,istat)
    if(istat==0)usek0=.true.
    stress0_stat=1
    cycle
  endif

  ! read traction information
  if (trim(token)=='traction:')then
    if(traction_stat==1)then
      write(errtag,*)'ERROR: copy of line type traction: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    call seek_string('trfile',strval,args,narg)
    if(.not. isblank(strval))then
      istraction=.true.
      trfile=trim(strval)
    endif
    call seek_string('trfsfile',strval,args,narg)
    if(.not. isblank(strval))then
      isfstraction=.true.
      trfsfile=trim(strval)
      trfspath=get_string('trfspath',args,narg)
    endif
    call seek_string('trfsfile0',strval,args,narg)
    if(.not. isblank(strval))then
      isfstraction0=.true.
      trfsfile0=trim(strval)
      trfspath0=get_string('trfspath0',args,narg)
    endif
    
    call seek_string('case',strval,args,narg)
    if(.not. isblank(strval))then
      if(trim(strval).eq.'external'.or.trim(strval).eq.'External')then
        trcase=0
      elseif(trim(strval).eq.'internal'.or.trim(strval).eq.'Internal')then
        trcase=1
      elseif(trim(strval).eq.'both'.or.trim(strval).eq.'Both')then
        trcase=2
      else
        write(errtag,*)'ERROR: inrecognized traction case: '//trim(strval)
        return
      endif
    endif
    traction_stat=1
    !istraction=.true.
    !print*,trfile
    cycle
  endif

  ! read magnetic traction information
  if (trim(token)=='mtraction:')then
    if(mtraction_stat==1)then
      write(errtag,*)'ERROR: copy of line type mtraction: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    mtrfile=get_string('trfile',args,narg)
    mtraction_stat=1
    ismtraction=.true.
    !print*,mtrfile
    cycle
  endif

  ! read material list
  if (trim(token)=='material:')then
    if(material_stat==1)then
      write(errtag,*)'ERROR: copy of line type material: not permitted!'
      return
    endif
    material_stat=-1
    call split_string(tag,',',args,narg)

    ! default path for material file
    mat_path=trim(inp_path)
    call seek_string('matpath',strval,args,narg)
    if (.not. isblank(strval))mat_path=trim(strval)
    slen=len_trim(mat_path)
    if(mat_path(slen:slen)/='/')mat_path=trim(mat_path)//'/'
    matfile=get_string('matfile',args,narg)
    call seek_integer('allelastic',iselastic,args,narg,istat)
    if(istat==0 .and. iselastic==1)allelastic=.true.
    call seek_integer('density',ival,args,narg,istat)
    if(istat==0 .and. ival==1)isdensity=.true.
    call seek_string('model',strval,args,narg)
    if (.not. isblank(strval))cmodel=trim(strval)
    if(trim(cmodel)=='chakravarthi')then
      cmodel_zref=get_real('zref',args,narg)
      cmodel_alpha=get_real('alpha',args,narg)
      cmodel_drho0=get_real('drho0',args,narg)
    endif

    material_stat=1
    cycle
  endif

  ! read earthquake loading information
  if (trim(token)=='eqload:')then
    if(eqload_stat==1)then
      write(errtag,*)'ERROR: copy of line type eqload: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    eqkx=get_real('eqkx',args,narg)
    eqky=get_real('eqky',args,narg)
    eqkz=get_real('eqkz',args,narg)
    eqload_stat=1
    iseqload=.true.
    cycle
  endif

  ! read earthquake source (slip/CMT) information
  if (trim(token)=='eqsource:')then
    if(eqsource_stat==1)then
      write(errtag,*)'ERROR: copy of line type eqsource: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    eqsource_type=get_integer('type',args,narg)
    if(eqsource_type==0)then
      slipfile=get_string('slipfile',args,narg)
    elseif(eqsource_type==1)then
      cmtfile=get_string('cmtfile',args,narg)
    elseif(eqsource_type==2)then
      faultfile=get_string('faultfile',args,narg)
      faultmetafile=get_string('faultmetafile',args,narg)
    elseif(eqsource_type.eq.3 .or. eqsource_type.eq.4)then
      itaper_slip=0       ! Default: No tapering of slip
      divide_slip=.false. ! Default: No
      srate=.false. ! Default: No
      isgsplit=.false. ! Default: No
      isnoslip=.false. ! Default: No
      faultslipfile_plus=get_string('faultslipfile_plus',args,narg)
      faultslipfile_minus=get_string('faultslipfile_minus',args,narg)
      call seek_integer('taper',ival,args,narg,istat)
      if(istat==0)itaper_slip=ival
      call seek_integer('shalf',ival,args,narg,istat)
      if(istat==0 .and. ival.eq.1)divide_slip=.true.
      call seek_integer('srate',ival,args,narg,istat)
      if(istat==0 .and. ival.eq.1)srate=.true.
      call seek_integer('gsplit',ival,args,narg,istat)
      if(istat==0 .and. ival.eq.1)isgsplit=.true.
      call seek_integer('noslip',ival,args,narg,istat)
      if(istat==0 .and. ival.eq.1)isnoslip=.true.
    else
      write(errtag,'(a,i2,a)')'ERROR: eqsource_type: ',eqsource_type,' is not &
      &supported!'
      return
    endif

    eqsource_stat=1
    iseqsource=.true.
    cycle
  endif

  !read benchmark information.
  !For now (9/16) only okada benchmarking is supported
  if (trim(token)=='benchmark:')then
      if(benchmark_stat==1)then
          write(errtag,*)'ERROR: copy of line type okada: not permitted!'
          return
      endif
      call split_string(tag,',',args,narg)
      call seek_integer('okada',ival,args,narg,istat)
      if(istat==0.and.ival.gt.0)benchmark_okada=.true.
      call seek_integer('error',ival,args,narg,istat)
      if(istat==0.and.ival.gt.0)benchmark_error=.true.
      benchmark_stat=1
      cycle
  endif

  ! read water surface information
  if (trim(token)=='water:')then
    if(water_stat==1)then
      write(errtag,*)'ERROR: copy of line type water: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    wsfile=get_string('wsfile',args,narg)
    water_stat=1
    iswater=.true.
    cycle
  endif

  ! read control information
  if (trim(token)=='control:')then
    if(control_stat==1)then
      write(errtag,*)'ERROR: copy of line type control: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    call seek_real('ksp_tol',rval,args,narg,istat)
    if(istat==0)KSP_RTOL=rval
    call seek_integer('ksp_maxiter',ival,args,narg,istat)
    if(istat==0)KSP_MAXITER=ival
    call seek_real('nl_tol',rval,args,narg,istat)
    if(istat==0)nl_tol=rval
    call seek_integer('nl_maxiter',ival,args,narg,istat)
    if(istat==0)nl_maxiter=ival
    ! time step variables
    call seek_real('dt',rval,args,narg,istat)
    if(istat==0)dtstep=rval
    call seek_string('tunit',strval,args,narg)
    if (.not. isblank(strval))tunit=trim(strval)
    call seek_integer('ntstep',ival,args,narg,istat)
    if(istat==0)ntstep=ival
    ! strength reduction variables
    call seek_integer('nsrf',ival,args,narg,istat)
    if(istat==0)nsrf=ival
    allocate(srf(nsrf),rvect(nsrf))
    srf=1.0_kreal
    call seek_real_vect('srf',rvect,nsrf,args,narg,istat)
    if(nsrf>1 .and. istat/=0)then
      write(errtag,'(a)')'ERROR: argument "srf" not found or insufficient list &
      &for "srf"!'
      return
    endif
    if(istat==0)srf=rvect
    deallocate(rvect)
    if(minval(srf)<=zero)then
      write(errtag,'(a)')'ERROR: "srf" must be positive!'
      return
    endif
    ! excavation variables
    call seek_integer('nexcav',ival,args,narg,istat)
    if(istat==0)nexcav=ival
    if(nexcav>0)then
      allocate(nexcavid(nexcav))
      nexcavid=1 ! default -> 1 region in each stage
      allocate(ivect(nexcav))
      call seek_integer_vect('nexcavid',ivect,nexcav,args,narg,istat)
      if(istat==0)nexcavid=ivect
      deallocate(ivect)
      nexcavid_all=sum(nexcavid)
      allocate(excavid(nexcavid_all))
      excavid=get_integer_vect('excavid',nexcavid_all,args,narg)
    endif
    if(nsrf>1 .and. nexcav>1)then
      write(errtag,'(a)')'ERROR: cannot run slope stabiliy and excavation &
      &simultaneously!'
      return
    endif
    !---------------------------

    call seek_integer('ninc',ival,args,narg,istat)
    if(istat==0)ninc=ival
    call seek_integer('phinu',ival,args,narg,istat)
    if(istat==0.and.ival==1)phinu=.true.
    control_stat=1
    cycle
  endif

  ! read station information
  if (trim(token)=='station:')then
    if(station_stat==1)then
      write(errtag,*)'ERROR: copy of line type station: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    stationfile=get_string('sfile',args,narg)

    station_stat=1
    isstation=.true.
    cycle
  endif

  ! read save options
  if (trim(token)=='save:')then
    if(save_stat==1)then
      write(errtag,*)'ERROR: copy of line type save: not permitted!'
      return
    endif
    save_stat=-1
    call split_string(tag,',',args,narg)
    call seek_integer('model',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%model=.true.
    call seek_integer('cell',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%model_cell=.true.
    if(ISDISP_DOF)then
      call seek_integer('disp',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%disp=.true.
      call seek_integer('stress',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%stress=.true.
      call seek_integer('strain',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%strain=.true.
      call seek_integer('porep',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%porep=.true.
      call seek_integer('psigma',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%psigma=.true.
      call seek_integer('nsigma',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%nsigma=.true.
      call seek_integer('scf',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%scf=.true.
      call seek_integer('vmeps',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%vmeps=.true.
    endif
    if(ISPOT_DOF)then
      if(POT_TYPE==PGRAVITY)then
        call seek_integer('gpot',issave,args,narg,istat)
        if(istat==0 .and. issave==1)savedata%gpot=.true.
        call seek_integer('agrav',issave,args,narg,istat)
        if(istat==0 .and. issave==1)savedata%agrav=.true.
      endif
      if(POT_TYPE==PMAGNETIC)then
        call seek_integer('mpot',issave,args,narg,istat)
        if(istat==0 .and. issave==1)savedata%mpot=.true.
        call seek_integer('magb',issave,args,narg,istat)
        if(istat==0 .and. issave==1)savedata%magb=.true.
      endif
    endif
    if(infbc)then
      call seek_integer('inf',issave,args,narg,istat)
      if(istat==0 .and. issave==1)savedata%infinite=.true.
    endif
    call seek_integer('fsplot',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%fsplot=.true.
    call seek_integer('fsplot_plane',issave,args,narg,istat)
    if(istat==0 .and. issave==1)savedata%fsplot_plane=.true.

    save_stat=1
    cycle
  endif

  ! read development vaiables if any
  if (trim(token)=='devel:')then
    if(devel_stat==1)then
      write(errtag,*)'ERROR: copy of line type eqsource: not permitted!'
      return
    endif
    call split_string(tag,',',args,narg)
    
    call seek_integer('nondim',ival,args,narg,istat)
    if(istat==0 .and. ival.eq.0)devel_nondim=.false.
    if(devel_nondim .and. POT_TYPE==PMAGNETIC)then
      devel_nondim=.false.
      if(myrank==0)then
        write(logunit,'(a)')'WARNING: nondimentionalize option for magnetic &
        &DOF is currently disabled! The option "devel-nondim" was reset to &
        &".false.".'
      endif
    endif
    call seek_integer('mgll',ival,args,narg,istat)
    if(istat==0 .and. ival.eq.1)devel_mgll=.true.
    call seek_real('gaminf',rval,args,narg,istat)
    if(istat==0)devel_gaminf=rval
    call seek_string('example',strval,args,narg)
    if (.not. isblank(strval))devel_example=trim(strval)
    call seek_string('solver',strval,args,narg)
    if (.not. isblank(strval))then
      if(trim(strval).eq.'builtin' .or. trim(strval).eq.'BUILTIN')then
        solver_type=builtin_solver
      elseif(trim(strval).eq.'petsc' .or. trim(strval).eq.'PETSC')then
        solver_type=petsc_solver
      else
        write(errtag,'(a)')'ERROR: unsupported solver: '//trim(strval)//'!'
        return
      endif
    endif

    call seek_real('rtfac',rval,args,narg,istat)
    if(istat==0)devel_rtfac=rval
    if(devel_rtfac.ne.ONE .and. devel_rtfac.ne.TWO)then
      write(errtag,'(a)')'ERROR: invalid "rtfac" value! &
      &It must be either 1.0 or 2.0!'
      return
    endif
    if(myrank==0)then
      write(logunit,'(a,f0.6)')' relaxation time factor: ',devel_rtfac
      flush(logunit)
    endif

    devel_stat=1
    cycle
  endif
  write(errtag,'(a)')'ERROR: invalid line type: "'//trim(token)//'"!'
  return

enddo ! do

if(.not.iswater)savedata%porep=.false.

if(.not.ISDISP_DOF .and. .not.ISPOT_DOF)then
  write(errtag,'(a)')'ERROR: at least one DOF type must be active!'
  return
endif

! check input status
if (preinfo_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read pre information! make sure the line &
  &with "preinfo:" token is correct.'
  return
endif

! check input status
if (mesh_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read mesh information! make sure the line &
  &with "mesh:" token is correct.'
  return
endif

! check output status
if (bc_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read BC information! make sure the line &
  &with "bc:" token is correct.'
  return
endif

! check material status
if (material_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read material information! make sure the &
  &line with "material:" token is correct.'
  return
endif

! check control status
if (control_stat /= 1)then
  write(errtag,'(a)')'ERROR: cannot read control information! make sure the &
  &line with "control:" token is correct.'
  return
endif
if(myrank==0)then
  write(logunit,'(a)')'complete!'
  flush(logunit)
endif
!--------------------------------------------------------
! set data path
if(ismpi.and.nproc.gt.1)then
  data_path=trim(part_path)
else
  data_path=trim(inp_path)
endif

if(myrank==0)then
  write(logunit,'(a)')'reading mesh & material IDs...'
  flush(logunit)
endif
! read coordinates information
do i=1,ndim
  fname=trim(data_path)//trim(coordfile(i))//trim(ptail_inp)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif
  read(11,*)tmp_nnode
  if(i==1)then
    nnode=tmp_nnode
    allocate(g_coord(ndim,nnode))
  endif
  if(tmp_nnode/=nnode)then
    write(errtag,'(a)')'ERROR: total number of nodes mismatch!'
    return
  endif
  if(ismpi.and.nproc.gt.1)then
    do i_node=1,nnode
      read(11,*)inode,g_coord(i,inode)
    enddo
  else
    do i_node=1,nnode
      read(11,*)g_coord(i,i_node)
    enddo
  endif
enddo
close(11)
! Read connectivity
fname=trim(data_path)//trim(confile)//trim(ptail_inp)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)nelmt
allocate(g_num(nenode,nelmt))
if(ismpi.and.nproc.gt.1)then
  do i=1,nelmt
    read(11,*)ielmt,g_num(:,ielmt)
  enddo
else
  do i=1,nelmt
    read(11,*)g_num(:,i)
  enddo
endif
close(11)

! Recreate faultslipfile 
if(iseqsource.and.eqsource_type.eq.3)then
  if(isfrom_partmesh .or. .not.ismpi .or. (ismpi.and.nproc.eq.1))then
    call recreate_faultslip_file(faultslipfile_plus)
    call recreate_faultslip_file(faultslipfile_minus)
  endif
endif

! Read material id
fname=trim(data_path)//trim(idfile)//trim(ptail_inp)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)tmp_nelmt
if(tmp_nelmt/=nelmt)then
  write(errtag,'(a)')'ERROR: total number of elements mismatch!'
  return
endif
allocate(mat_id(nelmt))
if(ismpi.and.nproc.gt.1)then
  do i=1,nelmt
    read(11,*)ielmt,mat_id(ielmt)
  enddo
else
  do i=1,nelmt
    read(11,*)mat_id(i)
  enddo
endif
close(11)

if(myrank==0)then
  write(logunit,'(a)')'complete!'
  flush(logunit)
endif

if(myrank==0)then
  write(logunit,'(a)')'reading material properties...'
  flush(logunit)
endif
! read material lists
! for partmesh library following information is read in the library itself
! set data path
fname=trim(mat_path)//trim(matfile)
if(ismpi.and.nproc.gt.1)then
  if(trim(matinf_type).eq.'inherit')then
    fname=trim(mat_path)//trim(matfile)//'_inherit'
  endif
endif
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR1: file "'//trim(fname)//'" cannot be opened!'
  return
endif
read(11,*)
read(11,*)nmatblk
allocate(isempty_blk(nmatblk),mat_domain(nmatblk),type_blk(nmatblk), &
gam_blk(nmatblk),rho_blk(nmatblk),ym_blk(nmatblk),coh_blk(nmatblk),  &
nu_blk(nmatblk),phi_blk(nmatblk),psi_blk(nmatblk),water(nmatblk))
allocate(mfile_blk(nmatblk))
! initilize
isempty_blk=.false.
type_blk=-10000
mfile_blk=""
gam_blk=-inftol
rho_blk=-inftol
ym_blk=-inftol
nu_blk=-inftol
coh_blk=-inftol
phi_blk=-inftol
psi_blk=-inftol
allocate(ismat(nmatblk))
ismat=.false.
do i=1,nmatblk
  ! This will read a line and proceed to next line
  ! if the input line is long enough only the part of the lineword will be
  ! filled
  lineword=""
  read(11,'(a)',iostat=ios)line
  read(line,*,iostat=ios)lineword
  imat=str2int(lineword(1))
  mat_domain(imat)=str2int(lineword(2))
  type_blk(imat)=str2int(lineword(3))
  if(type_blk(imat).eq.0)then
    ! block material properties
    if(isdensity)then
      rho_blk(imat)=str2real(lineword(4))
      gam_blk(imat)=agrav*rho_blk(imat)
    else
      gam_blk(imat)=str2real(lineword(4))
      rho_blk(imat)=gam_blk(imat)/agrav
    endif

    ym_blk(imat)=str2real(lineword(5))
    nu_blk(imat)=str2real(lineword(6))
    phi_blk(imat)=str2real(lineword(7))
    coh_blk(imat)=str2real(lineword(8))
    psi_blk(imat)=str2real(lineword(9))

    if(rho_blk(imat).eq.ZERO .and. ym_blk(imat).eq.ZERO)then
      isempty_blk(imat)=.true.
    endif
    ismat(imat)=.true.
  elseif(type_blk(imat).eq.-1)then
    ! tomographic model defined on regular structured grid
    mfile_blk(imat)=trim(lineword(4))
    ismat(imat)=.true.
    ! NOTE: this block is not EMPTY!
  else
    print*,'ERROR: type_blk:',type_blk(imat),' is unsupported!'
    stop
  endif
enddo
if(count(.not.ismat).gt.0)then
  write(errtag,'(a,a)')'ERROR: some material blocks are undefined!'//new_line('a'), &
                         'HINT: please check material list file!'
  return
endif
deallocate(ismat)
nmatblk_viscoelas=count(mat_domain.eq.VISCOELASTIC_DOMAIN.or. &
                        mat_domain.eq.VISCOELASTIC_TRINFDOMAIN.or. &
                        mat_domain.eq.VISCOELASTIC_INFDOMAIN)
if(myrank==0)then
  write(logunit,'(a,i0)')' Number of empty blocks: ',count(isempty_blk)
  write(logunit,'(a,i0)')' Number of viscoelastic blocks: ',nmatblk_viscoelas
  flush(logunit)
endif
allocate(muratio_blk(nmaxwell,nmatblk_viscoelas), &
viscosity_blk(nmaxwell,nmatblk_viscoelas))
! read visoelastic properties
if(nmatblk_viscoelas.gt.0)then
  read(11,*)visco_rheology
  if(trim(visco_rheology)=='Maxwell' .or. & 
     trim(visco_rheology)=='maxwell')then
    visco_model=VISCO_MAXWELL
  elseif(trim(visco_rheology)=='Zener' .or. & 
         trim(visco_rheology)=='zener')then
    visco_model=VISCO_ZENER
  elseif(trim(visco_rheology)=='General_Maxwell' .or. & 
         trim(visco_rheology)=='general_maxwell')then
    visco_model=VISCO_GENMAXWELL
  else
    print*,'ERROR: invalid viscoelastic rheology: ',trim(visco_rheology)
    stop
  endif
  allocate(imat_to_imatve(nmatblk),imatve_to_imat(nmatblk_viscoelas))
  imat_to_imatve=-9999
  imatve_to_imat=-9999
  do i=1,nmatblk_viscoelas
    read(11,*)imat,muratio_blk(:,i),viscosity_blk(:,i)
    imat_to_imatve(imat)=i
    imatve_to_imat(i)=imat
    if(mat_domain(imat).ne.VISCOELASTIC_DOMAIN .and. &
      mat_domain(imat).ne.VISCOELASTIC_TRINFDOMAIN .and. &
      mat_domain(imat).ne.VISCOELASTIC_INFDOMAIN)then
      write(errtag,'(a,i0,a,a)')'ERROR: material block ID: ',imat,&
        ' is not viscoelastic!'//new_line('a'), &
        'HINT: please define viscoelastic properties only for viscoelastic domains!'
      return

    endif

  enddo
endif
! test IDs
if(minval(mat_id)<1 .or. maxval(mat_id)>nmatblk)then
  write(errtag,'(a)')'ERROR: material IDs must be consistent with the &
  &defined material regions!'
  print*,minval(mat_id),maxval(mat_id),nmatblk
  return
endif
! water table
water=.false.
if(iswater)then
  read(11,*,iostat=ios)nwmat
  if( ios /= 0 ) then
    write(errtag,'(a)')'ERROR: water IDs cannot be read from material list!'
    return
  endif
  allocate(waterid(nwmat))
  do i=1,nwmat
    read(11,*)id
    water(id)=.true.
    waterid(i)=id
  enddo
endif

! read magnetization information
if(POT_TYPE==PMAGNETIC)then
  do i_line=1,NMAXLINE
    read(11,'(a)',iostat=ios)line
    ! This will read a line and proceed to next line
    if (ios/=0)exit
    ! check for blank and comment line
    if (isblank(line) .or. iscomment(line,'#'))cycle

    call first_token(line,token)
    if (trim(token)=='magnetization:')then
      if(mag_stat==1)then
        write(errtag,*)'ERROR: copy of line type "magnetization:" not permitted!'
        return
      endif
      read(11,*)nmatblk_magnet
      allocate(imat_to_imatmag(nmatblk),imatmag_to_imat(nmatblk_magnet))
      allocate(Mmag_blk(nmatblk_magnet),magnetization_blk(NDIM,nmatblk_magnet))
      allocate(ismagnet_blk(nmatblk))
      Mmag_blk=ZERO
      imat_to_imatmag=-9999
      imatmag_to_imat=-9999
      ismagnet_blk=.false.
      do i=1,nmatblk_magnet
        read(11,*)imat,M0,incORlat,rval,dec,azim
        imat_to_imatmag(imat)=i
        imatmag_to_imat(i)=imat
        if(incORlat==0)then
          inc=rval*DEG2RAD
        elseif(incORlat==1)then
          inc=atan(TWO*tan(rval*DEG2RAD))
        else
          write(errtag,'(a,i0)')'ERROR: incORlat must be either 0 or 1 but has a &
          &value: ',incORlat
          return
        endif
        dec=dec*DEG2RAD
        azim=azim*DEG2RAD
        magnetization_blk(:,i)=M0*magnetic_unitvec(inc,dec,azim)
        Mmag_blk(i)=M0
        ismagnet_blk(imat)=.true.
      enddo
      
      mag_stat=1
      cycle
    endif

  enddo
  ! check magnetic status
  if (mag_stat /= 1)then
    write(errtag,'(a)')'ERROR: cannot read magnetization information! make sure &
    &the "magnetization:" information is added in the material list file.'
    return
  endif
endif

! infinite elements 
if(infbc)then
  if(isfrom_partmesh .or. .not.ismpi .or. (ismpi.and.nproc.eq.1))then
    if(add_infmesh)then
      if(trim(pole0)=='center')then
        write(errtag,'(a)')'ERROR: "center" pole cannot be used for internal &
                            &infinite meshing! Use "origin" or "user" option!'
        return
      endif
      ! Add infinite element layer
      call add_infinite_elements(errcode,errtag)
    else
      ! Reorder infinite elements 
      call reorder_infinite_elements(errcode,errtag)
    endif
    if(errcode.lt.0)return
  endif
endif

! Compute bulk and shear modulus.
! It is done after creating infinite elements because there may be inherited
! blocks for infinite elements.
allocate(bulkmod_blk(nmatblk),shearmod_blk(nmatblk))
bulkmod_blk=ZERO
shearmod_blk=ZERO
do i=1,nmatblk
  if(type_blk(i).eq.0)then
    ! convert to bulk modulus and shear modulus
    bulkmod_blk(i)=ym_blk(i)/(THREE*(one-two*nu_blk(i)))
    shearmod_blk(i)=half*ym_blk(i)/(one+nu_blk(i))
    !if(myrank==0)print*,'bulkmod & shearmod: ',bulkmod_blk(i),shearmod_blk(i)
  endif
enddo

errcode=0
if(myrank==0)then
  write(logunit,'(a)')'complete!'
  flush(logunit)
endif

end subroutine read_input
!===============================================================================
end module input
!===============================================================================
