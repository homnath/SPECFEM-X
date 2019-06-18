module earthquake
contains
!-------------------------------------------------------------------------------

! DESCRIPTION
!  This routine computes the load contributed by the earthquake moment-tensor
!  source. The moement tensor is read either directly from CMTSOLUTION or
!  computed from the slip information file.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
!  Leah Langer, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO
!  - read and implement CMTSOLUTION
!  - check for the sources shared among elements/processors
subroutine earthquake_load(lneq,load,errcode,errtag)
use dimensionless
use global
use element,only:hex8_gnode
use math_constants
use conversion_constants
use math_library,only:cross_product,determinant,sqdistance,invert,  &
IsPointInHexahedron,norm,vector_rotateZ
use string_library
use shape_library,only : dshape_function_hex8p
use gll_library,only : gll_lagrange3d_point,zwgljd
use cmtsolution
use map_location
use elastic,only:compute_cmat_elastic
#if (USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif
implicit none
integer,intent(in) :: lneq
real(kind=kreal),intent(inout) :: load(0:lneq)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer,parameter :: nmax_line=100 ! maximum number of lines in the eqsource
!file
integer :: i_elmt,i_line,ielmt,imid,inum,ios,istat,ix
integer :: num(nenode),egdofu(nedofu)

real(kind=kreal) :: coord(ndim,8),jac(ndim,ndim),xp(ndim),xip(ndim)
real(kind=kreal) :: detjac

integer :: i_w,i_l,i_p,i_src,i0,i1,i2,i3,i4,i0_nextrow
integer :: npatch_length,npatch_width,npatch
real(kind=kreal) :: pfac,rk,rl,rw
real(kind=kreal),dimension(3) :: wx1,wx2
real(kind=kreal) :: fault_area,fault_length,fault_width
! components of vectors: center, strike, normal
real(kind=kreal),dimension(3) :: fault_center
real(kind=kreal),dimension(3,1) :: fault_normal,fault_slip,fault_strike
real(kind=kreal),dimension(3) :: fault_strike90,fault_dip,vecN
real(kind=kreal) :: proj
real(kind=kreal) :: located_x(ndim),source_x(ndim),source_xi(ndim)
real(kind=kreal) :: Mmat(ndim,ndim),M(1,9),M_V(nst,1), &
local_M(3,3)
! M0_V=(n s) is an unsymmetric matrix. Therefore, we need to store all, i.e.,
! 9 components NOT 6!
real(kind=kreal) :: M0_V(NDIM2,1)
real(kind=kreal) :: deriv(ndim,ngll),derivmat(9,nedofu),fload(1,nedofu)
real(kind=kreal) :: mzz,myz,mxz,mxy,myy,mxx
logical :: is_fcenter,is_fstrike,is_fnormal,is_fslip,is_flength,is_fwidth,     &
is_fnplength,is_fnpwidth
logical :: is_mzz,is_myy,is_mxx,is_myz,is_mxz,is_mxy
logical :: is_fnpatch

integer :: mdomain
real(kind=kreal) :: bulkmod,shearmod
real(kind=kreal) :: cmat(nst,nst)
real(kind=kreal) :: cmatU(nst,NDIM2)

logical :: is_located
integer :: this_src_located,total_src_located
integer,allocatable:: isrc_located(:)
integer :: niter
real(kind=kreal) :: errx,errxd,minerr

integer :: np,np_l,np_w
real(kind=kreal),dimension(3) :: vecOA,vecOB,vecOC,vecOD,vecOC_norm
real(kind=kreal),dimension(3) :: f_cx1,f_cx2,f_cx3,f_cx4

real(kind=kreal),allocatable :: source_coord(:,:)
real(kind=kreal),allocatable :: M_cmt(:,:)

integer,allocatable :: p_connect(:,:)
real(kind=kreal),allocatable :: p_coord(:,:),fault_info(:,:), &
finite_slip(:,:,:),Mfin(:,:,:)
real(kind=kreal) :: Rmat(ndim,ndim),slip_vec(ndim,1),normal_vec(ndim,1)
real(kind=kreal) :: slip,strike,dp,rake
logical :: is_dip,is_rake

integer :: f_elmt,n_felmt

real(kind=kreal),dimension(:,:),allocatable :: dshape_hex8
real(kind=kreal),dimension(:),allocatable :: lagrange_gll
real(kind=kreal),dimension(:,:),allocatable :: dlagrange_gll

real(kind=kreal) :: sload(0:lneq)

integer :: i_gll,ielmt_min,inode_min
integer :: nelmt_srctry
integer,allocatable :: ielmt_srctry(:)
logical,allocatable :: isnode(:),iselmt(:)
real(kind=kreal) :: minsqdist,sqdist

integer :: ielmt_src
integer :: indmin(1),src_inrank
real(kind=kreal) :: gminerr,gmaxerr
real(kind=kreal) :: gminerr_src,gmaxerr_src
real(kind=kreal) :: all_minerr(1,0:nproc-1)
integer :: ipass_strict,nfail_strict
logical :: isinside

character(len=1) :: tchar
character(len=250) :: line,tag
character(len=80) :: token
character(len=80) :: fname
character(len=80) :: data_path
character(len=250) :: pfile

errtag="ERROR: unknown!"
errcode=-1

! set data path
data_path=trim(inp_path)

if(eqsource_type==0)then
! slip distribution on the fault
  fname=trim(data_path)//trim(slipfile)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif

  !slip input
  is_fcenter   = .false.
  is_fstrike   = .false.
  is_fnormal   = .false.
  is_fslip     = .false.
  is_flength   = .false.
  is_fwidth    = .false.
  is_fnplength = .false.
  is_fnpwidth  = .false.
  do i_line=1,nmax_line
    read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
    if (ios/=0)exit
    ! check for blank and comment line
    if (isblank(line) .or. iscomment(line,'#'))cycle

    !first token
    tag=trim(line)
    call first_token(tag,token)

    ! fault center
    if (trim(token)=='fault_center')then
      if(is_fcenter)then
        write(errtag,*)'WARNING: copy of fault_center found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_center
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_center!'
        stop
      endif
      fault_center=NONDIM_L*fault_center
      is_fcenter=.true.
    endif
    ! fault strike
    if (trim(token)=='fault_strike')then
      if(is_fstrike)then
        write(errtag,*)'WARNING: copy of fault_strike found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_strike(:,1)
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_strike!'
        stop
      endif
      is_fstrike=.true.
    endif
    ! fault normal
    if (trim(token)=='fault_normal')then
      if(is_fnormal)then
        write(errtag,*)'WARNING: copy of fault_normal found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_normal(:,1)
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_normal!'
        stop
      endif
      is_fnormal=.true.
    endif
    ! fault slip
    if (trim(token)=='fault_slip')then
      if(is_fslip)then
        write(errtag,*)'WARNING: copy of fault_slip found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_slip
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_slip!'
        stop
      endif
      fault_slip=NONDIM_L*fault_slip
      is_fslip=.true.
    endif
    ! fault length
    if (trim(token)=='fault_length')then
      if(is_flength)then
        write(errtag,*)'WARNING: copy of fault_length found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_length
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_length!'
        stop
      endif
      fault_length=NONDIM_L*fault_length
      is_flength=.true.
    endif
    ! fault width
    if (trim(token)=='fault_width')then
      if(is_fwidth)then
        write(errtag,*)'WARNING: copy of fault_width found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_width
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_width!'
        stop
      endif
      fault_width=NONDIM_L*fault_width
      is_fwidth=.true.
    endif
    ! number of patches along the length
    if (trim(token)=='npatch_length')then
      if(is_fnplength)then
        write(errtag,*)'WARNING: copy of npatch_length found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,npatch_length
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read npatch_length!'
        stop
      endif
      is_fnplength=.true.
    endif
    ! number of patches along the width
    if (trim(token)=='npatch_width')then
      if(is_fnpwidth)then
        write(errtag,*)'WARNING: copy of npatch_width found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,npatch_width
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read npatch_width!'
        stop
      endif
      is_fnpwidth=.true.
    endif
  enddo

  ! check slip information reading status
  if(.not.is_fcenter)then
    write(*,*)'ERROR: "fault_center" is missing in slip information file!'
    stop
  endif
  if(.not.is_fstrike)then
    write(*,*)'ERROR: "fault_strike" is missing in slip information file!'
    stop
  endif
  if(.not.is_fnormal)then
    write(*,*)'ERROR: "fault_normal" is missing in slip information file!'
    stop
  endif
  if(.not.is_fslip)then
    write(*,*)'ERROR: "fault_slip" is missing in slip information file!'
    stop
  endif
  if(.not.is_flength)then
    write(*,*)'ERROR: "fault_length" is missing in slip information file!'
    stop
  endif
  if(.not.is_fwidth)then
    write(*,*)'ERROR: "fault_width" is missing in slip information file!'
    stop
  endif
  if(.not.is_fnplength)then
    write(*,*)'ERROR: "napatch_length" is missing in slip information file!'
    stop
  endif
  if(.not.is_fnplength)then
    write(*,*)'ERROR: "napatch_length" is missing in slip information file!'
    stop
  endif

  if(myrank==0)then
    write(logunit,*)'------------------------------------------------------------'
    write(logunit,'(a,3(g0.6,1x))')'fault center :',DIM_L*fault_center
    write(logunit,'(a,3(g0.6,1x))')'fault strike :',fault_strike
    write(logunit,'(a,3(g0.6,1x))')'fault normal :',fault_normal
    write(logunit,'(a,3(g0.6,1x))')'fault slip   :',DIM_L*fault_slip
    write(logunit,'(a,g0.6)')'fault length :',DIM_L*fault_length
    write(logunit,'(a,g0.6)')'fault width  :',DIM_L*fault_width
    write(logunit,'(a,i0)')'npatch length:',npatch_length
    write(logunit,'(a,i0)')'npatch width :',npatch_width
    write(logunit,*)'------------------------------------------------------------'
    flush(logunit)
  endif

  !compute fault area and patches
  ! make unit vector if not already
  fault_strike=fault_strike/norm(fault_strike(:,1))
  fault_normal=fault_normal/norm(fault_normal(:,1))

  !-------parameters for Okada solution-----------------------------------------
  ! only for eqsource_type==0
  if(benchmark_okada)then
    ! Okada origin is directly above the fault center on the free surface
    okada_origin=fault_center
    okada_origin(3)=ZERO
    okada_depth=-fault_center(NDIM) ! Positive down
    okada_aw2=HALF*fault_width
    okada_al2=HALF*fault_length
    okada_aw1=-okada_aw2
    okada_al1=-okada_al2

    ! Dip angle = angle made by a dipping plane with the horizontal
    !           = angle made by a fault normal with the vertical
    ! Note: Positive vertical is always along positive Z direction
    !okada_dip=acos(dot_product(vecOC_norm,(/1,0,0/)))*RAD2DEG
    okada_dip=acos(dot_product(fault_normal(:,1),(/ ZERO,ZERO,ONE /)))*RAD2DEG

    ! vector along the dip
    ! first rotate strike vector clockwise by 90 degree
    fault_strike90=vector_rotateZ(fault_strike(:,1),-90.0d0)
    fault_strike90=fault_strike90/norm(fault_strike90)

    ! calculate projection of strike90 onto fault normal
    proj=dot_product(fault_normal(:,1),fault_strike90)
    ! Note: both fault_noraml and fault_strike90 are unit vectors
    if(proj.eq.ONE)then
      ! Dip is vertical
      fault_dip=(/ ZERO,ZERO,-ONE /)
    elseif(proj.eq.-ONE)then
      write(errtag,*)'ERROR: dot product of fault_normal and fault_strike90 is &
      &-1.0!'
      return
    else
      ! method 1
      fault_dip=fault_strike90- &
              fault_normal(:,1)*dot_product(fault_normal(:,1),fault_strike90)
      ! method 2
      !! vector perpendicular to the plane formed by fault noraml and strike90
      !vecN=cross_product(fault_strike90,fault_normal(:,1))
      !! vector along dip
      !fault_dip=cross_product(fault_normal(:,1),vecN)
    endif
    ! unit vector along dip
    fault_dip=fault_dip/norm(fault_dip)
    !print*,'fault_dip1:',fault_dip

    ! slip component along fault strike
    okada_disl1=dot_product(fault_slip(:,1),fault_strike(:,1))
    ! slip component along fault dip
    okada_disl2=dot_product(fault_slip(:,1),fault_dip)
    !slip componenet along the fault normal
    okada_disl3=dot_product(fault_slip(:,1), &
                            fault_normal(:,1)/norm(fault_normal(:,1)))

    if(myrank==0)then
      write(logunit,*)'NOTE: strike direction must be along X axis for Okada!'
      write(logunit,*)'Okada parameters: slip along strike, dip, and normal'
      write(logunit,*)'okada_disl',okada_disl1,okada_disl2,okada_disl3
      flush(logunit)
    endif
  endif ! benchmark_okada
  !-----------------------------------------------------------------------------
  
! find corner points
! 4-------------C-------------3
! |                           |
! |                           |
! A             O------>s     B
! |                           |
! |                           |
! 1-------------D-------------2

  ! vectors
  vecOB=half*fault_length*fault_strike(:,1)
  vecOA=-vecOB
  vecOC=cross_product(fault_normal(:,1),fault_strike(:,1))
  vecOC_norm=vecOC/norm(vecOC) ! unit vector
  vecOC=half*fault_width*vecOC_norm
  vecOD=-vecOC

  ! four corners
  ! related vectors represent the coordinates of the related corners
  f_cx1=fault_center+(vecOA+vecOD)
  f_cx2=fault_center+(vecOB+vecOD)
  f_cx3=fault_center+(vecOB+vecOC)
  f_cx4=fault_center+(vecOA+vecOC)

  ! rectangular patches
  npatch=npatch_length*npatch_width

  ! number of patch corner points (NOT center)
  np_l=npatch_length+1
  np_w=npatch_width+1
  np=np_l*np_w

  ! connectivity of patches
  allocate(p_connect(4,npatch))
  inum=0
  i0=1 ! starting
  do i_w=1,npatch_width

    do i_l=1,npatch_length
      inum=inum+1
      i1=i0
      i2=i1+1
      i4=i1+np_l
      if(i_l.eq.1)i0_nextrow=i4
      i3=i4+1
      i0=i2
      p_connect(:,inum)=(/ i1,i2,i3,i4 /)
    enddo
    i0=i0_nextrow
  enddo

  rw=one/real(npatch_width,kreal)
  rl=one/real(npatch_length,kreal)
  ! coordinates of patches
  allocate(p_coord(3,np_l*np_w))
  inum=0
  do i_w=1,np_w
    ! find two end points
    if(i_w.eq.1)then
      wx1=f_cx1
      wx2=f_cx2
    elseif(i_w.eq.np_w)then
      wx1=f_cx4
      wx2=f_cx3
    else
      ! division formula
      rk=real(i_w-1,kreal)*rw
      wx1=rk*f_cx4+(one-rk)*f_cx1
      wx2=rk*f_cx3+(one-rk)*f_cx2
    endif

    do i_l=1,np_l
      inum=inum+1
      if(i_l.eq.1)then
        ! starting point
        p_coord(:,inum)=wx1
      elseif(i_l.eq.np_l)then
        ! end point
        p_coord(:,inum)=wx2
      else
        ! division formula
        rk=real(i_l-1,kreal)*rl
        p_coord(:,inum)=rk*wx2+(one-rk)*wx1
      endif
    enddo
  enddo

  ! plot VTK file: patches
  pfile=trim(out_path)//trim(file_head)//'_patch'//trim(ptail)//'.vtk'
  open(100,file=trim(pfile),action='write',status='replace')

  write(100,'(a)')'# vtk DataFile Version 2.0'
  write(100,'(a)')'Unstructured Grid Example'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(100,'(a,1x,i0,1x,a)')'POINTS',np,'float'
  do i_p=1,np
    write(100,'(3(f14.6,1x))')DIM_L*p_coord(:,i_p)
  enddo
  write(100,*)
  write(100,'(a,1x,i0,1x,i0)')'CELLS',npatch,5*npatch
  do i_p=1,npatch
    ! VTK format indexing starts from 0
    write(100,'(5(i0,1x))')4,p_connect(:,i_p)-1
  enddo
  write(100,*)
  write(100,'(a,1x,i0)')'CELL_TYPES',npatch
  do i_p=1,npatch
    write(100,'(i0)')9
  enddo
  close(100)

  ! patch center points
  ! source coordinates
  allocate(source_coord(3,npatch))
  do i_p=1,npatch
    ! mid point of a diagonal is the center point of the rectangle
    i1=p_connect(1,i_p)
    i3=p_connect(3,i_p)
    source_coord(:,i_p)=0.5_kreal*(p_coord(:,i1)+p_coord(:,i3))
  enddo
  ! plot VTK file: source locations
  pfile=trim(out_path)//trim(file_head)//'_sources'//trim(ptail)//'.vtk'
  open(100,file=trim(pfile),action='write',status='replace')

  write(100,'(a)')'# vtk DataFile Version 2.0'
  write(100,'(a)')'Unstructured Grid Example'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(100,'(a,i5,a)')'POINTS',npatch,' float'
  do i_p=1,npatch
    write(100,'(3(f14.6,1x))')DIM_L*source_coord(:,i_p)
  enddo
  write(100,*)
  write(100,'(a,i5,a,i5)')'CELLS',npatch,' ',2*npatch
  do i_p=1,npatch
    ! VTK format indexing starts from 0
    write(100,'(5(i5,1x))')1,i_p-1
  enddo
  write(100,*)
  write(100,'(a,i5)')'CELL_TYPES',npatch
  do i_p=1,npatch
    write(100,'(i1)')1
  enddo
  close(100)

  pfac=one/real(npatch,kreal) ! patch factor

  fault_area=fault_length*fault_width

  ! Compute moment-density tensor: only the material independent portion
  ! it will be later multiplied by appropriate shear modulus
  ! fault_slip=slip*slip_vector
  ! The general relation for moment-density tensor
  ! m=C: (n s)
  ! The matrix (n s) is not symmetric!
  ! See Hom Nath Gharti's note on 'fault' or Dahlen & Tromp PP 157 Eq 5.44
  ! NEED to check why the following doesn't work.
  Mmat=fault_area*matmul(fault_normal,transpose(fault_slip))

  ! store in 1D array
  !M0(1,1)=Mmat(1,1)
  !M0(1,2)=Mmat(1,2)
  !M0(1,3)=Mmat(1,3)
  !M0(1,4)=Mmat(2,1)
  !M0(1,5)=Mmat(2,2)
  !M0(1,6)=Mmat(2,3)
  !M0(1,7)=Mmat(3,1)
  !M0(1,8)=Mmat(3,2)
  !M0(1,9)=Mmat(3,3)
  ! store in Voigt notation order
  !xx,yy,zz,xy,yz,zx
  !11,22,33,12,23,31
  M0_V(1,1)=Mmat(1,1)
  M0_V(2,1)=Mmat(2,2)
  M0_V(3,1)=Mmat(3,3)
  M0_V(4,1)=Mmat(1,2)
  M0_V(5,1)=Mmat(2,3)
  M0_V(6,1)=Mmat(3,1)
  M0_V(7,1)=Mmat(2,1)
  M0_V(8,1)=Mmat(3,2)
  M0_V(9,1)=Mmat(1,3)

  nsource = npatch

  close(11)

elseif(eqsource_type==1)then
! CMT solution
  ! read CMTSOLUTION
  ! This CMTSOLUTION file is not the actual CMTSOLUTION. Moment tensor is stored
  ! in a spherical coordinates in the actual CMTSOLUTION. In order to preserve
  ! the format we use CMTSOLUTION convention according to following assumptions:
  !
  ! latitude is y, longitude is x, depth is z (all three in meters)
  ! r -> z, theta -> -y, phi -> x
  ! moment tensor components are in CGS units
  ! NOTE: we could add an option to read CMTSOLUTION in a global format as well!
  !
  !  Mrr =  Mzz
  !  Mtt =  Myy
  !  Mpp =  Mxx
  !  Mrt = -Myz
  !  Mrp =  Mxz
  !  Mtp = -Mxy

  ! count CMT sources
  call count_cmtsolution(errcode,errtag)
  if(errcode.ne.0)return
  if(myrank==0)then
    write(logunit,'(a,i0)')'Total number of CMT sources: ',ncmt_source
    flush(logunit)
  endif
  allocate(source_coord(3,ncmt_source),M_cmt(6,ncmt_source))
  ! Store 6 unique moment-tensor components in the order
  ! 1: Mxx
  ! 2: Myy
  ! 3: Mzz
  ! 4: Mxy
  ! 5: Myz
  ! 6: Mzx
  ! read CMT sources
  call read_cmtsolution(source_coord,M_cmt,errcode,errtag)
  if(errcode.ne.0)return
  ! NOTE: all quanitites are nondimensionalized within read_cmtsolution routine

  ! plot VTK file: CMT sources
  pfile=trim(out_path)//trim(file_head)//'_cmt_sources'//trim(ptail)//'.vtk'
  open(100,file=trim(pfile),action='write',status='replace')

  write(100,'(a)')'# vtk DataFile Version 2.0'
  write(100,'(a)')'Unstructured Grid Example'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(100,'(a,1x,i0,1x,a)')'POINTS',ncmt_source,' float'
  do i_p=1,ncmt_source
    write(100,'(3(f14.6,1x))')DIM_L*source_coord(:,i_p)
  enddo
  write(100,*)
  write(100,'(a,1x,i0,1x,i0)')'CELLS',ncmt_source,2*ncmt_source
  do i_p=0,ncmt_source-1 ! VTK format indexing starts from 0
    write(100,'(5(i5,1x))')1,i_p
  enddo
  write(100,*)
  write(100,'(a,i0)')'CELL_TYPES ',ncmt_source
  do i_p=1,ncmt_source
    write(100,'(i0)')1
  enddo
  close(100)

  nsource=ncmt_source

elseif(eqsource_type==11)then
! Customized CMT solution format
  ! read CMTSOLUTION
  ! This CMTSOLUTION file is not the actual CMTSOLUTION. Moment tensor is stored
  ! in a spherical coordinates in the actual CMTSOLUTION. In order to preserve
  ! the format we use CMTSOLUTION convention according to following assumptions:
  !
  ! latitude is y, longitude is x, depth is z (all three in meters)
  ! r -> z, theta -> -y, phi -> x
  ! moment tensor components are in CGS units
  ! NOTE: we could add an option to read CMTSOLUTION in a global format as well!
  !
  !  Mrr =  Mzz
  !  Mtt =  Myy
  !  Mpp =  Mxx
  !  Mrt = -Myz
  !  Mrp =  Mxz
  !  Mtp = -Mxy

  ! need to add something here to convert to cartesian system
  ! above relation should be sufficient for cartesian code!

  ! convert to CGS to SI units
  fname=trim(data_path)//trim(cmtfile)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif

!  is_flat= .false.
!  is_flon = .false.
  is_fcenter   = .false.
  is_mzz       = .false.
  is_myy = .false.
  is_mxx = .false.
  is_myz = .false.
  is_mxz = .false.
  is_mxy = .false.
  do i_line=1,nmax_line
    read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
    if (ios/=0)exit
    ! check for blank and comment line
    if (isblank(line) .or. iscomment(line,'#'))cycle

    !first token
    tag=trim(line)
    call first_token(tag,token)

    ! fault center
    if (trim(token)=='fault_center')then
      if(is_fcenter)then
        write(errtag,*)'WARNING: copy of fault_center found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_center
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_center!'
        stop
      endif
      is_fcenter=.true.
    endif
! fault latitude
!    if (trim(token)=='latitude')then
!      if(is_flat)then
!        write(errtag,*)'WARNING: copy of latitude found! Copies are discarded!'
!        cycle
!      endif
!      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar, cmt_lat
!      if(istat.ne.0)then
!        write(*,*)'ERROR: cannot read latitude!'
!        stop
!      endif
!      is_flat=.true.
!    endif
! fault longitude
!    if (trim(token)=='longitude')then
!      if(is_flon)then
!        write(errtag,*)'WARNING: copy of longitude found! Copies are discarded!'
!        cycle
!      endif
!      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar, cmt_lon
!      if(istat.ne.0)then
!        write(*,*)'ERROR: cannot read longitude!'
!        stop
!      endif
!      is_flon=.true.
!    endif
    ! moment tensor parameters
    !mzz
    if (trim(token)=='mzz')then
      if(is_mzz)then
        write(errtag,*)'WARNING: copy of mzz found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,mzz
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read mzz!'
        stop
      endif
      is_mzz=.true.
    endif
    !myy
    if (trim(token)=='myy')then
      if(is_myy)then
        write(errtag,*)'WARNING: copy of myy found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,myy
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read myy!'
        stop
      endif
      is_myy=.true.
    endif
    !mxx
    if (trim(token)=='mxx')then
      if(is_mxx)then
        write(errtag,*)'WARNING: copy of mxx found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,mxx
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read mxx!'
        stop
      endif
      is_mxx=.true.
    endif
    !myz
    if (trim(token)=='myz')then
      if(is_myz)then
        write(errtag,*)'WARNING: copy of myz found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,myz
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read myz!'
        stop
      endif
      is_myz=.true.
    endif
    !mxz
    if (trim(token)=='mxz')then
      if(is_mxz)then
        write(errtag,*)'WARNING: copy of mxz found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,mxz
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read mxz!'
        stop
      endif
      is_mxz=.true.
    endif
    !mxy
    if (trim(token)=='mxy')then
      if(is_mxy)then
        write(errtag,*)'WARNING: copy of mxy found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,mxy
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read mxy!'
        stop
      endif
      is_mxy=.true.
    endif
  enddo

  ! check CMT information reading status
  if(.not.is_fcenter)then
    write(*,*)'ERROR: "fault_center" is missing in slip information file!'
    stop
  endif
  if(.not.is_mzz)then
    write(*,*)'ERROR: "mzz" is missing in slip information file!'
    stop
  endif
  if(.not.is_myy)then
    write(*,*)'ERROR: "myy" is missing in slip information file!'
    stop
  endif
  if(.not.is_mxx)then
    write(*,*)'ERROR: "mxx" is missing in slip information file!'
    stop
  endif
  if(.not.is_myz)then
    write(*,*)'ERROR: "myz" is missing in slip information file!'
    stop
  endif
  if(.not.is_mxz)then
    write(*,*)'ERROR: "mxz" is missing in slip information file!'
    stop
  endif
  if(.not.is_mxy)then
    write(*,*)'ERROR: "mxy" is missing in slip information file!'
    stop
  endif

  if(myrank==0)then
    write(logunit,*)'------------------------------------------------------------'
    write(logunit,'(a,3(e14.6,1x))')'fault center :',fault_center
 !   write(*,'(a,i6)')'mzz,myy,mxx,myz,mxz,mxy :'mzz,myy,mxx,myz,mxz,mxy
    write(logunit,*)'------------------------------------------------------------'
    flush(logunit)
  endif

!  cmt_M(1,1)=mrr
!  cmt_M(1,2)=mrt
!  cmt_M(1,3)=mrp
!  cmt_M(2,1)=mrt
!  cmt_M(2,2)=mtt
!  cmt_M(2,3)=mtp
!  cmt_M(3,1)=mrp
!  cmt_M(3,2)=mtp
!  cmt_M(3,3)=mpp

      ! convert geographic latitude lat (degrees) to geocentric colatitude theta (radians)
!      call lat_2_geocentric_colat_dble(cmt_lat,theta)                       
      
      ! convert latitude to radians                                                                             
!      phi = cmt_lon*DEG2RAD                                     
!      call reduce(theta,phi)

      ! compute spherical to cartesian transformation matrix
!      call compute_Tmat_spher2cart(theta,phi,Tmat)

      ! convert spherical moment tensor to cartesian
!      xyz_M=matmul(Tmat,matmul(cmt_M,transpose(Tmat)))

      !compute rotation matrix from global to local xyz coords
!      call compute_local_rot_mat(theta,phi,Lmat)
!      local_M=matmul(Lmat,matmul(xyz_M,transpose(Lmat)))

  local_M(1,1)=mxx
  local_M(1,2)=mxy
  local_M(1,3)=mxz
  local_M(2,1)=mxy
  local_M(2,2)=myy
  local_M(2,3)=myz
  local_M(3,1)=mxz
  local_M(3,2)=myz
  local_M(3,3)=mzz
 
  M(1,1)=local_M(1,1)
  M(1,2)=local_M(1,2)
  M(1,3)=local_M(1,3)
  M(1,4)=local_M(2,1)
  M(1,5)=local_M(2,2)
  M(1,6)=local_M(2,3)
  M(1,7)=local_M(3,1)
  M(1,8)=local_M(3,2)
  M(1,9)=local_M(3,3)

  nsource=1
  allocate(source_coord(3,1))
  source_coord(:,1)=fault_center
  pfac=one

  close(11)

  ! dimensionlessalize
  M=M*NONDIM_MTENS
  source_coord=source_coord*NONDIM_L 

elseif(eqsource_type==2)then
! finite fault with slip information
  ! read finite fault metadata (npatches, faultnormal) from faultmetafile
  fname=trim(data_path)//trim(faultmetafile)
  open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif

  is_fstrike = .false.
  is_dip = .false.
  is_rake = .false.
  is_fnpatch  = .false.
  is_flength = .false.
  is_fwidth = .false.

  do i_line=1,nmax_line
    read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
    if (ios/=0)exit
    ! check for blank and comment line
    if (isblank(line) .or. iscomment(line,'#'))cycle
    !first token
    tag=trim(line)
    call first_token(tag,token)
    ! strike
    if (trim(token)=='strike')then
      if(is_fstrike)then
        write(errtag,*)'WARNING: copy of strike found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,strike
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read strike!'
        stop
      endif
      is_fstrike=.true.
    endif
    ! dip
    if (trim(token)=='dip')then
      if(is_dip)then
        write(errtag,*)'WARNING: copy of dip found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,dp
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read dip!'
        stop
      endif
      is_dip=.true.
    endif
    ! rake
    if (trim(token)=='rake')then
      if(is_rake)then
        write(errtag,*)'WARNING: copy of rake found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,rake
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read rake!'
        stop
      endif
      is_rake=.true.
    endif
    ! fault npatch
    if (trim(token)=='npatch')then
      if(is_fnpatch)then
        write(errtag,*)'WARNING: copy of npatch found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,npatch
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read npatch!'
        stop
      endif
      is_fnpatch=.true.
    endif
    ! fault length
    if (trim(token)=='length')then
      if(is_flength)then
        write(errtag,*)'WARNING: copy of fault_length found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_length
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_length!'
        stop
      endif
      is_flength=.true.
    endif
    ! fault width
    if (trim(token)=='width')then
      if(is_fwidth)then
        write(errtag,*)'WARNING: copy of fault_width found! Copies are discarded!'
        cycle
      endif
      read(line(len_trim(token)+1:len_trim(line)),*,iostat=istat)tchar,fault_width
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read fault_width!'
        stop
      endif
      is_fwidth=.true.
    endif
  enddo ! do i_line=1,nmax_line

  ! check fault information reading status
  if(.not.is_fstrike)then
    write(*,*)'ERROR: "fault_strike" is missing in fault information file!'
    stop
  endif
  if(.not.is_dip)then
    write(*,*)'ERROR: "fault_dip" is missing in fault information file!'
    stop
  endif
  if(.not.is_rake)then
    write(*,*)'ERROR: "fault_rake" is missing in fault information file!'
    stop
  endif
  if(.not.is_fnpatch)then
    write(*,*)'ERROR: "fault_npatch" is missing in fault information file!'
    stop
  endif
  if(.not.is_flength)then
    write(*,*)'ERROR: "fault_length" is missing in slip information file!'
    stop
  endif
  if(.not.is_fwidth)then
    write(*,*)'ERROR: "fault_width" is missing in slip information file!'
    stop
  endif
 
  if(myrank==0)then
    write(logunit,*)'------------------------------------------------------------'
    write(logunit,'(a,e14.6)')'fault length :',fault_length
    write(logunit,'(a,e14.6)')'fault width  :',fault_width
    write(logunit,'(a,i6)')'npatch:',npatch
    write(logunit,'(a,e14.6)')'strike:',strike
    write(logunit,'(a,e14.6)')'dip:',dp
    write(logunit,'(a,e14.6)')'rake:',rake
    write(logunit,*)'------------------------------------------------------------'
    flush(logunit)
  endif

  ! read finite fault data from gridded file. columns are x,y,z coords and slip
  ! value
  fname=trim(data_path)//trim(faultfile)
  open(unit=12,file=trim(fname),status='old',action='read',iostat = ios)
  if( ios /= 0 ) then
    write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
    return
  endif

  fault_normal(1,1)=-sin(dp)*sin(strike)
  fault_normal(2,1)=-sin(dp)*cos(strike)
  fault_normal(3,1)=cos(dp)

  fault_slip(1,1)=cos(rake)*cos(strike) + sin(rake)*cos(dp)*sin(strike)
  fault_slip(2,1)=-cos(rake)*sin(strike)+sin(rake)*cos(dp)*cos(strike)
  fault_slip(3,1)=sin(rake)*sin(dp)

  !rotate fault normal and slip vectors to our coordinate system
  Rmat(1,1)=cos(strike-90.0)
  Rmat(1,2)=sin(strike-90.0)
  Rmat(3,3)=0.0
  Rmat(2,1)=-sin(strike-90.0)
  Rmat(2,2)=cos(strike-90.0)
  Rmat(2,3)=0.0
  Rmat(3,1)=0.0
  Rmat(3,2)=0.0
  Rmat(3,3)=1.0

  slip_vec=matmul(Rmat,fault_slip)
  slip_vec=slip_vec/norm(slip_vec(:,1))
  normal_vec=matmul(Rmat,fault_normal)
  normal_vec=normal_vec/norm(normal_vec(:,1))

  pfac=one/real(npatch,kreal) ! patch factor
  fault_area=fault_length*fault_width

  allocate(fault_info(4,npatch))
  allocate(finite_slip(ndim,1,npatch))
  allocate(Mfin(ndim,ndim,npatch))
  nsource=npatch
  read(12,*,iostat=ios)fault_info
  allocate(source_coord(3,npatch))

  do ix = 1,nsource
    source_coord(:,ix)=fault_info(1:3,ix)
    slip=fault_info(4,ix)
    finite_slip(:,:,ix)= slip*slip_vec
  enddo    

  ! plot VTK file: source locations
  pfile=trim(out_path)//trim(file_head)//'_sources'//trim(ptail)//'.vtk'
  open(100,file=trim(pfile),action='write',status='replace')

  write(100,'(a)')'# vtk DataFile Version 2.0'
  write(100,'(a)')'Unstructured Grid Example'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(100,'(a,i5,a)')'POINTS',npatch,' float'
  do i_p=1,npatch
    write(100,'(3(f14.6,1x))')source_coord(:,i_p)
  enddo
  write(100,*)
  write(100,'(a,i5,i5)')'CELLS',npatch,2*npatch
  do i_p=1,npatch
    ! VTK format indexing starts from 0
    write(100,'(5(i5,1x))')1,i_p-1
  enddo
  write(100,*)
  write(100,'(a,i5)')'CELL_TYPES',npatch
  do i_p=1,npatch
    write(100,'(i1)')1
  enddo
  close(100)

  ! dimensionlessalize
  finite_slip=finite_slip*NONDIM_L
  fault_area=fault_area*NONDIM_L*NONDIM_L
  source_coord=source_coord*NONDIM_L

  do ix = 1,nsource
    ! compute moment tensor: only the material independent portion
    ! it will be later multiplied by appropriate shear modulus
    ! fault_slip=slip*slip_vector
    Mfin(:,:,ix)=fault_area*(matmul(finite_slip(:,:,ix),transpose(fault_normal))+    &
       matmul(fault_normal,transpose(finite_slip(:,:,ix))))
    ! See the comment above for general relation
    !Mfin(:,:,ix)=fault_area*matmul(fault_normal,transpose(finite_slip(:,:,ix)))

  enddo
  close(11)
  close(12)
else
  write(errtag,'(a,i0)')'ERROR:unrecognized eqsource type: ',eqsource_type
  return
endif

!-------------------------------------------------------------------------------
!! TEST for IsPointInHexahedron function
!if(myrank==0)then
!coord(:,1)=(/zero,zero,zero/)
!coord(:,2)=(/one,zero,zero/)
!coord(:,3)=(/one,one,zero/)
!coord(:,4)=(/zero,one,zero/)
!coord(:,5)=(/zero,zero,one/)
!coord(:,6)=(/one,zero,one/)
!coord(:,7)=(/one,one,one/)
!coord(:,8)=(/zero,one,one/)
!xp=(/one,one,one/)
!call IsPointInHexahedron(coord,xp,isinside)
!print*,'TEST:',xp,isinside
!xp=(/half,half,half/)
!call IsPointInHexahedron(coord,xp,isinside)
!print*,'TEST:',xp,isinside
!xp=(/1.02_kreal,one,one/)
!call IsPointInHexahedron(coord,xp,isinside)
!print*,'TEST:',xp,isinside
!xp=(/half,0.05_kreal,0.001_kreal/)
!call IsPointInHexahedron(coord,xp,isinside)
!print*,'TEST:',xp,isinside
!endif
!call sync_process
!stop
!-------------------------------------------------------------------------------

! Get derivatives of shape functions for 8-noded hex.
allocate(dshape_hex8(ngnode,ndim))
allocate(lagrange_gll(ngll),dlagrange_gll(ndim,ngll))

imid=(ngll+1)/2

allocate(isrc_located(nsource))
isrc_located=0
! Compute contribution of the earthquake source/s to the elemental loads.
allocate(isnode(nnode),iselmt(nelmt))
gminerr_src=INFTOL
gmaxerr_src=ZERO
nfail_strict=0
source: do i_src=1,nsource
  sload=zero
  source_x=source_coord(:,i_src)
  ! Find the element which contains this earthquake source.
  is_located=.false.
  
  !print*,myrank, source_x(1),pmodel_minx,pmodel_maxx
  !print*,myrank, source_x(2),pmodel_miny,pmodel_maxy
  !print*,myrank, source_x(3),pmodel_minz,pmodel_maxz

  ! Check if the source is within the model range.
  prange:if(source_x(1).lt.NONDIM_L*pmodel_minx .or. source_x(1).gt.NONDIM_L*pmodel_maxx .or. & 
     source_x(2).lt.NONDIM_L*pmodel_miny .or. source_x(2).gt.NONDIM_L*pmodel_maxy .or. & 
     source_x(3).lt.NONDIM_L*pmodel_minz .or. source_x(3).gt.NONDIM_L*pmodel_maxz)then
    nelmt_srctry=0
    ! NOTE: we cannot cycle the loop here otherwise the process hangs for ever
    ! because inbetween there are some statements which require all procesors!
  else
 
    ! First, find the element which is not too far and has the nearest GLL point
    ! from the target point.
    ielmt_min=1
    inode_min=g_num(1,ielmt_min)
    minsqdist=INFTOL
    do i_elmt=1,nelmt
      ! Do not locate sources in transition or infinite elements.
      mdomain=mat_domain(mat_id(i_elmt))
      if(mdomain==ELASTIC_TRINFDOMAIN .or.      &
         mdomain==ELASTIC_INFDOMAIN   .or.      &
         mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
         mdomain==VISCOELASTIC_INFDOMAIN)cycle
      ! Preliminary test. Discard far enough elements.
      num=g_num(:,i_elmt)
      xp=g_coord(:,num(imid))
      if(sqdistance(source_x,xp,ndim).gt.sqmaxsize_elmt)cycle

      do i_gll=1,ngll
        xp=g_coord(:,num(i_gll))
        sqdist=sqdistance(source_x,xp,ndim)
        if(sqdist.lt.minsqdist)then
          minsqdist=sqdist
          ielmt_min=i_elmt
          inode_min=num(i_gll)
        endif
      enddo
    enddo

    ! Count the number of and tag elements that share the GLL point just found.
    isnode=.false.
    isnode(g_num(hex8_gnode,ielmt_min))=.true.
    iselmt=.false.
    iselmt(ielmt_min)=.true.
    nelmt_srctry=1
    do i_elmt=1,nelmt
      if(i_elmt.eq.ielmt_min)cycle
      if(any(isnode(g_num(hex8_gnode,i_elmt))))then
        nelmt_srctry=nelmt_srctry+1
        iselmt(i_elmt)=.true.
      endif
    enddo
    ! Find and assign the list of elements to try.
    allocate(ielmt_srctry(nelmt_srctry))
    ielmt_srctry=-9999
    inum=0
    do i_elmt=1,nelmt
      if(iselmt(i_elmt))then
        inum=inum+1
        ielmt_srctry(inum)=i_elmt
      endif
    enddo
    if(inum.ne.nelmt_srctry)then
      write(*,*)'ERROR: nelmt_srctry mismatch!'
      stop
    endif
  endif prange
  !print*,'try elements:',i_src,myrank,nelmt_srctry,count(iselmt) 

  n_felmt=0
  minerr=INFTOL
  ! Find the actual element that contains the source.
  ! Loop through the list of elements
  ! STRICT TEST: Element must be a proper hexagon, i.e., must have six faces. 
  ipass_strict=0
  element_try1: do i_elmt=1,nelmt_srctry
    ielmt=ielmt_srctry(i_elmt)

    num=g_num(:,ielmt)
    coord=g_coord(:,num(hex8_gnode))
    ! Check if the source is in this element.
    call IsPointInHexahedron(coord,source_x,isinside)
    if(isinside)then
      n_felmt=n_felmt+1
      ! Map source location to natural coordinates.
      call  map_point2naturalhex8(coord,source_x,xip,located_x,niter,errx)
      errxd=errx*DIM_L
      if(n_felmt==1)then
        ! Initialize errors
        minerr=errxd
        ! Set source element 
        ielmt_src=ielmt
      elseif(n_felmt.gt.1)then
        if(errxd.lt.minerr)then
          minerr=errxd 
          ! Reset source element 
          ielmt_src=ielmt
        endif
      endif
      ipass_strict=1
    endif !(isinside)
  enddo element_try1 ! i_elmt
  
  ipass_strict=sumscal(ipass_strict)
  ! FLEXIBLE TEST: Element may be a proper hexagon, i.e., may have more than 
  !                six faces. 
  ! There will be a WARNING in the log file. 
  ! This test will be performed only if the STRICT TEST fails.
  if(ipass_strict.lt.1)then
  ! Source point was not located on the element list
    nfail_strict=nfail_strict+1
    element_try2: do i_elmt=1,nelmt_srctry
      ielmt=ielmt_srctry(i_elmt)

      num=g_num(:,ielmt)
      coord=g_coord(:,num(hex8_gnode))
      n_felmt=n_felmt+1
      ! map source location to natural coordinates
      call  map_point2naturalhex8(coord,source_x,xip,located_x,niter,errx)
      errxd=errx*DIM_L
      if(n_felmt==1)then
        ! initialize errors
        minerr=errxd
        ! Set source element 
        ielmt_src=ielmt
      elseif(n_felmt.gt.1)then
        if(errxd.lt.minerr)then
          minerr=errxd 
          ! Reset source element 
          ielmt_src=ielmt
        endif
      endif
    enddo element_try2 ! i_elmt
  endif !(ipass_strict.lt.1)
  if(allocated(ielmt_srctry))deallocate(ielmt_srctry)
  ! Only the element with the least error was chosen
  n_felmt=sumscal(n_felmt)
  if(n_felmt.lt.1)then
    if(myrank==0)then
      write(logunit,'(a,3(g0.6,1x))')'WARNING: source cannot be located: ',i_src,DIM_L*source_x
      flush(logunit)
    endif
    cycle source
  else 
    isrc_located(i_src)=1
  endif
  gminerr=minscal(minerr)

  if(gminerr.lt.gminerr_src)gminerr_src=gminerr
  if(gminerr.gt.gmaxerr_src)gmaxerr_src=gminerr

  if(myrank==0)then
    if(gminerr.gt.DIM_L*maxsize_elmt)then
      write(logunit,'(a)')'WARNING: this source location may be inaccurate!'
      write(logunit,'(3(g0.6,1x))')xip
      write(logunit,'(2(g0.6,1x))')gminerr,maxsize_elmt
      flush(logunit)
    endif
  endif

  ! Determine the processor this source belongs to.
  all_minerr=allgather_scal(minerr)
  !if(myrank==0)write(*,'(g0)')all_minerr
  indmin=minloc(all_minerr,2)
  ! NOTE: minloc gives the position starting from 1
  src_inrank=indmin(1)-1
  !print*,i_src,src_inrank
  ! Compute source contribution in the appropriate rank
  if(myrank==src_inrank)then
    !print*,'myrank:',myrank,i_src,minerr,gminerr
    f_elmt=ielmt_src
    source_xi=xip
    is_located=.true.
    isrc_located(i_src)=1
    ! compute load
    ! get deriative of shape function at source location
    call dshape_function_hex8p(ngnode,source_xi(1),source_xi(2), &
                               source_xi(3),dshape_hex8)

    ! compute jacobian
    !WAARNING: need to double check jacobian
    jac=matmul(transpose(dshape_hex8),transpose(coord))
    !jac=matmul(coord,dshape_hex8)
    detjac=determinant(jac)
    call invert(jac)

    ! compute GLL lagrange function and their derivative at source location
    call gll_lagrange3d_point(ndim,ngllx,nglly,ngllz,ngll,source_xi(1), &
                              source_xi(2),source_xi(3),       &
                              lagrange_gll,dlagrange_gll)

    deriv=matmul(jac,dlagrange_gll) ! use der for gll
    !\nabla w
    derivmat=zero
    derivmat(1,1:nedofu:3)=deriv(1,:) !dNx/dx
    derivmat(2,2:nedofu:3)=deriv(1,:) !dNy/dx
    derivmat(3,3:nedofu:3)=deriv(1,:) !dNz/dx
    derivmat(4,1:nedofu:3)=deriv(2,:) !dNx/dy
    derivmat(5,2:nedofu:3)=deriv(2,:) !dNy/dy
    derivmat(6,3:nedofu:3)=deriv(2,:) !dNz/dy
    derivmat(7,1:nedofu:3)=deriv(3,:) !dNx/dz
    derivmat(8,2:nedofu:3)=deriv(3,:) !dNy/dz
    derivmat(9,3:nedofu:3)=deriv(3,:) !dNz/dz

    ! it may be a better idea to define material properties for the fault
    ! separately from the model properpties
    ! however, in the statement below, model properties are interpolated on
    ! the fault surface
    ! interpolate material properties at a source point
    bulkmod=dot_product(bulkmod_elmt(:,f_elmt),lagrange_gll)
    shearmod=dot_product(shearmod_elmt(:,f_elmt),lagrange_gll)
    call compute_cmat_elastic(bulkmod,shearmod,cmat)

    egdofu=gdof_elmt(edofu,f_elmt)

    ! Set M0_V: M0 in Voigt notation
    ! if(eqsource_type.eq.0)then
    !    Moment-density tensor patches from the slip information. 
    !    M0_V is already set and constant for all sources.
    ! if(eqsource_type.eq.1.or.eqsource_type.eq.11)then
    !   This is a CMT solution, therfore the M is directly known. We do not
    !   need to compute M0
    if(eqsource_type.eq.2)then
      !pull out M0 for this source element in Voigt notation
      Mmat = Mfin(:,:,i_src)
      M0_V(1,1)=Mmat(1,1)
      M0_V(2,1)=Mmat(2,2)
      M0_V(3,1)=Mmat(3,3)
      M0_V(4,1)=Mmat(1,2)
      M0_V(5,1)=Mmat(2,3)
      M0_V(6,1)=Mmat(3,1)
      M0_V(7,1)=Mmat(2,1)
      M0_V(8,1)=Mmat(3,2)
      M0_V(9,1)=Mmat(1,3)
    endif
    !-------------------------------------------------------------------------

    ! Compute/set full moemnt-tensor M
    if(eqsource_type.eq.0 .or. eqsource_type.eq.2)then
      ! See the comment above for general relation
      ! For C_ijkl n_k s_l or C:(n s) where (n s) is unsymmertic and hence has 
      ! 9 components, we need to modify cmat
      cmatU(1:nst,1:nst)=cmat
      cmatU(:,7)=cmat(:,4)
      cmatU(:,8)=cmat(:,5)
      cmatU(:,9)=cmat(:,6)
      M_V=pfac*matmul(cmatU,M0_V)
      
      ! unpack to 1 X 9 array equivalent to 3 X 3 array
      M(1,1)=M_V(1,1) !Mmat(1,1)
      M(1,2)=M_V(4,1) !Mmat(1,2)
      M(1,3)=M_V(6,1) !Mmat(1,3)
      M(1,4)=M_V(4,1) !Mmat(2,1)
      M(1,5)=M_V(2,1) !Mmat(2,2)
      M(1,6)=M_V(5,1) !Mmat(2,3)
      M(1,7)=M_V(6,1) !Mmat(3,1)
      M(1,8)=M_V(5,1) !Mmat(3,2)
      M(1,9)=M_V(3,1) !Mmat(3,3)
    elseif(eqsource_type.eq.1.or.eqsource_type.eq.11)then
      ! M is already know for this case
      M(1,1:6)=M_cmt(:,i_src)
      M(1,7)=M_cmt(6,i_src) !Mmat(3,1)
      M(1,8)=M_cmt(5,i_src) !Mmat(3,2)
      M(1,9)=M_cmt(3,i_src) !Mmat(3,3)
    else
      write(errtag,'(a,i0)')'ERROR:unrecognized eqsource type: ',eqsource_type
      return
    endif
    ! compute M:\nabla w
    fload=matmul(M,derivmat)
    sload(egdofu)=sload(egdofu)+fload(1,:)
    !exit element ! this will place the source in only one element
  endif !(myrank==src_inrank)
  
  call sync_process
  ! NOTE: We allow only a single element to have the source point. 
  ! Therefore, no averaging is necessary.
  ! However, this may have to be modified for the instance where
  ! the source point lies in the interface of two different material domains!
  !n_felmt=sumscal(n_felmt)
  !if(n_felmt.gt.1)then
  !  ! add average load per source
  !  load=load+sload/real(n_felmt,kreal)
  !else
    load=load+sload
  !endif
enddo source ! i_src
deallocate(isnode,iselmt)

call sync_process
isrc_located=maxvec(isrc_located,nsource)
where(isrc_located.gt.1)isrc_located=1
total_src_located=sum(isrc_located)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0)')'Total defined sources: ',nsource, &
                           'Total located sources: ',total_src_located
  write(logunit,'(a,i0)')'Total sources fail strict test: ',nfail_strict
  if(nfail_strict.gt.0)then
    write(logunit,'(a,i0)')'NOTE: failed strict test indicates that some of &
                        the elements may NOT be proper hexahedra!'
  endif
  flush(logunit)
endif
if(total_src_located.gt.0)then
  if(myrank==0)then
    write(logunit,'(a,g0.6,a,g0.6)')'source location errors (m) => min: ', &
    gminerr_src, ' max: ',gmaxerr_src
    flush(logunit)
  endif
endif
deallocate(isrc_located)
deallocate(dshape_hex8)
deallocate(lagrange_gll,dlagrange_gll)

if(eqsource_type.eq.0)then
  deallocate(p_connect,p_coord)
endif

deallocate(source_coord)
if(eqsource_type==2)then
  deallocate(fault_info)
  deallocate(Mfin)
  deallocate(finite_slip)
endif

errcode=0

return

end subroutine earthquake_load
!===============================================================================

end module earthquake
!===============================================================================
