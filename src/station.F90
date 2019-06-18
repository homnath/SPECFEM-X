module station
use set_precision
implicit none
integer,allocatable :: station_id(:)
integer,allocatable :: station_elmt(:)
real(kind=kreal),allocatable :: station_x(:,:)
real(kind=kreal),allocatable :: station_xi(:,:)

contains
!-------------------------------------------------------------------------------

! DESCRIPTION
!  This routine locates the stations and store their coordinates in a natural
!  coordinates system.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
!  Leah Langer, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO
subroutine locate_station(lneq,errcode,errtag)
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
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

! maximum number of lines in the station file
integer,parameter :: nmax_line=100
integer :: i_elmt,i_line,nline,ielmt,imid,inum,ios,istat,istation
integer :: num(nenode),egdofu(nedofu)

real(kind=kreal) :: coord(ndim,8),jac(ndim,ndim),xp(ndim),xip(ndim)
real(kind=kreal) :: detjac

integer :: i_w,i_l,i_p,i_rec,i0,i1,i2,i3,i4,i0_nextrow
real(kind=kreal) :: located_x(ndim),recx(ndim),recxi(ndim)
real(kind=kreal) :: deriv(ndim,ngll)

integer :: mdomain

logical :: is_located
integer :: this_rec_located,total_rec_located
integer,allocatable:: irec_located(:)
integer :: niter
real(kind=kreal) :: errx,errxd,minerr,maxerr,minerrg,maxerrg

integer :: r_elmt,n_relmt

real(kind=kreal),dimension(:,:),allocatable :: dshape_hex8
real(kind=kreal),dimension(:),allocatable :: lagrange_gll
real(kind=kreal),dimension(:,:),allocatable :: dlagrange_gll

integer :: i_gll,ielmt_min,inode_min
integer :: nelmt_rectry
integer,allocatable :: ielmt_rectry(:)
logical,allocatable :: isnode(:),iselmt(:)
real(kind=kreal) :: minsqdist,sqdist

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

! station files
fname=trim(data_path)//trim(stationfile)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif

nstation=0
nline=0
! count the number of stations
do i_line=1,nmax_line
  ! This will read a line and proceed to next line
  read(11,'(a)',iostat=ios)line
  if (ios/=0)exit
  nline=nline+1
  ! check for blank and comment line
  if (isblank(line) .or. iscomment(line,'#'))cycle
  nstation=nstation+1
enddo
if(myrank==0)then
  write(logunit,'(a,i0)')'total number stations: ',nstation
  flush(logunit)
endif
allocate(station_id(nstation),station_x(NDIM,nstation))

istation=0
do i_line=1,nline
  ! This will read a line and proceed to next line
  read(11,'(a)',iostat=ios)line
  if (ios/=0)exit
  ! check for blank and comment line
  if (isblank(line) .or. iscomment(line,'#'))cycle
  istation=istation+1
  read(line,*)station_id,station_x
enddo

if(myrank==0)then
  ! plot VTK file: stations
  pfile=trim(out_path)//trim(file_head)//'_station'//'.vtk'
  open(100,file=trim(pfile),action='write',status='replace')

  write(100,'(a)')'# vtk DataFile Version 2.0'
  write(100,'(a)')'Unstructured Grid Example'
  write(100,'(a)')'ASCII'
  write(100,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(100,'(a,i5,a)')'POINTS',nstation,' float'
  do i_p=1,nstation
    write(100,'(3(f14.6,1x))')DIM_L*station_x(:,i_p)
  enddo
  write(100,*)
  write(100,'(a,i5,a,i5)')'CELLS',nstation,' ',2*nstation
  do i_p=1,nstation
    ! VTK format indexing starts from 0
    write(100,'(5(i5,1x))')1,i_p-1
  enddo
  write(100,*)
  write(100,'(a,i5)')'CELL_TYPES',nstation
  do i_p=1,nstation
    write(100,'(i1)')1
  enddo
  close(100)
  flush(100)
endif

! get derivatives of shape functions for 8-noded hex
allocate(dshape_hex8(ngnode,ndim))
allocate(lagrange_gll(ngll),dlagrange_gll(ndim,ngll))

imid=(ngll+1)/2

allocate(irec_located(nstation))
irec_located=0
! compute contribution of the earthquake station/s to the elemental loads
! number of stations
allocate(isnode(nnode),iselmt(nelmt))
maxerr=1e-30_kreal
minerr=INFTOL
station: do i_rec=1,nstation
  recx=station_x(:,i_rec)
  ! find the element which contains this earthquake station
  is_located=.false.
  
  ! First, find the element which is not too far and has the nearest GLL point
  ! from the target point.
  ielmt_min=1
  inode_min=g_num(1,ielmt_min)
  minsqdist=INFTOL
  do i_elmt=1,nelmt
    ! Do not locate stations in transition or infinite elements.
    mdomain=mat_domain(mat_id(i_elmt))
    if(mdomain==ELASTIC_TRINFDOMAIN .or.      &
       mdomain==ELASTIC_INFDOMAIN   .or.      &
       mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
       mdomain==VISCOELASTIC_INFDOMAIN)cycle
    ! Preliminary test. Discard far enough elements.
    num=g_num(:,i_elmt)
    xp=g_coord(:,num(imid))
    if(sqdistance(recx,xp,ndim).gt.sqmaxsize_elmt)cycle

    do i_gll=1,ngll
      xp=g_coord(:,num(i_gll))
      sqdist=sqdistance(recx,xp,ndim)
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
  nelmt_rectry=1
  do i_elmt=1,nelmt
    if(i_elmt.eq.ielmt_min)cycle
    if(any(isnode(g_num(hex8_gnode,i_elmt))))then
      nelmt_rectry=nelmt_rectry+1
      iselmt(i_elmt)=.true.
    endif
  enddo
  ! Find and assign the list of elements to try.
  allocate(ielmt_rectry(nelmt_rectry))
  ielmt_rectry=-9999
  inum=0
  do i_elmt=1,nelmt
    if(iselmt(i_elmt))then
      inum=inum+1
      ielmt_rectry(inum)=i_elmt
    endif
  enddo
  if(inum.ne.nelmt_rectry)then
    write(*,*)'ERROR: nelmt_rectry mismatch!'
    stop
  endif
  
  ! Loop throught the list of elements
  n_relmt=0
  minsqdist=INFTOL
  element_try: do i_elmt=1,nelmt_rectry
    ielmt=ielmt_rectry(i_elmt)

    ! preliminary test. discard far enough elements
    num=g_num(:,ielmt)
    coord=g_coord(:,num(hex8_gnode))
    ! check if the station is in this element
    call IsPointInHexahedron(coord,recx,isinside)
    if(isinside)then
      n_relmt=n_relmt+1
      ! map station location to natural coordinates
      call  map_point2naturalhex8(coord,recx,xip,located_x,niter,errx)
      errxd=errx*DIM_L
      if(n_relmt==1)then
        ! initialize errors
        maxerr=errxd
        minerr=errxd
      elseif(n_relmt.gt.1)then
        if(errxd.gt.maxerr)maxerr=errxd 
        if(errxd.lt.minerr)minerr=errxd 
      endif
      !write(*,'(a,i0,1x,g0.6)')'station location => niter and error (m): ',niter,errxd
      if(errxd.gt.maxsize_elmt)then
        write(logunit,'(a)')'WARNING: this station location may be inaccurate!'
        write(logunit,'(3(g0.6,1x))')xip
        flush(logunit)
      endif
      
      sqdist=sqdistance(recx,located_x,ndim)
      if(sqdist.lt.minsqdist)then
        r_elmt=ielmt
        recxi=xip
      endif
      is_located=.true.
      irec_located(i_rec)=1
      ! get deriative of shape function at station location
      call dshape_function_hex8p(ngnode,recxi(1),recxi(2), &
                                 recxi(3),dshape_hex8)
    endif
  enddo element_try ! i_elmt
  deallocate(ielmt_rectry)
  call sync_process
  n_relmt=sumscal(n_relmt)

  if(is_located)then
    if(n_relmt.gt.8)then
      write(logunit,'(a,i0,1x,i0)')'WARNING: number of elements shared by the &
      &station: ',i_rec,n_relmt
      flush(logunit)
    endif
  endif
  this_rec_located=sumscal(irec_located(i_rec))
  !if(this_rec_located.gt.0)then
  !  if(myrank==0)then
  !    write(*,'(a,i0,1x,i0)')'number of elements shared by the &
  !    &station: ',i_rec,n_relmt
  !  endif
  !endif
  if(this_rec_located.lt.1)then
    write(logunit,'(a,3(g0.6,1x))')'WARNING: station cannot be located: ',recx
    flush(logunit)
  endif

  ! find a unique element and a unique processor for this station
  ! and broadcast the results so that all processes are aware of this.

enddo station ! i_rec
deallocate(isnode,iselmt)
irec_located=maxvec(irec_located,nstation)
where(irec_located.gt.1)irec_located=1
total_rec_located=sum(irec_located)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0)')'total defined stations: ',nstation, &
                           'total located stations: ',total_rec_located
  flush(logunit)
endif
if(total_rec_located.gt.0)then
  maxerrg=maxscal(maxerr)
  minerrg=minscal(minerr)
  if(myrank==0)then
    write(logunit,'(a,g0.6,a,g0.6)')'station location errors (m) => min: ', &
    minerrg, ' max: ',maxerrg
    flush(logunit)
  endif
endif
deallocate(dshape_hex8)
deallocate(lagrange_gll,dlagrange_gll)

deallocate(station_x)

errcode=0

return

end subroutine locate_station
!===============================================================================

subroutine compute_station(lneq,errcode,errtag)
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
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt,i_line,ielmt,imid,inum,ios,istat
integer :: num(nenode),egdofu(nedofu)

real(kind=kreal) :: coord(ndim,8),jac(ndim,ndim),xp(ndim),xip(ndim)
real(kind=kreal) :: detjac

integer :: i_w,i_l,i_p,i_rec,i0,i1,i2,i3,i4,i0_nextrow
real(kind=kreal) :: pfac,rk,rl,rw
real(kind=kreal) :: located_x(ndim),recx(ndim),recxi(ndim)
real(kind=kreal) :: deriv(ndim,ngll)
logical :: is_fnpatch

integer :: mdomain

logical :: is_located
integer :: this_rec_located,total_rec_located
integer,allocatable:: irec_located(:)
integer :: niter
real(kind=kreal) :: errx,errxd,minerr,maxerr,minerrg,maxerrg

integer :: np,np_l,np_w

integer :: r_elmt,n_relmt

real(kind=kreal),dimension(:,:),allocatable :: dshape_hex8
real(kind=kreal),dimension(:),allocatable :: lagrange_gll
real(kind=kreal),dimension(:,:),allocatable :: dlagrange_gll

integer :: i_gll,ielmt_min,inode_min
integer :: nelmt_rectry
integer,allocatable :: ielmt_rectry(:)
logical,allocatable :: isnode(:),iselmt(:)
real(kind=kreal) :: minsqdist,sqdist

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

! get derivatives of shape functions for 8-noded hex
allocate(dshape_hex8(ngnode,ndim))
allocate(lagrange_gll(ngll),dlagrange_gll(ndim,ngll))

imid=(ngll+1)/2

allocate(irec_located(nstation))
irec_located=0
! compute contribution of the earthquake station/s to the elemental loads
! number of stations
allocate(isnode(nnode),iselmt(nelmt))
maxerr=1e-30_kreal
minerr=INFTOL
station: do i_rec=1,nstation
  recx=station_x(:,i_rec)
  ! find the element which contains this earthquake station
  is_located=.false.
  
  ! First, find the element which is not too far and has the nearest GLL point
  ! from the target point.
  ielmt_min=1
  inode_min=g_num(1,ielmt_min)
  minsqdist=INFTOL
  do i_elmt=1,nelmt
    ! Do not locate stations in transition or infinite elements.
    mdomain=mat_domain(mat_id(i_elmt))
    if(mdomain==ELASTIC_TRINFDOMAIN .or.      &
       mdomain==ELASTIC_INFDOMAIN   .or.      &
       mdomain==VISCOELASTIC_TRINFDOMAIN .or. &
       mdomain==VISCOELASTIC_INFDOMAIN)cycle
    ! Preliminary test. Discard far enough elements.
    num=g_num(:,i_elmt)
    xp=g_coord(:,num(imid))
    if(sqdistance(recx,xp,ndim).gt.sqmaxsize_elmt)cycle

    do i_gll=1,ngll
      xp=g_coord(:,num(i_gll))
      sqdist=sqdistance(recx,xp,ndim)
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
  nelmt_rectry=1
  do i_elmt=1,nelmt
    if(i_elmt.eq.ielmt_min)cycle
    if(any(isnode(g_num(hex8_gnode,i_elmt))))then
      nelmt_rectry=nelmt_rectry+1
      iselmt(i_elmt)=.true.
    endif
  enddo
  ! Find and assign the list of elements to try.
  allocate(ielmt_rectry(nelmt_rectry))
  ielmt_rectry=-9999
  inum=0
  do i_elmt=1,nelmt
    if(iselmt(i_elmt))then
      inum=inum+1
      ielmt_rectry(inum)=i_elmt
    endif
  enddo
  if(inum.ne.nelmt_rectry)then
    write(*,*)'ERROR: nelmt_rectry mismatch!'
    stop
  endif
  
  ! Loop throught the list of elements
  n_relmt=0
  element_try: do i_elmt=1,nelmt_rectry
    ielmt=ielmt_rectry(i_elmt)

    ! preliminary test. discard far enough elements
    num=g_num(:,ielmt)
    coord=g_coord(:,num(hex8_gnode))
    ! check if the station is in this element
    call IsPointInHexahedron(coord,recx,isinside)
    if(isinside)then
      n_relmt=n_relmt+1
      ! map station location to natural coordinates
      call  map_point2naturalhex8(coord,recx,xip,located_x,niter,errx)
      errxd=errx*DIM_L
      if(n_relmt==1)then
        ! initialize errors
        maxerr=errxd
        minerr=errxd
      elseif(n_relmt.gt.1)then
        if(errxd.gt.maxerr)maxerr=errxd 
        if(errxd.lt.minerr)minerr=errxd 
      endif
      !write(*,'(a,i0,1x,g0.6)')'station location => niter and error (m): ',niter,errxd
      if(errxd.gt.maxsize_elmt)then
        write(logunit,'(a)')'WARNING: this station location may be inaccurate!'
        write(logunit,'(3(g0.6,1x))')xip
        flush(logunit)
      endif

      r_elmt=ielmt
      recxi=xip
      is_located=.true.
      irec_located(i_rec)=1
      ! get deriative of shape function at station location
      call dshape_function_hex8p(ngnode,recxi(1),recxi(2), &
                                 recxi(3),dshape_hex8)

      ! compute jacobian
      !WAARNING: need to double check jacobian
      jac=matmul(transpose(dshape_hex8),transpose(coord))
      !jac=matmul(coord,dshape_hex8)
      detjac=determinant(jac)
      call invert(jac)

      ! compute GLL lagrange function and their derivative at station location
      call gll_lagrange3d_point(ndim,ngllx,nglly,ngllz,ngll, &
                                recxi(1),recxi(2),recxi(3),  &
                                lagrange_gll,dlagrange_gll)

      deriv=matmul(jac,dlagrange_gll) ! use der for gll

      egdofu=gdof_elmt(edofu,r_elmt)

      !exit element ! this will place the station in only one element
    endif
  enddo element_try ! i_elmt
  deallocate(ielmt_rectry)
  call sync_process
  n_relmt=sumscal(n_relmt)

  if(is_located)then
    if(n_relmt.gt.8)then
      write(logunit,'(a,i0,1x,i0)')'WARNING: number of elements shared by the &
      &station: ',i_rec,n_relmt
      flush(logunit)
    endif
  endif
  this_rec_located=sumscal(irec_located(i_rec))
  !if(this_rec_located.gt.0)then
  !  if(myrank==0)then
  !    write(*,'(a,i0,1x,i0)')'number of elements shared by the &
  !    &station: ',i_rec,n_relmt
  !  endif
  !endif
  if(this_rec_located.lt.1)then
    write(logunit,'(a,3(g0.6,1x))')'WARNING: station cannot be located: ',recx
    flush(logunit)
  endif
enddo station ! i_rec
deallocate(isnode,iselmt)
irec_located=maxvec(irec_located,nstation)
where(irec_located.gt.1)irec_located=1
total_rec_located=sum(irec_located)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0)')'total defined stations: ',nstation, &
                           'total located stations: ',total_rec_located
  flush(logunit)
endif
if(total_rec_located.gt.0)then
  maxerrg=maxscal(maxerr)
  minerrg=minscal(minerr)
  if(myrank==0)then
    write(logunit,'(a,g0.6,a,g0.6)')'station location errors (m) => min: ', &
    minerrg, ' max: ',maxerrg
    flush(logunit)
  endif
endif
deallocate(dshape_hex8)
deallocate(lagrange_gll,dlagrange_gll)

deallocate(station_x)

errcode=0

return

end subroutine compute_station
!===============================================================================

end module station
!===============================================================================
