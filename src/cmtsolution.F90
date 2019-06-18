module cmtsolution
implicit none
! Number of lines per CMT source in a CMTSOLUTION file
integer,parameter :: nline_cmtsolution=13
! Number of CMT sources in a CMTSOLUTION file
integer :: ncmt_source

contains
!-------------------------------------------------------------------------------

! DESCRIPTION
!  This routine counts the number of CMT sources in the CMTSOLUTION file.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
!  Leah Langer, Princeton University
! REVISION
!  HNG, Oct 04, 2018
! TODO
subroutine count_cmtsolution(errcode,errtag)
use global,only:inp_path,cmtfile
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer,parameter :: nmax_line=100 ! maximum number of lines in the eqsource
!file
integer :: ios
integer :: nline

character(len=1) :: tchar
character(len=250) :: line
character(len=80) :: fname
character(len=80) :: data_path

errtag="ERROR: unknown!"
errcode=-1

! set data path
data_path=trim(inp_path)

! open CMTSOLUTION file
fname=trim(data_path)//trim(cmtfile)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif

nline=0
do
  read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
  if (ios/=0)exit
  nline=nline+1
enddo
close(11)
if(mod(nline,nline_cmtsolution).ne.0)then
  write(errtag,*)'ERROR: invalid number of lines in the CMTSOLUTION file "'//trim(cmtfile)//'!'
  return
endif
ncmt_source=nline/nline_cmtsolution

errcode=0

return

end subroutine count_cmtsolution
!===============================================================================

! DESCRIPTION
!  This routine counts the number of CMT sources in the CMTSOLUTION file.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
!  Leah Langer, Princeton University
! REVISION
!  HNG, Oct 04, 2018
! TODO
subroutine read_cmtsolution(source_coord,M_cmt,errcode,errtag)
use dimensionless
use global
use math_constants
use conversion_constants
use math_library,only:IsPointInPolygon
use shape_library,only:shape_function_quad4p
use map_location,only:map_point2naturalquad4
use utmgeo
use free_surface
use string_library
#if (USE_MPI)
use mpi_library
use math_library_mpi
#else
use serial_library
use math_library_serial
#endif
implicit none
real(kind=kreal),intent(out) :: source_coord(:,:)
real(kind=kreal),intent(out) :: M_cmt(:,:)
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer,parameter :: nmax_line=100 ! maximum number of lines in the eqsource
!file
integer :: i_elmt,i_face,i_line,i_src,iface,iface_all,ind,ios,istat
integer :: niter,nline
integer :: nsrc

real(kind=kreal) :: lat,long,depth,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
real(kind=kreal),dimension(3,3) :: Mcmt

real(kind=kreal) :: mzz,myz,mxz,mxy,myy,mxx
real(kind=kreal) :: utmx(2)
real(kind=kreal) :: coord(2,4),vx(4),vy(4),vz(4)
real(kind=kreal) :: located_x(2),xip(2),errx
real(kind=kreal) :: shape_quad4(4)
real(kind=kreal) :: elevation

integer :: this_src_located,total_src_located
integer,allocatable :: isrc_located(:)
logical :: isinside
character(len=1) :: tchar
character(len=250) :: line,tag
character(len=80) :: token
character(len=80) :: fname
character(len=80) :: data_path
character(len=250) :: pfile

logical :: islat,islong,isdepth,ismrr,ismtt,ismpp,ismrt,ismrp,ismtp
errtag="ERROR: unknown!"
errcode=-1

! Set data path
data_path=trim(inp_path)

! open CMTSOLUTION file
fname=trim(data_path)//trim(cmtfile)
open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
if( ios /= 0 ) then
  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  return
endif

!! TESTING
!if(myrank==0)then
!  lat=34.213_kreal
!  long=-118.537_kreal
!  call geodetic2utm(long,lat,utmx(1),utmx(2))
!  write(*,*)'Northridge epicenter:',utmx
!endif
allocate(isrc_located(ncmt_source))
isrc_located=0
nsrc=0
src:do i_src=1,ncmt_source
  islat=.false.
  islong=.false.
  isdepth=.false.
  ismrr=.false.
  ismtt=.false.
  ismpp=.false.
  ismrt=.false.
  ismrp=.false.
  ismtp=.false.
  do i_line=1,nline_cmtsolution
    read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
    if (ios/=0)exit src
    !first token
    tag=trim(line)
    call first_token(tag,token)

    ! latitude
    if (trim(token)=='latitude:')then
      if(islat)then
        write(errtag,*)'WARNING: copy of latitude/latorUTM found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "latitude:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)lat
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read latitude!',ind,trim(line)
        stop
      endif
      islat=.true.
    endif
    if (trim(token)=='latorUTM:')then
      if(islat)then
        write(errtag,*)'WARNING: copy of latitude/latorUTM found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "latorUTM:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)lat
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read latorUTM!',ind,trim(line)
        stop
      endif
      islat=.true.
    endif
    ! longitude
    if (trim(token)=='longitude:')then
      if(islong)then
        write(errtag,*)'WARNING: copy of longitude/longorUTM found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "longitude:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)long
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read longitude!'
        stop
      endif
      islong=.true.
    endif
    if (trim(token)=='longorUTM:')then
      if(islong)then
        write(errtag,*)'WARNING: copy of longitude/longorUTM found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "langorUTM:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)long
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read longitude!'
        stop
      endif
      islong=.true.
    endif
    ! depth
    if (trim(token)=='depth:')then
      if(isdepth)then
        write(errtag,*)'WARNING: copy of depth found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "depth:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)depth
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read depth!'
        stop
      endif
      isdepth=.true.
    endif
    ! Mrr
    if (trim(token)=='Mrr:')then
      if(ismrr)then
        write(errtag,*)'WARNING: copy of Mrr found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "Mrr:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)Mrr
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read Mrr!'
        stop
      endif
      ismrr=.true.
    endif
    ! Mtt
    if (trim(token)=='Mtt:')then
      if(ismtt)then
        write(errtag,*)'WARNING: copy of Mtt found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "Mtt:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)Mtt
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read Mtt!'
        stop
      endif
      ismtt=.true.
    endif
    ! Mpp
    if (trim(token)=='Mpp:')then
      if(ismpp)then
        write(errtag,*)'WARNING: copy of Mpp found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "Mpp:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)Mpp
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read Mpp!'
        stop
      endif
      ismpp=.true.
    endif
    ! Mrt
    if (trim(token)=='Mrt:')then
      if(ismrt)then
        write(errtag,*)'WARNING: copy of Mrt found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "Mrt:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)Mrt
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read Mrt!'
        stop
      endif
      ismrt=.true.
    endif
    ! Mrp
    if (trim(token)=='Mrp:')then
      if(ismrp)then
        write(errtag,*)'WARNING: copy of Mrp found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "Mrp:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)Mrp
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read Mrp!'
        stop
      endif
      ismrp=.true.
    endif
    ! Mtp
    if (trim(token)=='Mtp:')then
      if(ismtp)then
        write(errtag,*)'WARNING: copy of Mtp found! Copies are discarded!'
        cycle
      endif
      ind=index(line,':')
      if(ind.le.0)then
        write(errtag,*)'ERROR: invalid "Mtp:" token!',trim(line)
        return
      endif
      read(line(ind+1:len_trim(line)),*,iostat=istat)Mtp
      if(istat.ne.0)then
        write(*,*)'ERROR: cannot read Mtp!'
        stop
      endif
      ismtp=.true.
    endif
  enddo
  ! Check status
  if(.not.islat)then
    write(errtag,*)'ERROR: cannot read latitude!'
    return
  endif
  if(.not.islong)then
    write(errtag,*)'ERROR: cannot read longitude!'
    return
  endif
  if(.not.isdepth)then
    write(errtag,*)'ERROR: cannot read depth!'
    return
  endif
  if(.not.ismrr)then
    write(errtag,*)'ERROR: cannot read Mrr!'
    return
  endif
  if(.not.ismtt)then
    write(errtag,*)'ERROR: cannot read Mtt!'
    return
  endif
  if(.not.ismpp)then
    write(errtag,*)'ERROR: cannot read Mpp!'
    return
  endif
  if(.not.ismrt)then
    write(errtag,*)'ERROR: cannot read Mrt!'
    return
  endif
  if(.not.ismrp)then
    write(errtag,*)'ERROR: cannot read Mrp!'
    return
  endif
  if(.not.ismtp)then
    write(errtag,*)'ERROR: cannot read Mtp!'
    return
  endif
  nsrc=nsrc+1

  ! Nondimentionalize M
  ! NOTE: unit of M is in CGS units. We should convert them to SI unit.
  Mpp=NONDIM_MTENS*CGS2SI_MOMENT*Mpp
  Mtt=NONDIM_MTENS*CGS2SI_MOMENT*Mtt
  Mrr=NONDIM_MTENS*CGS2SI_MOMENT*Mrr
  Mtp=NONDIM_MTENS*CGS2SI_MOMENT*Mtp
  Mrt=NONDIM_MTENS*CGS2SI_MOMENT*Mrt
  Mrp=NONDIM_MTENS*CGS2SI_MOMENT*Mrp
  
  !  Mrr =  Mzz
  !  Mtt =  Myy
  !  Mpp =  Mxx
  !  Mrt = -Myz
  !  Mrp =  Mxz
  !  Mtp = -Mxy

  ! Store 6 unique moment-tensor components in the order
  ! 1: Mxx
  ! 2: Myy
  ! 3: Mzz
  ! 4: Mxy
  ! 5: Myz
  ! 6: Mzx
  M_cmt(1,i_src)= Mpp !Mxx
  M_cmt(2,i_src)= Mtt !Myy
  M_cmt(3,i_src)= Mrr !Mzz
  M_cmt(4,i_src)=-Mtp !Mxy
  M_cmt(5,i_src)=-Mrt !Myz
  M_cmt(6,i_src)= Mrp !Mzx
 
  ! Compute UTM coordinates
  call geodetic2utm(long,lat,utmx(1),utmx(2))
 
  ! Nondimentionalize coordinates
  utmx=NONDIM_L*utmx
  depth=NONDIM_L*KM2M*depth

  source_coord(1,i_src)=utmx(1)
  source_coord(2,i_src)=utmx(2)
  !source_coord(3,i_src)=-depth

  ! Find Z coordinate
  ! step 1: find elevation
  ! NOTE: All the processors will try to find elevation, but only a single
  ! processor or processors shared by the point will find the correct elevation.
  ! Otherwise the elevation is set to ZERO.
  call free_surface_elevation(utmx,elevation,isrc_located(i_src))

  source_coord(3,i_src)=elevation-depth
  
  !iface=0
  !do i_face=1,nelmt_fs
  !  vx=g_coord(1,gnum4_fs(:,i_face))
  !  vy=g_coord(2,gnum4_fs(:,i_face))
  !  coord(1,:)=vx
  !  coord(2,:)=vy
  !  vz=g_coord(3,gnum4_fs(:,i_face))
  !  isinside=IsPointInPolygon(vx,vy,utmx(1),utmx(2))
  !  if(isinside)then
  !    iface=iface+1
  !    call  map_point2naturalquad4(coord,utmx,xip,located_x,niter,errx)
  !    call shape_function_quad4p(4,xip(1),xip(2),shape_quad4)
  !    elevation=sum(vz*shape_quad4)
  !    
  !    source_coord(3,i_src)=elevation-depth
  !    isrc_located(i_src)=1
  !  endif
  !enddo
  call sync_process
  this_src_located=sumscal(isrc_located(i_src))
  iface_all=sumscal(iface) 
  if(this_src_located.lt.1)then
    write(logunit,'(a,2(g0.6,1x),i0,1x,i0,1x,i0,1x,i0))')'WARNING: cmt source cannot be projected &
                                   &on the free surface:',utmx,i_src,myrank,iface_all,this_src_located
    flush(logunit)
  endif

enddo src
close(11)
call sync_process

isrc_located=maxvec(isrc_located,ncmt_source)
where(isrc_located.gt.1)isrc_located=1
total_src_located=sum(isrc_located)
if(myrank==0)then
  write(logunit,'(a,i0,1x,a,i0)')'Total defined CMT sources: ',ncmt_source, &
                           'Total projected CMT sources: ',total_src_located
  flush(logunit)
endif
if(nsrc.ne.ncmt_source)then
  write(errtag,*)'ERROR: number of CMT sources mismatch!'
  return
endif

deallocate(isrc_located)

errcode=0

return

end subroutine read_cmtsolution
!===============================================================================

end module cmtsolution
!===============================================================================
