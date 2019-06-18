! AUTHOR
!   Hom Nath Gharti
! REVISION:
!  HNG, Mar 11,2011; HNG, Apr 09,2010
! TODO:
!  define the array in main program such that entire array can be written
!  without transpose in this routine (done!)
module visual
use global,only:NDIM
implicit none
character(len=20), parameter :: wild_char='********************'
! Element types
character(len=20) :: ensight_hex8='hexa8'
character(len=20) :: ensight_quad4='quad4'

! private
character(len=250),private :: myfname=' => visual.f90'
character(len=500),private :: errsrc

contains
!-------------------------------------------------------------------------------

subroutine write_ensight_casefile(case_file,geo_file,add_tag,errcode,errtag)
use global,only:file_head,ptail,savedata,ismassdens,ismagnetization
implicit none
character(len=250),intent(in) :: case_file,geo_file
character(len=60),intent(in) :: add_tag 
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

character(len=250) :: file_tag
integer :: ios

errtag=""
errcode=0

open(unit=11,file=trim(case_file),status='replace',action='write',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(case_file)//'" cannot be opened!'
  errcode=-1
  return
endif

write(11,'(a)')'FORMAT'
write(11,'(a,/)')'type:  ensight gold'

write(11,'(a)')'GEOMETRY'
write(11,'(a,a,/)')'model:    ',trim(geo_file)

if(savedata%model_cell)then
  write(11,'(a)')'VARIABLE'

  ! file tag
  file_tag=trim(file_head)//trim(add_tag)

  if(ismassdens)then
    write(11,'(a/)')'scalar per element: mass_density '//trim(file_tag)//   &
    trim(ptail)//'.rho'
  endif

  if(ismagnetization)then
    write(11,'(a/)')'scalar per element: magnetization '//trim(file_tag)//   &
    trim(ptail)//'.mag'
  endif
endif
close(11)
end subroutine write_ensight_casefile
!===============================================================================

subroutine write_ensight_casefile_long(case_file,geo_file,add_tag,isgeo_change,&
ts,ns,fs,fi,twidth,errcode,errtag,freesurf,isplane)
use global,only:file_head,ptail,savedata,benchmark_okada,dtstep
implicit none
character(len=250),intent(in) :: case_file,geo_file
character(len=60),intent(in) :: add_tag 
integer,intent(in) :: ts,ns,fs,fi,twidth
logical,intent(in) :: isgeo_change
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag
logical,intent(in),optional :: freesurf,isplane
character(len=250) :: file_tag
integer :: i,ios
logical :: isfreesurf

errtag=""
errcode=0

! Determine if it is for the free surface
isfreesurf=.false.
if(present(freesurf))then
  if(freesurf)isfreesurf=.true.
endif

! file tag
file_tag=trim(file_head)//trim(add_tag)

open(unit=11,file=trim(case_file),status='replace',action='write',iostat = ios)
if( ios /= 0 ) then
  write(errtag,'(a)')'ERROR: file "'//trim(case_file)//'" cannot be opened!'
  errcode=-1
  return
endif

write(11,'(a)')'FORMAT'
write(11,'(a,/)')'type:  ensight gold'

write(11,'(a)')'GEOMETRY'

if(isgeo_change)then
  write(11,'(a,i0,a,a/)')'model:    ',ts,' ',trim(geo_file)
else
  write(11,'(a,a/)')'model:    ',trim(geo_file)
endif

write(11,'(a)')'VARIABLE'

if(present(isplane))then
  if(isplane)then
    write(11,'(a/)')'scalar per node: z_coordinate '//trim(file_tag)//   &
    trim(ptail)//'.z'
  endif
endif

! Do not plot model for the free surface
if(savedata%model .and. .not.isfreesurf)then
  write(11,'(a/)')'scalar per node: bulk_modulus '//trim(file_tag)//   &
  trim(ptail)//'.kappa'
  write(11,'(a/)')'scalar per node: shear_modulus '//trim(file_tag)//  &
  trim(ptail)//'.mu'
  write(11,'(a/)')'scalar per node: mass_density '//trim(file_tag)//   &
  trim(ptail)//'.rho'
endif

if(savedata%disp)then
  write(11,'(a,i0,a,a,a,a,/)')'vector per node: ',ts,' ','displacement',' ',  &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.dis'
  if(benchmark_okada)then
    write(11,'(a,i0,a,a,a,a,/)')'vector per node: ',ts,' ','displacement_okada',' ',  &
    trim(file_tag)//'_okada'//trim(ptail)//'.dis'
  endif
endif
if(savedata%stress)then
  write(11,'(a,i0,a,a,a,a,/)')'tensor symm per node: ',ts,' ','stress',' ',   &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.sig'
endif
if(savedata%strain)then
  write(11,'(a,i0,a,a,a,a,/)')'tensor symm per node: ',ts,' ','strain',' ',   &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.eps'
endif
if(savedata%psigma)then
  write(11,'(a,i0,a,a,a,a,/)')'vector per node: ',ts,' ','principal_stress',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.psig'
endif
if(savedata%porep)then
   write(11,'(a,a,a,a,a,/)')'scalar per node: ',' ','pore_pressure',' ', &
   trim(file_tag)//trim(ptail)//'.por'
endif
if(savedata%vmeps)then
  write(11,'(a,i0,a,a,a,a,/)')'scalar per node: ',ts,' ','plastic_strain',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.eps'
endif
if(savedata%scf)then
  write(11,'(a,i0,a,a,a,a,/)')'scalar per node: ',ts,' ','stress_concentration_factor',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.scf'
endif
if(savedata%maxtau)then
  write(11,'(a,i0,a,a,a,a,/)')'scalar per node: ',ts,' ','max_shear_stress',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.mtau'
endif
if(savedata%nsigma)then
  write(11,'(a,i0,a,a,a,a,/)')'scalar per node: ',ts,' ','normal_stress',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.nsig'
endif
if(savedata%gpot)then
  write(11,'(a,i0,a,a,a,a,/)')'scalar per node: ',ts,' ','gravity_potential',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.gpot'
endif
if(savedata%agrav)then
  write(11,'(a,i0,a,a,a,a,/)')'vector per node: ',ts,' ','gravity_acceleration',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.grav'
endif
if(savedata%mpot)then
  write(11,'(a,i0,a,a,a,a,/)')'scalar per node: ',ts,' ','magnetic_potential',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.mpot'
endif
if(savedata%magb)then
  write(11,'(a,i0,a,a,a,a,/)')'vector per node: ',ts,' ','magnetic_field',' ', &
  trim(file_tag)//'_step'//wild_char(1:twidth)//trim(ptail)//'.magb'
endif
write(11,'(a)')'TIME'
write(11,'(a,i0)')'time set: ',ts
write(11,'(a,i0)')'number of steps: ',ns
write(11,'(a,i0)')'filename start number: ',fs
write(11,'(a,i0)')'filename increment: ',fi
write(11,'(a)',advance='no')'time values: '

!if(nexcav==0)then
!  do i=1,ns
!    write(11,'(es12.5)',advance='yes')srf(i)
!  enddo
!else
  do i=0,ns-1
    write(11,'(es12.5)',advance='yes')real(i)*dtstep
  enddo
!endif
close(11)
end subroutine write_ensight_casefile_long
!===============================================================================

! This subroutine writes an ensight geofile only upto coordinates and returns
! the file unit to the calling program so that the calling program can writes
! the remaining part (connectivity) of the geo file and close it.
subroutine write_ensight_geocoord(out_fname,ipart,spart,nnode,coord,funit)
implicit none
character(len=250),intent(in) :: out_fname
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: nnode
real,dimension(:,:),intent(in) :: coord !3D
integer,intent(out) :: funit

character(len=80) :: buffer ! this must be 80 characters long
character(len=250) :: myname
integer :: ios

myname='=>write_ensight_geocoord'

if(ubound(coord,1).ne.3 .or. ubound(coord,2).ne.nnode)then
  write(*,*)'ERROR: wrong coord dimensions!'
  stop
endif

funit=11
open(unit=funit,file=trim(out_fname),access='stream',form='unformatted',       &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer='C Binary'
write(funit)buffer
buffer='Created by write_ensight Routine'
write(funit)buffer
buffer='semfem3d'
write(funit)buffer
buffer='node id off'
write(funit)buffer
buffer='element id off'
write(funit)buffer
buffer='part'
write(funit)buffer
write(funit)ipart(1)
buffer=spart(1)
write(funit)buffer
buffer='coordinates'
write(funit)buffer
write(funit)nnode
write(funit)transpose(coord)
! do not close funit here
return
end subroutine write_ensight_geocoord
!===============================================================================

! This subroutine writes an ensight geo file that consists of the mesh
! information
subroutine write_ensight_geo(out_fname,etype,ipart,spart,nelmt,nnode,coord,   &
connect)
implicit none
character(len=250),intent(in) :: out_fname
character(len=20),intent(in) :: etype
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: nelmt,nnode
real,dimension(:,:),intent(in) :: coord !3D
integer,dimension(:,:),intent(in) :: connect !hexahedral elements

character(len=80) :: buffer ! this must be 80 characters long
character(len=250) :: myname
integer :: ios

myname='=>write_ensight_geo'

if(ubound(coord,1).ne.3 .or. ubound(coord,2).ne.nnode)then
  write(*,*)'ERROR: wrong "coord" dimensions!'//trim(myname),ubound(coord)
  stop
endif
if(ubound(connect,1).ne.8 .or. ubound(connect,2).ne.nelmt)then
  write(*,*)'ERROR: wrong "connect" dimensions!'
  stop
endif

open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer='C Binary'
write(11)buffer
buffer='Created by write_ensight Routine'
write(11)buffer
buffer='semfem3d'
write(11)buffer
buffer='node id off'
write(11)buffer
buffer='element id off'
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(1)
buffer=spart(1)
write(11)buffer
buffer='coordinates'
write(11)buffer
write(11)nnode
write(11)transpose(coord)
! writes element information
buffer=etype
write(11)buffer
write(11)nelmt

! do not substract 1 for ensight file
write(11)connect
close(11)
return
end subroutine write_ensight_geo
!===============================================================================

! This subroutine writes an ensight geo file that consists of the mesh
! information in multi-blocks (parts). This subroutine write only a particular
! block/part defined by the parameters. 
subroutine write_ensight_geocoord_part1(out_fname,ipart,spart,pindex, &
nnode1,node1,nnode,coord,funit)
implicit none
character(len=250),intent(in) :: out_fname
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: pindex
integer,intent(in) :: nnode1
integer,intent(in) :: node1(:)
integer,intent(in) :: nnode
real,dimension(:,:),intent(in) :: coord !3D
integer,intent(out) :: funit

character(len=80) :: buffer ! this must be 80 characters long
character(len=250) :: myname
integer :: ios

myname='=>write_ensight_geo_part2'

if(ubound(coord,1).ne.3 .or. ubound(coord,2).ne.nnode)then
  write(*,*)'ERROR: wrong "coord" dimensions!'//trim(myname),ubound(coord)
  stop
endif

open(newunit=funit,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer='C Binary'
write(funit)buffer
buffer='Created by write_ensight Routine'
write(funit)buffer
buffer='semfem3d'
write(funit)buffer
buffer='node id off'
write(funit)buffer
buffer='element id off'
write(funit)buffer
buffer='part'
write(funit)buffer
write(funit)ipart(pindex)
buffer=spart(pindex)
write(funit)buffer
buffer='coordinates'
write(funit)buffer
write(funit)nnode1
if(nnode1.gt.0)then
  write(funit)transpose(coord(:,node1))
endif
! do not close funit here
return
end subroutine write_ensight_geocoord_part1
!===============================================================================

! This subroutine writes an ensight geo file that consists of the mesh
! information in multi-blocks (parts). This subroutine write only a particular
! block/part defined by the parameters. 
subroutine write_ensight_geocoord_plane_part1(out_fname,ipart,spart,pindex, &
iplane,nnode1,node1,nnode,coord,funit)
implicit none
character(len=250),intent(in) :: out_fname
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: pindex
integer,intent(in) :: iplane
integer,intent(in) :: nnode1
integer,intent(in) :: node1(:)
integer,intent(in) :: nnode
real,dimension(:,:),intent(in) :: coord !3D
integer,intent(out) :: funit

real,parameter :: RZERO=0.0
real,dimension(:,:),allocatable :: pcoord !3D
character(len=80) :: buffer ! this must be 80 characters long
character(len=250) :: myname
integer :: ios

myname='=>write_ensight_geo_part2'

if(ubound(coord,1).ne.3 .or. ubound(coord,2).ne.nnode)then
  write(*,*)'ERROR: wrong "coord" dimensions!'//trim(myname),ubound(coord)
  stop
endif

open(newunit=funit,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer='C Binary'
write(funit)buffer
buffer='Created by write_ensight Routine'
write(funit)buffer
buffer='semfem3d'
write(funit)buffer
buffer='node id off'
write(funit)buffer
buffer='element id off'
write(funit)buffer
buffer='part'
write(funit)buffer
write(funit)ipart(pindex)
buffer=spart(pindex)
write(funit)buffer
buffer='coordinates'
write(funit)buffer
write(funit)nnode1
if(nnode1.gt.0)then
  allocate(pcoord(3,nnode1))
  pcoord=coord(:,node1)
  ! Set the coordinate for the iplane to ZERO
  pcoord(iplane,:)=RZERO
  ! Write coordinates
  write(funit)transpose(pcoord)
  deallocate(pcoord)
endif
! do not close funit here
return
end subroutine write_ensight_geocoord_plane_part1
!===============================================================================

! This subroutine writes an ensight geo file that consists of the mesh
! information in multi-blocks (parts). This subroutine write only a particular
! block/part defined by the parameters. 
subroutine write_ensight_geo_part1(out_fname,etype,ipart,spart,pindex, &
nelmt1,nnode1,node1,nnode,coord,connect1)
implicit none
character(len=250),intent(in) :: out_fname
character(len=20),intent(in) :: etype
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: pindex
integer,intent(in) :: nelmt1,nnode1
integer,intent(in) :: node1(:)
integer,intent(in) :: nnode
real,dimension(:,:),intent(in) :: coord !3D
integer,dimension(:,:),intent(in) :: connect1 !hexahedral elements

character(len=80) :: buffer ! this must be 80 characters long
character(len=250) :: myname
integer :: ios

myname='=>write_ensight_geo_part1'

if(ubound(coord,1).ne.3 .or. ubound(coord,2).ne.nnode)then
  write(*,*)'ERROR: wrong "coord" dimensions!'//trim(myname),ubound(coord)
  stop
endif
if(ubound(connect1,1).ne.8 .or. ubound(connect1,2).ne.nelmt1)then
  write(*,*)'ERROR: wrong "connect1" dimensions!'
  stop
endif

open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer='C Binary'
write(11)buffer
buffer='Created by write_ensight Routine'
write(11)buffer
buffer='semfem3d'
write(11)buffer
buffer='node id off'
write(11)buffer
buffer='element id off'
write(11)buffer
!call write_string(extent_str//char(0),fd)
!do j=1,3
!  do i=1,2
!    call write_float(real(extent(i,j)),fd)
!  enddo
!enddo
buffer='part'
write(11)buffer
write(11)ipart(pindex)
buffer=spart(pindex)
write(11)buffer
buffer='coordinates'
write(11)buffer
write(11)nnode1
write(11)transpose(coord(:,node1))
! writes element information
buffer=etype
write(11)buffer
write(11)nelmt1
! do not substract 1 for ensight file
write(11)connect1
close(11)
return
end subroutine write_ensight_geo_part1
!===============================================================================

! This subroutine writes an ensight geo file that consists of the mesh
! information in two-blocks (parts)
subroutine write_ensight_geo_part2(out_fname,etype,npart,ipart,spart, &
nelmt1,nelmt2,nnode1,nnode2,node1,node2,nnode,coord,connect1,connect2)
implicit none
character(len=250),intent(in) :: out_fname
character(len=20),intent(in) :: etype
integer,intent(in) :: npart
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: ipart(:)
integer,intent(in) :: nelmt1,nelmt2,nnode1,nnode2
integer,intent(in) :: node1(:),node2(:)
integer,intent(in) :: nnode
real,dimension(:,:),intent(in) :: coord !3D
integer,dimension(:,:),intent(in) :: connect1,connect2 !hexahedral elements

character(len=80) :: buffer ! this must be 80 characters long
character(len=250) :: myname
integer :: i_part,ios

myname='=>write_ensight_geo_part2'

if(ubound(coord,1).ne.3 .or. ubound(coord,2).ne.nnode)then
  write(*,*)'ERROR: wrong "coord" dimensions!'//trim(myname),ubound(coord)
  stop
endif
if(ubound(connect1,1).ne.8 .or. ubound(connect1,2).ne.nelmt1)then
  write(*,*)'ERROR: wrong "connect1" dimensions!'
  stop
endif
if(ubound(connect2,1).ne.8 .or. ubound(connect2,2).ne.nelmt2)then
  write(*,*)'ERROR: wrong "connect2" dimensions!'
  stop
endif

open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer='C Binary'
write(11)buffer
buffer='Created by write_ensight Routine'
write(11)buffer
buffer='semfem3d'
write(11)buffer
buffer='node id off'
write(11)buffer
buffer='element id off'
write(11)buffer
!call write_string(extent_str//char(0),fd)
!do j=1,3
!  do i=1,2
!    call write_float(real(extent(i,j)),fd)
!  enddo
!enddo
do i_part=1,npart
  buffer='part'
  write(11)buffer
  write(11)ipart(i_part)
  buffer=spart(i_part) !'unstructured meshes'
  write(11)buffer
  buffer='coordinates'
  write(11)buffer
  if(i_part==1)then
    write(11)nnode1
    write(11)transpose(coord(:,node1))
  else
    write(11)nnode2
    write(11)transpose(coord(:,node2))
  endif
  ! writes element information
  buffer=etype
  write(11)buffer
  if(i_part==1)then
    write(11)nelmt1
    ! do not substract 1 for ensight file
    write(11)connect1
  else
    write(11)nelmt2
    ! do not substract 1 for ensight file
    write(11)connect2
  endif
enddo
close(11)
return
end subroutine write_ensight_geo_part2
!===============================================================================

! This subroutines writes ensight gold per-node variable
! ncomp:1 = scalar, 3 = vector and 6 = symmetric tensor
subroutine write_ensight_pernode(out_fname,ipart,spart,ncomp,n,var)
implicit none
character(len=250),intent(in) :: out_fname
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: ncomp,n
real,intent(in) :: var(ncomp,n)

character(len=80) :: buffer ! this must be 80 characters long
integer :: ios

if(ncomp/=1 .and. ncomp/=3 .and. ncomp/=6)then
  write(*,'(a)')'ERROR: invalid ncomp for ensight per_node variable!'
  print*,ncomp
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(1)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(1)
buffer='coordinates'
write(11)buffer
!do i=1,ncomp
!  write(11)var(i,:)
!enddo
write(11)transpose(var)

close(11)
return
end subroutine write_ensight_pernode
!===============================================================================

! This subroutines writes ensight gold per-node variable
! ncomp:1 = scalar, 3 = vector and 6 = symmetric tensor
subroutine write_ensight_pernodeVECAS(out_fname,ipart,spart,ncomp,n,var)
implicit none
character(len=250),intent(in) :: out_fname
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: ncomp,n
real,intent(in) :: var(:,:)

character(len=80) :: buffer ! this must be 80 characters long
integer :: ios

if(ncomp/=1 .and. ncomp/=3 .and. ncomp/=6)then
  write(*,'(a)')'ERROR: invalid ncomp for ensight per_node variable!'
  print*,ncomp
  stop
endif
if(ubound(var,1).ne.ncomp .or. ubound(var,2).ne.n)then
  write(*,'(a)')'ERROR: wrong bounds of "var"!'
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(1)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(1)
buffer='coordinates'
write(11)buffer
!do i=1,ncomp
!  write(11)var(i,:)
!enddo
if(n.gt.0)then
  write(11)transpose(var)
endif

close(11)
return
end subroutine write_ensight_pernodeVECAS
!===============================================================================

! This subroutines writes ensight gold per-node variable
! ncomp:1 = scalar, 3 = vector and 6 = symmetric tensor
subroutine write_ensight_pernodeVECAS_part1(out_fname,ipart,spart,pindex, &
nnode1,node1,ncomp,n,var)
implicit none
character(len=250),intent(in) :: out_fname
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
! plot index
integer,intent(in) :: pindex
integer,intent(in) :: nnode1
integer,intent(in) :: node1(:)
integer,intent(in) :: ncomp,n
real,intent(in) :: var(:,:)

character(len=80) :: buffer ! this must be 80 characters long
integer :: ios

if(ncomp/=1 .and. ncomp/=3 .and. ncomp/=6)then
  write(*,'(a)')'ERROR: invalid ncomp for ensight per_node variable!'
  print*,ncomp
  stop
endif
if(ubound(var,1).ne.ncomp .or. ubound(var,2).ne.n)then
  write(*,'(a)')'ERROR: wrong bounds of "var"!'
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(pindex)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(pindex)
buffer='coordinates'
write(11)buffer
!do i=1,ncomp
!  write(11)var(i,:)
!enddo
if(nnode1.gt.0)then
  write(11)transpose(var(:,node1))
endif

close(11)
return
end subroutine write_ensight_pernodeVECAS_part1
!===============================================================================

! This subroutines writes ensight gold per-node variable
! ncomp:1 = scalar, 3 = vector and 6 = symmetric tensor
subroutine write_ensight_pernodeSCALAS(out_fname,ipart,spart,n,var)
implicit none
character(len=250),intent(in) :: out_fname
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: n
real,intent(in) :: var(:)

character(len=80) :: buffer ! this must be 80 characters long
integer :: ios

if(size(var).ne.n)then
  write(*,'(a)')'ERROR: wrong bounds of "var"!'
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(1)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(1)
buffer='coordinates'
write(11)buffer
!do i=1,ncomp
!  write(11)var(i,:)
!enddo
if(n.gt.0)then
  write(11)var
endif

close(11)
return
end subroutine write_ensight_pernodeSCALAS
!===============================================================================

! This subroutines writes ensight gold per-node variable
! ncomp:1 = scalar, 3 = vector and 6 = symmetric tensor
subroutine write_ensight_pernodeSCALAS_part1(out_fname,ipart,spart,pindex, &
nnode1,node1,n,var)
implicit none
character(len=250),intent(in) :: out_fname
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: pindex
integer,intent(in) :: nnode1
integer,intent(in) :: node1(:)
integer,intent(in) :: n
real,intent(in) :: var(:)

character(len=80) :: buffer ! this must be 80 characters long
integer :: i_node,ios

if(size(var).ne.n)then
  write(*,'(a)')'ERROR: wrong bounds of "var"!'
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(pindex)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(pindex)
buffer='coordinates'
write(11)buffer
if(nnode1.gt.0)then
  do i_node=1,nnode1
    write(11)var(node1(i_node))
  enddo
endif
! This will create temporary array
!write(11)var(node1)
close(11)
return
end subroutine write_ensight_pernodeSCALAS_part1
!===============================================================================

! This subroutines writes ensight gold per-elmt variable (SCALAR)
subroutine write_ensight_perelementSCALAS(out_fname,etype,ipart,spart,n,var)
implicit none
character(len=250),intent(in) :: out_fname
character(len=20),intent(in) :: etype
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: n
real,intent(in) :: var(:)

character(len=80) :: buffer ! this must be 80 characters long
integer :: ios

if(size(var).ne.n)then
  write(*,'(a)')'ERROR: wrong bounds of "var"!'
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(1)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(1)
buffer=etype
write(11)buffer
!do i=1,ncomp
!  write(11)var(i,:)
!enddo
write(11)var

close(11)
return
end subroutine write_ensight_perelementSCALAS
!===============================================================================

! This subroutines writes ensight gold per-elmt variable
! ncomp:1 = scalar, 3 = vector and 6 = symmetric tensor
subroutine write_ensight_perelementSCALAS_part1(out_fname,etype,ipart,spart,pindex, &
nelmt1,elmt1,n,var)
implicit none
character(len=250),intent(in) :: out_fname
character(len=20),intent(in) :: etype
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: pindex
integer,intent(in) :: nelmt1
integer,intent(in) :: elmt1(:)
integer,intent(in) :: n
real,intent(in) :: var(:)

character(len=80) :: buffer ! this must be 80 characters long
integer :: i_elmt,ios

if(size(var).ne.n)then
  write(*,'(a)')'ERROR: wrong bounds of "var"!'
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(pindex)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(pindex)
buffer=etype
write(11)buffer
do i_elmt=1,nelmt1
  write(11)var(elmt1(i_elmt))
enddo
! This will create temporary array
!write(11)var(elmt1)
close(11)
return
end subroutine write_ensight_perelementSCALAS_part1
!===============================================================================

! This subroutines writes ensight gold per-elmt variable (SCALAR)
subroutine write_ensight_perelementVECAS(out_fname,etype,ipart,spart,ncomp,n,var)
implicit none
character(len=250),intent(in) :: out_fname
character(len=20),intent(in) :: etype
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: ncomp,n
real,intent(in) :: var(:,:)

character(len=80) :: buffer ! this must be 80 characters long
integer :: ios

if(ncomp/=1 .and. ncomp/=3 .and. ncomp/=6)then
  write(*,'(a)')'ERROR: invalid ncomp for ensight per_node variable!'
  print*,ncomp
  stop
endif
if(ubound(var,1).ne.ncomp .or. ubound(var,2).ne.n)then
  write(*,'(a)')'ERROR: wrong bounds of "var"!'
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(1)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(1)
buffer=etype
write(11)buffer
!do i=1,ncomp
!  write(11)var(i,:)
!enddo
write(11)transpose(var)

close(11)
return
end subroutine write_ensight_perelementVECAS
!===============================================================================

! This subroutines writes ensight gold per-elmt variable
! ncomp:1 = scalar, 3 = vector and 6 = symmetric tensor
subroutine write_ensight_perelementVECAS_part1(out_fname,etype,ipart,spart,pindex, &
nelmt1,elmt1,ncomp,n,var)
implicit none
character(len=250),intent(in) :: out_fname
character(len=20),intent(in) :: etype
integer,intent(in) :: ipart(:)
character(len=80),intent(in) :: spart(:)
integer,intent(in) :: pindex
integer,intent(in) :: nelmt1
integer,intent(in) :: elmt1(:)
integer,intent(in) :: ncomp,n
real,intent(in) :: var(:,:)

character(len=80) :: buffer ! this must be 80 characters long
integer :: i_dim,i_elmt,ios

if(ncomp/=1 .and. ncomp/=3 .and. ncomp/=6)then
  write(*,'(a)')'ERROR: invalid ncomp for ensight per_node variable!'
  print*,ncomp
  stop
endif
if(ubound(var,1).ne.ncomp .or. ubound(var,2).ne.n)then
  write(*,'(a)')'ERROR: wrong bounds of "var"!'
  stop
endif
! open Ensight Gold data file to store data
open(unit=11,file=trim(out_fname),access='stream',form='unformatted',          &
status='replace',action='write',iostat=ios)
if (ios /= 0)then
  write(*,'(a)')'ERROR: output file "'//out_fname//'" cannot be opened!'
  stop
endif

buffer=spart(pindex)
write(11)buffer
buffer='part'
write(11)buffer
write(11)ipart(pindex)
buffer=etype
write(11)buffer
do i_dim=1,NDIM
  do i_elmt=1,nelmt1
    write(11)var(i_dim,elmt1(i_elmt))
  enddo
enddo
! This will create temporary array
!write(11)transpose(var(:,elmt1))
close(11)
return
end subroutine write_ensight_perelementVECAS_part1
!===============================================================================

end module visual
!===============================================================================
