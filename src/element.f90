! AUTHOR
!   Hom Nath Gharti
module element
use set_precision
implicit none
character(len=250),private :: myfname=" => element.f90"
character(len=500),private :: errsrc

! Local node numbering in each face CUBIT/EXODUS convention
integer,parameter,dimension(4,6) :: node_face=reshape( (/ 1,2,6,5, &
                                                          2,3,7,6, &
                                                          4,3,7,8, &
                                                          1,4,8,5, &
                                                          1,2,3,4, &
                                                          5,6,7,8/),(/4,6/) )

! Local node numbering in each edge CUBIT/EXODUS convention
integer,parameter,dimension(2,12) :: node_edge=reshape( (/ 1,2, &
                                                           2,3, &
                                                           3,4, &
                                                           4,1, &
                                                           1,5, &
                                                           2,6, &
                                                           3,7, &
                                                           4,8, &
                                                           5,6, &
                                                           6,7, &
                                                           7,8, &
                                                           8,5 /),(/2,12/) )

! Map sequential node numbering to exodus/cubit order for 8-noded hexahedra.
integer,parameter :: map2exodus_hex8(8)=(/ 1,2,4,3,5,6,8,7 /)
! Map sequential node numbering to exodus/cubit order for 4-noded quadrilateral.
integer,parameter :: map2exodus_quad4(4)=(/ 1,2,4,3 /)
integer :: hex8_gnode(8)
real(kind=kreal) :: hexface_sign(6)
! face sign or normal orientation (outward +, inward -)

type hex_face
  integer,allocatable :: node(:)
  integer :: gnode(4) ! geometric (corner) nodes only
  integer,allocatable :: edof(:)
end type hex_face
type (hex_face) :: hexface(6)

type hex_face_edge
  ! node index in ngllx * nglly * ngllz
  integer,allocatable :: node(:)
  ! node index in ngllx * nglly or nglly * ngllz, etc.
  integer,allocatable :: fnode(:)
end type hex_face_edge
type (hex_face_edge) :: hexface_edge(6,4) ! each of 6 HEX faces has 4 edges 
contains

!-------------------------------------------------------------------------------
subroutine prepare_hex(errcode,errtag)
use global,only:ngllx,ngllz,ngll,ngllxy
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

errtag="ERROR: unknown!"
errcode=-1

! geometrical nodes (corner nodes) in EXODUS/CUBIT order
! bottom nodes
hex8_gnode(1)=1;
hex8_gnode(2)=ngllx
hex8_gnode(3)=ngllxy;
hex8_gnode(4)=hex8_gnode(3)-ngllx+1
! top nodes
hex8_gnode(5)=(ngllz-1)*ngllxy+1;
hex8_gnode(6)=hex8_gnode(5)+ngllx-1
hex8_gnode(7)=ngll;
hex8_gnode(8)=hex8_gnode(7)-ngllx+1

errcode=0
return

end subroutine prepare_hex
!===============================================================================

subroutine prepare_hexface(errcode,errtag)
use global,only:ngllx,nglly,ngllz,ngllxy,ngllyz,ngllzx,nndof
use math_constants,only:ONE
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i,i1,i2,i3,i4,i5,i6,j,k
integer :: i_gll,idof1,idof2
integer :: i_face,inode
integer :: iedge,iface
integer :: jm1,nx,ny
integer :: nfgll
integer,allocatable :: indx(:),indy(:),indz(:)

errtag="ERROR: unknown!"
errcode=-1

allocate(hexface(1)%node(ngllzx),hexface(3)%node(ngllzx))
allocate(hexface(2)%node(ngllyz),hexface(4)%node(ngllyz))
allocate(hexface(5)%node(ngllxy),hexface(6)%node(ngllxy))

! local node numbers for the faces (faces are numbered in exodus/CUBIT
! convention)
inode=0
i1=0; i2=0; i3=0; i4=0; i5=0; i6=0
do k=1,ngllz
  do j=1,nglly
    do i=1,ngllx
      inode=inode+1
      if (i==1)then
        ! face 4
        i4=i4+1
        hexface(4)%node(i4)=inode
      endif
      if (i==ngllx)then
        ! face 2
        i2=i2+1
        hexface(2)%node(i2)=inode
      endif
      if (j==1)then
        ! face 1
        i1=i1+1
        hexface(1)%node(i1)=inode
      endif
      if (j==nglly)then
        ! face 3
        i3=i3+1
        hexface(3)%node(i3)=inode
      endif
      if (k==1)then
        ! face 5
        i5=i5+1
        hexface(5)%node(i5)=inode
      endif
      if (k==ngllz)then
        ! face 6
        i6=i6+1
        hexface(6)%node(i6)=inode
      endif
    enddo
  enddo
enddo

! find geometric corners nodes
do i_face=1,6 ! there are 6 faces in a hexahedron
  if(i_face==1 .or. i_face==3)then ! ZX plane
    hexface(i_face)%gnode(1)=hexface(i_face)%node(1)
    hexface(i_face)%gnode(2)=hexface(i_face)%node(ngllx)
    hexface(i_face)%gnode(3)=hexface(i_face)%node(ngllzx)
    hexface(i_face)%gnode(4)=hexface(i_face)%node(ngllzx-ngllx+1)
  elseif(i_face==2 .or. i_face==4)then ! YZ plane
    hexface(i_face)%gnode(1)=hexface(i_face)%node(1)
    hexface(i_face)%gnode(2)=hexface(i_face)%node(nglly)
    hexface(i_face)%gnode(3)=hexface(i_face)%node(ngllyz)
    hexface(i_face)%gnode(4)=hexface(i_face)%node(ngllyz-nglly+1)
  elseif(i_face==5 .or. i_face==6)then ! XY plane
    hexface(i_face)%gnode(1)=hexface(i_face)%node(1)
    hexface(i_face)%gnode(2)=hexface(i_face)%node(ngllx)
    hexface(i_face)%gnode(3)=hexface(i_face)%node(ngllxy)
    hexface(i_face)%gnode(4)=hexface(i_face)%node(ngllxy-ngllx+1)
  else
    write(errtag,'(a)')'ERROR: wrong face ID for traction!'
    return
  endif
enddo

! orientation of the normals
hexface_sign(1)=one
hexface_sign(2)=one
hexface_sign(6)=one
hexface_sign(3)=-one
hexface_sign(4)=-one
hexface_sign(5)=-one

! local degrees of freedoms on face
allocate(hexface(1)%edof(nndof*ngllzx),hexface(3)%edof(nndof*ngllzx))
allocate(hexface(2)%edof(nndof*ngllyz),hexface(4)%edof(nndof*ngllyz))
allocate(hexface(5)%edof(nndof*ngllxy),hexface(6)%edof(nndof*ngllxy))

do i_face=1,6
  if(i_face==1.or.i_face==3)then
    nfgll=ngllzx
  elseif(i_face==2.or.i_face==4)then
    nfgll=ngllyz
  elseif(i_face==5.or.i_face==6)then
    nfgll=ngllxy
  else
    write(errtag,'(a)')'ERROR: wrong hexface id!'
    return
  endif
  idof1=1; idof2=nndof
  do i_gll=1,nfgll
    !print*,idof1,idof2     
    hexface(i_face)%edof(idof1:idof2)=(/(i,i=1,nndof)/)+(hexface(i_face)%node(i_gll)-1)*nndof
    idof1=idof1+nndof; idof2=idof2+nndof
  enddo
enddo

! face_edge_nodes
allocate(indx(ngllx),indy(nglly),indz(ngllz))

!-------------------------------------------------------------------------------
! face 1 => [1,2,6,5]
iface=1
allocate(hexface_edge(1,1)%node(ngllx),hexface_edge(1,2)%node(ngllz), &
         hexface_edge(1,3)%node(ngllx),hexface_edge(1,4)%node(ngllz))

nx=ngllx
ny=ngllz

! edge 1 => [1,2]
iedge=1
do j=1,1
  jm1=j-1
  do i=1,nx
    indx(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indx
hexface_edge(iface,iedge)%node=hexface(iface)%node(indx)

! edge 2 => [2,6]
iedge=2
do j=1,ny
  jm1=j-1
  do i=nx,nx
    indy(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)

! edge 3 => [6,5]
iedge=3
do j=ny,ny
  jm1=j-1
  do i=1,nx
    indx(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indx
hexface_edge(iface,iedge)%node=hexface(iface)%node(indx)

! edge 4 => [5,1]
iedge=4
do j=1,ny
  jm1=j-1
  do i=1,1
    indy(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)
!-------------------------------------------------------------------------------

! face:2 => [2,3,7,6]
iface=2
allocate(hexface_edge(2,1)%node(nglly),hexface_edge(2,2)%node(ngllz), &
         hexface_edge(2,3)%node(nglly),hexface_edge(2,4)%node(ngllz))

nx=nglly
ny=ngllz

! edge 1 => [2,3]
iedge=1
do j=1,1
  jm1=j-1
  do i=1,nx
    indy(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)

! edge 2 => [3,7]
iedge=2
do j=1,ny
  jm1=j-1
  do i=nx,nx
    indz(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indz
hexface_edge(iface,iedge)%node=hexface(iface)%node(indz)

! edge 3 => [7,6]
iedge=3
do j=ny,ny
  jm1=j-1
  do i=1,nx
    indy(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)

! edge 4 => [6,2]
iedge=4
do j=1,ny
  jm1=j-1
  do i=1,1
    indz(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indz
hexface_edge(iface,iedge)%node=hexface(iface)%node(indz)
!-------------------------------------------------------------------------------

! face:3 => [4,3,7,8]
iface=3
allocate(hexface_edge(3,1)%node(ngllx),hexface_edge(3,2)%node(ngllz), &
         hexface_edge(3,3)%node(ngllx),hexface_edge(3,4)%node(ngllz))

nx=ngllx
ny=ngllz

! edge 1 => [4,3]
iedge=1
do j=1,1
  jm1=j-1
  do i=1,nx
    indx(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indx
hexface_edge(iface,iedge)%node=hexface(iface)%node(indx)

! edge 2 => [3,7]
iedge=2
do j=1,ny
  jm1=j-1
  do i=nx,nx
    indz(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indz
hexface_edge(iface,iedge)%node=hexface(iface)%node(indz)

! edge 3 => [7,8]
iedge=3
do j=ny,ny
  jm1=j-1
  do i=1,nx
    indx(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indx
hexface_edge(iface,iedge)%node=hexface(iface)%node(indx)

! edge 4 => [8,4]
iedge=4
do j=1,ny
  jm1=j-1
  do i=1,1
    indz(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indz
hexface_edge(iface,iedge)%node=hexface(iface)%node(indz)
!-------------------------------------------------------------------------------

! face:4 => [1,4,8,5]
iface=4
allocate(hexface_edge(4,1)%node(nglly),hexface_edge(4,2)%node(ngllz), &
         hexface_edge(4,3)%node(nglly),hexface_edge(4,4)%node(ngllz))

nx=nglly
ny=ngllz

! edge 1 => [1,4]
iedge=1
do j=1,1
  jm1=j-1
  do i=1,nx
    indy(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)

! edge 2 => [4,8]
iedge=2
do j=1,ny
  jm1=j-1
  do i=nx,nx
    indz(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indz
hexface_edge(iface,iedge)%node=hexface(iface)%node(indz)

! edge 3 => [8,5]
iedge=3
do j=ny,ny
  jm1=j-1
  do i=1,nx
    indy(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)

! edge 4 => [5,1]
iedge=4
do j=1,ny
  jm1=j-1
  do i=1,1
    indz(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indz
hexface_edge(iface,iedge)%node=hexface(iface)%node(indz)
!-------------------------------------------------------------------------------

! face:5 => [1,2,3,4]
iface=5
allocate(hexface_edge(5,1)%node(ngllx),hexface_edge(5,2)%node(nglly), &
         hexface_edge(5,3)%node(ngllx),hexface_edge(5,4)%node(nglly))

nx=ngllx
ny=nglly

! edge 1 => [1,4]
iedge=1
do j=1,1
  jm1=j-1
  do i=1,nx
    indx(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indx
hexface_edge(iface,iedge)%node=hexface(iface)%node(indx)

! edge 2 => [4,8]
iedge=2
do j=1,ny
  jm1=j-1
  do i=nx,nx
    indy(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)

! edge 3 => [8,5]
iedge=3
do j=ny,ny
  jm1=j-1
  do i=1,nx
    indx(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indx
hexface_edge(iface,iedge)%node=hexface(iface)%node(indx)

! edge 4 => [5,1]
iedge=4
do j=1,ny
  jm1=j-1
  do i=1,1
    indy(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)
!-------------------------------------------------------------------------------

! face:6 => [5,6,7,8]
iface=6
allocate(hexface_edge(6,1)%node(ngllx),hexface_edge(6,2)%node(nglly), &
         hexface_edge(6,3)%node(ngllx),hexface_edge(6,4)%node(nglly))

nx=ngllx
ny=nglly

! edge 1 => [5,6]
iedge=1
do j=1,1
  jm1=j-1
  do i=1,nx
    indx(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indx
hexface_edge(iface,iedge)%node=hexface(iface)%node(indx)

! edge 2 => [6,7]
iedge=2
do j=1,ny
  jm1=j-1
  do i=nx,nx
    indy(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)

! edge 3 => [7,8]
iedge=3
do j=ny,ny
  jm1=j-1
  do i=1,nx
    indx(i)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indx
hexface_edge(iface,iedge)%node=hexface(iface)%node(indx)

! edge 4 => [8,5]
iedge=4
do j=1,ny
  jm1=j-1
  do i=1,1
    indy(j)=jm1*nx+i
  enddo
enddo
hexface_edge(iface,iedge)%fnode=indy
hexface_edge(iface,iedge)%node=hexface(iface)%node(indy)
!-------------------------------------------------------------------------------

deallocate(indx,indy,indz)
errcode=0
return

end subroutine prepare_hexface
!===============================================================================

subroutine cleanup_hexface(errcode,errtag)
implicit none
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_edge,i_face

errtag="ERROR: unknown!"
errcode=-1

do i_face=1,6
  deallocate(hexface(i_face)%node)
  deallocate(hexface(i_face)%edof)
enddo
do i_face=1,6
  do i_edge=1,4
    deallocate(hexface_edge(i_face,i_edge)%node)
  enddo
enddo

errcode=0

return

end subroutine cleanup_hexface
!===============================================================================

! This function finds the node ID of a true node in isnode
function get_nodeidHEX8(inode1,gnodes) result(id)
implicit none
integer,intent(in) :: inode1
integer,dimension(8),intent(in) :: gnodes
integer :: i,id
do i=1,8
  if(inode1==gnodes(i))then
    id=i
    return
  endif
enddo
write(*,*)'ERROR: no common node for node ID!'
stop
end function get_nodeidHEX8
!===============================================================================

! This function finds the edge ID according to node_edge of 
! a set of true nodes in isnode
function get_edgeidHEX8(inode2,gnodes) result(id)
use global,only:myrank
implicit none
integer,dimension(2),intent(in) :: inode2
integer,dimension(8),intent(in) :: gnodes
integer :: i,id
integer,dimension(2) :: gnode2

do i=1,12
  gnode2=gnodes(node_edge(:,i))
  if(all(isort2(inode2)==isort2(gnode2)))then
    id=i
    return
  endif
enddo
write(*,*)myrank,'ERROR: no common node for edge ID!',inode2,gnodes
stop
end function get_edgeidHEX8
!===============================================================================

! This function finds the face ID according to node_face of
! a set of true nodes in isnode
function get_faceidHEX8(inode4,gnodes) result(id)
implicit none
integer,dimension(4),intent(in) :: inode4
integer,dimension(8),intent(in) :: gnodes
integer :: i,id
integer,dimension(4) :: gnode4

do i=1,6
  gnode4=gnodes(node_face(:,i))
  if(all(isort4(inode4)==isort4(gnode4)))then
    id=i
    return
  endif
enddo
write(*,*)'ERROR: no common node for face ID!'
stop
end function get_faceidHEX8
!===============================================================================

! This function sorts two values in ascending order
function isort2(v2) result(s2)
implicit none
integer,dimension(2),intent(in) :: v2
integer,dimension(2) :: s2
integer :: i,id

if(v2(1).gt.v2(2))then
  s2(1)=v2(2)
  s2(2)=v2(1)
else
  s2(1)=v2(1)
  s2(2)=v2(2)
endif
end function isort2
!===============================================================================

! This function sorts four values in ascending order
function isort4(v4) result(s4)
implicit none
integer,dimension(4),intent(in) :: v4
integer,dimension(4) :: s4
integer :: itmp,i1,i2

! sort first pair
if(v4(1).gt.v4(2))then
  itmp=v4(1)
  s4(1)=v4(2)
  s4(2)=itmp
else
  ! already sorted
  s4(1)=v4(1)
  s4(2)=v4(2)
endif
! sort second pair
if(v4(3).gt.v4(4))then
  itmp=v4(3)
  s4(3)=v4(4)
  s4(4)=itmp
else
  ! already sorted
  s4(3)=v4(3)
  s4(4)=v4(4)
endif

! assign first by comparing two minimum values
if(s4(1).lt.s4(3))then
  i1=s4(3)
  !s4(1)=s4(1)
else
  i1=s4(1)
  s4(1)=s4(3)
endif

! assign last by comparing two maximum values
if(s4(2).gt.s4(4))then
  i2=s4(4)
  s4(4)=s4(2)
else
  i2=s4(2)
  s4(4)=s4(4)
endif

! compare two remaining values
if(i1.lt.i2)then
  s4(2)=i1
  s4(3)=i2
else
  s4(2)=i2
  s4(3)=i1
endif

end function isort4
!===============================================================================

end module element
!===============================================================================
