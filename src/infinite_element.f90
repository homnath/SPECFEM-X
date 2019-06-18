! DESCRIPTION
!  This module conains the infinite-element preprocessing routines.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Apr 11,2012; HNG, Jul 12,2011; HNG, Apr 09,2010
! TODO
!  - 
module infinite_element
use set_precision
use global,only : ndim,INFINITE_RADAU,INFINITE_GAUSS,logunit
contains
!-------------------------------------------------------------------------------

! This subroutine classfies the elements into a finite, a transition infinite, &
! and an infinite regions.
! This classification is mainly used for the visulization purpose.
subroutine classify_finite_infinite_elements(ipass)
use global,only:nenode,ngnode,nnode,nelmt,mat_domain,mat_id,g_num, &
ELASTIC_TRINFDOMAIN,ELASTIC_INFDOMAIN, &
VISCOELASTIC_TRINFDOMAIN,VISCOELASTIC_INFDOMAIN, &
nelmt_finite,nelmt_trinfinite,nelmt_infinite, &
elmt_finite,elmt_trinfinite,elmt_infinite, &
nnode_finite,nnode_trinfinite,nnode_infinite, &
node_finite,node_trinfinite,node_infinite, &
g_num_finite,g_num_trinfinite,g_num_infinite,myrank
implicit none
! ipass=0: first pass, 1: second pass
integer,intent(in) :: ipass
integer :: i_elmt,ielmt,i_node,inode,i_fin,i_trinf,i_inf,imat,ne
logical,allocatable :: isnode(:)
integer,allocatable :: nmir(:)

if(ipass==0)then
  nelmt_infinite=count(mat_domain(mat_id).ge.ELASTIC_INFDOMAIN)
  nelmt_trinfinite=count(mat_domain(mat_id).eq.ELASTIC_TRINFDOMAIN .or. &
                         mat_domain(mat_id).eq.VISCOELASTIC_TRINFDOMAIN)
  nelmt_finite=nelmt-nelmt_trinfinite-nelmt_infinite

  allocate(elmt_finite(nelmt_finite),elmt_trinfinite(nelmt_trinfinite), &
  elmt_infinite(nelmt_infinite))
  elmt_finite=-1
  elmt_trinfinite=-1
  elmt_infinite=-1
  i_fin=0
  i_trinf=0
  i_inf=0
  do i_elmt=1,nelmt
    imat=mat_id(i_elmt)
    if(mat_domain(imat).lt.ELASTIC_TRINFDOMAIN)then
      i_fin=i_fin+1
      elmt_finite(i_fin)=i_elmt
    elseif(mat_domain(imat).eq.ELASTIC_TRINFDOMAIN .or. &
           mat_domain(imat).eq.VISCOELASTIC_TRINFDOMAIN)then
      i_trinf=i_trinf+1
      elmt_trinfinite(i_trinf)=i_elmt
    else
      i_inf=i_inf+1
      elmt_infinite(i_inf)=i_elmt
      if(myrank==0)then
        if(i_elmt==9134)print*,i_elmt,i_inf
      endif
    endif
  enddo
  if(myrank==0)then
    write(logunit,'(a,i0,a,i0,a,i0)')' nelmt_finite: ',nelmt_finite, &
                                ' nelmt_trinfinite: ',nelmt_trinfinite, &
                                ' nelmt_infinite: ',nelmt_infinite
    flush(logunit)
  endif
  ! number of nodes per element is ngnode
  ne=ngnode
else
  ! number of nodes per element is ngll or nenode
  ne=nenode
endif
allocate(g_num_finite(ne,nelmt_finite),g_num_trinfinite(ne,nelmt_trinfinite), &
g_num_infinite(ne,nelmt_infinite))
g_num_finite=-1
g_num_trinfinite=-1
g_num_infinite=-1

allocate(isnode(nnode))
! classify trinfinite elements
isnode=.false.
do i_elmt=1,nelmt_trinfinite
   ielmt=elmt_trinfinite(i_elmt)
   isnode(g_num(:,ielmt))=.true.
enddo
nnode_trinfinite=count(isnode)
allocate(node_trinfinite(nnode_trinfinite))
allocate(nmir(nnode))
nmir=-1
inode=0
node_trinfinite=-1
do i_node=1,nnode
  if(isnode(i_node))then
    inode=inode+1
    nmir(i_node)=inode
    node_trinfinite(inode)=i_node
  endif
enddo
! store connectivity
do i_elmt=1,nelmt_trinfinite
   ielmt=elmt_trinfinite(i_elmt)
   g_num_trinfinite(:,i_elmt)=nmir(g_num(:,ielmt))
enddo

! classify infinite elements
isnode=.false.
do i_elmt=1,nelmt_infinite
   ielmt=elmt_infinite(i_elmt)
   isnode(g_num(:,ielmt))=.true.
enddo
nnode_infinite=count(isnode)
allocate(node_infinite(nnode_infinite))
nmir=-1
inode=0
node_infinite=-1
do i_node=1,nnode
  if(isnode(i_node))then
    inode=inode+1
    nmir(i_node)=inode
    node_infinite(inode)=i_node
  endif
enddo
! store connectivity
do i_elmt=1,nelmt_infinite
   ielmt=elmt_infinite(i_elmt)
   g_num_infinite(:,i_elmt)=nmir(g_num(:,ielmt))
enddo

! classify finite elements
isnode=.false.
do i_elmt=1,nelmt_finite
   ielmt=elmt_finite(i_elmt)
   isnode(g_num(:,ielmt))=.true.
enddo
nnode_finite=count(isnode)
allocate(node_finite(nnode_finite))
nmir=-1
inode=0
node_finite=-1
do i_node=1,nnode
  if(isnode(i_node))then
    inode=inode+1
    nmir(i_node)=inode
    node_finite(inode)=i_node
  endif
enddo
deallocate(isnode)
! store connectivity
do i_elmt=1,nelmt_finite
   ielmt=elmt_finite(i_elmt)
   g_num_finite(:,i_elmt)=nmir(g_num(:,ielmt))
enddo
if(allocated(nmir))deallocate(nmir)
return
end subroutine classify_finite_infinite_elements
!===============================================================================

! This subroutine computes GLL (along finite directions) and Radau
! (along infinite direction)
! quadrature points and weights for 3D infinite element with a decay 1/r
subroutine shape_function_infiniteGLHEX8ZW(inftype,ngllx,nglly,ngllz,ngll,nip,&
isfaces,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use math_constants,only:ZERO
use gll_library,only:lagrange1dGLL,lagrange1dGLLAS,lagrange1dGEN,              &
lagrange1dGENAS,zwgjd,zwgljd
implicit none
integer,intent(in) :: inftype ! type of integration along infinite direction
integer,intent(in) :: ngllx,nglly,ngllz,ngll,nip
logical,intent(in) :: isfaces(6)
real(kind=kreal),dimension(nip,8),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,8),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1,nipx(ndim)
real(kind=kreal) :: sdir(ndim),xi(ndim)
! gll points and weights
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx,igllpx,igllwx
real(kind=kreal),dimension(nglly) :: gllpy,gllwy,igllpy,igllwy
real(kind=kreal),dimension(ngllz) :: gllpz,gllwz,igllpz,igllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx

! x:1, y:2, z:3 direction for each face
integer,parameter :: idir_face(6)=(/ 2, 1, 2, 1, 3, 3 /)

! Sign of the direction for each face
integer,parameter :: sdir_face(6)=(/ -ONE, ONE, ONE, -ONE, -ONE, ONE /)
integer :: i_dir,i_face,idir
logical :: isface
logical :: isINF(ndim)

isINF=.false.
sdir=one
do i_face=1,6
  isface=isfaces(i_face)
  if(isface)then
    idir=idir_face(i_face)
    isINF(idir)=.true.
    sdir(idir)=sdir_face(i_face)
  endif
enddo

nipx(1)=ngllx
nipx(2)=nglly
nipx(3)=ngllz

! Compute everything in indexed order
! For alpha=beta=0, jacobi polynomial is legendre polynomial
! For ngllx=nglly=ngllz, need to call only once
! Get gll points and weights
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

! Integration points are the GLL points
! Initially set all to GLL quadrature, i.e., standard directions
igllpx=gllpx; igllpy=gllpy; igllpz=gllpz
igllwx=gllwx; igllwy=gllwy; igllwz=gllwz;

! Overwrite GLL points/weights with appropriate quadrature 
! along infinite direction/s
if(inftype==INFINITE_RADAU)then
  ! Radau
  if(isINF(1))call radau_quadrature(ngllx,igllpx,igllwx)
  if(isINF(2))call radau_quadrature(nglly,igllpy,igllwy)
  if(isINF(3))call radau_quadrature(ngllz,igllpz,igllwz)
elseif(inftype==INFINITE_GAUSS)then
  ! Gauss 
  if(isINF(1))call zwgjd(igllpx,igllwx,ngllx,jacobi_alpha,jacobi_beta)
  if(isINF(2))call zwgjd(igllpy,igllwy,nglly,jacobi_alpha,jacobi_beta)
  if(isINF(3))call zwgjd(igllpz,igllwz,ngllz,jacobi_alpha,jacobi_beta)
else
  write(*,*)'ERROR: invalid "inftype":',inftype
  stop
endif

ii=0
do k1=1,nipx(3)
  do j1=1,nipx(2)
    do i1=1,nipx(1)
      ii=ii+1

      ! Integration points
      xi(1)=igllpx(i1)
      xi(2)=igllpy(j1)
      xi(3)=igllpz(k1)

      ! Integration weights
      GLw(ii)=igllwx(i1)*igllwy(j1)*igllwz(k1)      
      
      call lagrange1dGLLAS(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLLAS(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLLAS(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))     

      ! Interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)= &
              lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)= &
              lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)= &
              lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)= &
              lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! Shape functions for HEX8
      ! Compute 1d lagrange polynomials
      call lagrange1dGENAS(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGENAS(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGENAS(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! Consider 3 nodes but compute only at 2 nodes
      do i_dir=1,NDIM
        if(isINF(i_dir))then
          if(sdir(i_dir).lt.ZERO)then
            ! Reverse direction
            call lagrange1d_infiniteZWRAS(3,xi(i_dir), &
                 lagrangeINF_x(i_dir,:),lagrangeINF_dx(i_dir,:))
          else
            ! Standard direction
            call lagrange1d_infiniteZWAS(3,xi(i_dir),    &
                 lagrangeINF_x(i_dir,:),lagrangeINF_dx(i_dir,:))
          endif
        endif
      enddo
      n=0
      do k=1,2
        do j=1,2
          do i=1,2
            n=n+1
            shape_infinite(ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)= &
              lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo
return
end subroutine shape_function_infiniteGLHEX8ZW
!===============================================================================

! This subroutine computes GLL (along finite directions) and Radau
! (along infinite direction)
! quadrature points and weights for 3D
subroutine shape_function_infiniteGLHEX8ZW_GLLRADAU(ngllx,nglly,ngllz,ngll,nip,&
isfaces,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use math_constants,only:ZERO
use gll_library,only:lagrange1dGLL,lagrange1dGLLAS,lagrange1dGEN,              &
lagrange1dGENAS,zwgljd
implicit none
integer,intent(in) :: ngllx,nglly,ngllz,ngll,nip
logical,intent(in) :: isfaces(6)
real(kind=kreal),intent(in) :: gam,a,nd !nd: order of decay
real(kind=kreal),dimension(nip,8),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,8),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1,nipx(ndim)
real(kind=kreal) :: sdir(ndim),xi(ndim)
! gll points and weights
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx,igllpx,igllwx
real(kind=kreal),dimension(nglly) :: gllpy,gllwy,igllpy,igllwy
real(kind=kreal),dimension(ngllz) :: gllpz,gllwz,igllpz,igllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx

! x:1, y:2, z:3 direction for each face
integer,parameter :: idir_face(6)=(/ 2, 1, 2, 1, 3, 3 /)

! Sign of the direction for each face
integer,parameter :: sdir_face(6)=(/ -ONE, ONE, ONE, -ONE, -ONE, ONE /)
integer :: i_dir,i_face,idir
logical :: isface
logical :: isINF(ndim)

isINF=.false.
sdir=one
do i_face=1,6
  isface=isfaces(i_face)
  if(isface)then
    idir=idir_face(i_face)
    isINF(idir)=.true.
    sdir(idir)=sdir_face(i_face)
  endif
enddo
!if(iface==1)then
!  isINF(2)=.true.; ddir=-one
!elseif(iface==2)then
!  isINF(1)=.true.; ddir=one
!elseif(iface==3)then
!  isINF(2)=.true.; ddir=one
!elseif(iface==4)then
!  isINF(1)=.true.; ddir=-one
!elseif(iface==5)then
!  isINF(3)=.true.; ddir=-one
!elseif(iface==6)then
!  isINF(3)=.true.; ddir=one
!endif

nipx(1)=ngllx
nipx(2)=nglly
nipx(3)=ngllz

! Compute everything in indexed order
! For alpha=beta=0, jacobi polynomial is legendre polynomial
! For ngllx=nglly=ngllz, need to call only once
! Get gll points and weights
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

! Integration points are the GLL points
! Initially set all to GLL quadrature, i.e., statndar directions
igllpx=gllpx; igllpy=gllpy; igllpz=gllpz
igllwx=gllwx; igllwy=gllwy; igllwz=gllwz;

! Overwrite GLL points/weights with radau counterpart along infinite direction/s 
if(isINF(1))call radau_quadrature(ngllx,igllpx,igllwx)
if(isINF(2))call radau_quadrature(nglly,igllpy,igllwy)
if(isINF(3))call radau_quadrature(ngllz,igllpz,igllwz)

ii=0
do k1=1,nipx(3)
  do j1=1,nipx(2)
    do i1=1,nipx(1)
      ii=ii+1

      ! Integration points
      xi(1)=igllpx(i1)
      xi(2)=igllpy(j1)
      xi(3)=igllpz(k1)

      ! Integration weights
      GLw(ii)=igllwx(i1)*igllwy(j1)*igllwz(k1)      
      
      call lagrange1dGLLAS(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLLAS(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLLAS(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))     

      ! Interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)= &
              lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)= &
              lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)= &
              lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)= &
              lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! Shape functions for HEX8
      ! Compute 1d lagrange polynomials
      call lagrange1dGENAS(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGENAS(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGENAS(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! Consider 3 nodes but compute only at 2 nodes
      do i_dir=1,NDIM
        if(isINF(i_dir))then
          if(sdir(i_dir).lt.ZERO)then
            ! Reverse direction
            call lagrange1d_infiniteZWRAS(3,xi(i_dir), &
                 lagrangeINF_x(i_dir,:),lagrangeINF_dx(i_dir,:))
          else
            ! Standard direction
            call lagrange1d_infiniteZWAS(3,xi(i_dir),    &
                 lagrangeINF_x(i_dir,:),lagrangeINF_dx(i_dir,:))
          endif
        endif
      enddo
      n=0
      do k=1,2
        do j=1,2
          do i=1,2
            n=n+1
            shape_infinite(ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)= &
              lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo
return
end subroutine shape_function_infiniteGLHEX8ZW_GLLRADAU
!===============================================================================

! This subroutine computes GLL (along finite directions) and Gauss
! (along infinite direction)
! quadrature points and weights for 3D
subroutine shape_function_infiniteGLHEX8ZW_GLLGAUSS(ngllx,nglly,ngllz,ngll,nip,&
isfaces,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use math_constants,only:ZERO
use gll_library,only:lagrange1dGLL,lagrange1dGLLAS,lagrange1dGEN,              &
lagrange1dGENAS,zwgjd,zwgljd
implicit none
integer,intent(in) :: ngllx,nglly,ngllz,ngll,nip
logical,intent(in) :: isfaces(6)
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,8),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,8),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1,nipx(ndim)
real(kind=kreal) :: sdir(ndim),xi(ndim)
! gll points and weights
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx,igllpx,igllwx
real(kind=kreal),dimension(nglly) :: gllpy,gllwy,igllpy,igllwy
real(kind=kreal),dimension(ngllz) :: gllpz,gllwz,igllpz,igllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx

! x:1, y:2, z:3 direction for each face
integer,parameter :: idir_face(6)=(/ 2, 1, 2, 1, 3, 3 /)

! Sign of the direction for each face
integer,parameter :: sdir_face(6)=(/ -ONE, ONE, ONE, -ONE, -ONE, ONE /)
integer :: i_dir,i_face,idir
logical :: isface
logical :: isINF(ndim)

isINF=.false.
sdir=one
do i_face=1,6
  isface=isfaces(i_face)
  if(isface)then
    idir=idir_face(i_face)
    isINF(idir)=.true.
    sdir(idir)=sdir_face(i_face)
  endif
enddo
!if(iface==1)then
!  isINF(2)=.true.; ddir=-one
!elseif(iface==2)then
!  isINF(1)=.true.; ddir=one
!elseif(iface==3)then
!  isINF(2)=.true.; ddir=one
!elseif(iface==4)then
!  isINF(1)=.true.; ddir=-one
!elseif(iface==5)then
!  isINF(3)=.true.; ddir=-one
!elseif(iface==6)then
!  isINF(3)=.true.; ddir=one
!endif

nipx(1)=ngllx
nipx(2)=nglly
nipx(3)=ngllz

! Compute everything in indexed order
! For alpha=beta=0, jacobi polynomial is legendre polynomial
! For ngllx=nglly=ngllz, need to call only once
! Get gll points and weights
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

! Integration points are the GLL points
! Initially set all to GLL quadrature, i.e., statndar directions
igllpx=gllpx; igllpy=gllpy; igllpz=gllpz
igllwx=gllwx; igllwy=gllwy; igllwz=gllwz;

! Overwrite GLL points/weights with Gauss counterpart along infinite direction/s 
if(isINF(1))call zwgjd(igllpx,igllwx,ngllx,jacobi_alpha,jacobi_beta)
if(isINF(2))call zwgjd(igllpy,igllwy,nglly,jacobi_alpha,jacobi_beta)
if(isINF(3))call zwgjd(igllpz,igllwz,ngllz,jacobi_alpha,jacobi_beta)

ii=0
do k1=1,nipx(3)
  do j1=1,nipx(2)
    do i1=1,nipx(1)
      ii=ii+1

      ! Integration points
      xi(1)=igllpx(i1)
      xi(2)=igllpy(j1)
      xi(3)=igllpz(k1)

      ! Integration weights
      GLw(ii)=igllwx(i1)*igllwy(j1)*igllwz(k1)      
      
      call lagrange1dGLLAS(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLLAS(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLLAS(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))     

      ! Interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)= &
              lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)= &
              lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)= &
              lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)= &
              lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! Shape functions for HEX8
      ! Compute 1d lagrange polynomials
      call lagrange1dGENAS(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGENAS(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGENAS(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! Consider 3 nodes but compute only at 2 nodes
      do i_dir=1,NDIM
        if(isINF(i_dir))then
          if(sdir(i_dir).lt.ZERO)then
            ! Reverse direction
            call lagrange1d_infiniteZWRAS(3,xi(i_dir), &
                 lagrangeINF_x(i_dir,:),lagrangeINF_dx(i_dir,:))
          else
            ! Standard direction
            call lagrange1d_infiniteZWAS(3,xi(i_dir),    &
                 lagrangeINF_x(i_dir,:),lagrangeINF_dx(i_dir,:))
          endif
        endif
      enddo
      n=0
      do k=1,2
        do j=1,2
          do i=1,2
            n=n+1
            shape_infinite(ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)= &
              lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)= &
              lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo
return
end subroutine shape_function_infiniteGLHEX8ZW_GLLGAUSS
!===============================================================================

! This subroutine computes Gauss quadrature points and weights for 3D
subroutine shape_function_infiniteGLHEX8ZW_GAUSS(ngllx,nglly,ngllz,ngll,nipx,nip,&
iface,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL,lagrange1dGEN,zwgjd,zwgljd
implicit none
integer,intent(in) :: ngllx,nglly,ngllz,ngll,nipx,nip,iface
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,8),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,8),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1
real(kind=kreal) :: ddir,xi(ndim)
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx
! gll points and weights
real(kind=kreal),dimension(nglly) :: gllpy,gllwy
real(kind=kreal),dimension(ngllz) :: gllpz,gllwz
real(kind=kreal),dimension(nipx) :: igllpx,igllwx
real(kind=kreal),dimension(nipx) :: igllpy,igllwy
real(kind=kreal),dimension(nipx) :: igllpz,igllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx
integer :: iINF

ddir=one
if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

! interpolation points
! compute everything in indexed order
! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz, need to call only once
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

! integration points
! gauss-jacobi or gauss-legendre points and weights
call zwgjd(igllpx,igllwx,nipx,jacobi_alpha,jacobi_beta)
if(nip.ne.8.and.nip.ne.27.and.nip.ne.64)then
  print*,'ERROR: illegal number of Gauss points:',nip,'!'
  stop
endif
igllpy=igllpx; igllpz=igllpx;
igllwy=igllwx; igllwz=igllwx;

ii=0
do k1=1,nipx
  do j1=1,nipx
    do i1=1,nipx
      ii=ii+1

      ! integration points
      xi(1)=igllpx(i1) !xi,   gll_points(1,ii)
      xi(2)=igllpy(j1) !eta,  gll_points(2,ii)
      xi(3)=igllpz(k1) !zeta, gll_points(3,ii)

      ! integration weights
      GLw(ii)=igllwx(i1)*igllwy(j1)*igllwz(k1)      
      
      call lagrange1dGLL(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLL(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLL(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! shape functions for HEX8
      ! compute 1d lagrange polynomials
      call lagrange1dGEN(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGEN(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGEN(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! consider 3 nodes but compute only at 2 nodes
      call lagrange1d_infiniteZW(3,xi(iINF),lagrangeINF_x(iINF,:),         &
      lagrangeINF_dx(iINF,:))      
      n=0
      do k=1,2
        do j=1,2
          do i=1,2
            n=n+1
            shape_infinite(ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)=lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo

return
end subroutine shape_function_infiniteGLHEX8ZW_GAUSS
!===============================================================================

! This subroutine computes GLL (along finite directions) and Radau
! (along infinite direction)
! quadrature points and weights for 3D
subroutine shape_function_infiniteGLHEX8ZW_GLLR(ngllx,nglly,ngllz,ngll,nip,   &
iface,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL,lagrange1dGEN,zwgljd
implicit none
integer,intent(in) :: ngllx,nglly,ngllz,ngll,nip,iface
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,8),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,8),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal!,five=5.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1,nipx(ndim)
real(kind=kreal) :: ddir,xi(ndim) !,eta,zeta
! gll points and weights
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx,igllpx,igllwx
real(kind=kreal),dimension(nglly) :: gllpy,gllwy,igllpy,igllwy
real(kind=kreal),dimension(ngllz) :: gllpz,gllwz,igllpz,igllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx
integer :: iINF

ddir=one
if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

nipx(1)=ngllx
nipx(2)=nglly
nipx(3)=ngllz
! compute everything in indexed order
! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz, need to call only once
! get gll points and weights
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

! integration points are the GLL points
igllpx=gllpx; igllpy=gllpy; igllpz=gllpz
igllwx=gllwx; igllwy=gllwy; igllwz=gllwz;

! overwrite GLL points/weights with radau counterpart along infinite direction 
if(iINF.eq.1)call radau_quadrature(ngllx,igllpx,igllwx)
if(iINF.eq.2)call radau_quadrature(nglly,igllpy,igllwy)
if(iINF.eq.3)call radau_quadrature(ngllz,igllpz,igllwz)

ii=0
do k1=1,nipx(3)
  do j1=1,nipx(2)
    do i1=1,nipx(1)
      ii=ii+1

      ! integration points
      xi(1)=igllpx(i1) !xi,   gll_points(1,ii)
      xi(2)=igllpy(j1) !eta,  gll_points(2,ii)
      xi(3)=igllpz(k1) !zeta, gll_points(3,ii)

      ! integration weights
      GLw(ii)=igllwx(i1)*igllwy(j1)*igllwz(k1)      
      
      call lagrange1dGLL(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLL(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLL(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))     

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! shape functions for HEX8
      ! compute 1d lagrange polynomials
      call lagrange1dGEN(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGEN(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGEN(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! consider 3 nodes but compute only at 2 nodes
      call lagrange1d_infiniteZW(3,xi(iINF),lagrangeINF_x(iINF,:),         &
      lagrangeINF_dx(iINF,:))   
      n=0
      do k=1,2
        do j=1,2
          do i=1,2
            n=n+1
            shape_infinite(ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)=lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo
return
end subroutine shape_function_infiniteGLHEX8ZW_GLLR
!===============================================================================

! This subroutine computes Gauss quadrature points and weights for 3D
subroutine shape_function_infiniteGLHEX8ZW_GQ(ngllx,nglly,ngllz,ngll,nipx,nip,&
iface,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL,lagrange1dGEN,zwgjd,zwgljd
implicit none
integer,intent(in) :: ngllx,nglly,ngllz,ngll,nipx,nip,iface
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,8),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,8),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1
real(kind=kreal) :: ddir,xi(ndim)
real(kind=kreal),dimension(ngllx) :: gllpx,gllwx
! gll points and weights
real(kind=kreal),dimension(nglly) :: gllpy,gllwy
real(kind=kreal),dimension(ngllz) :: gllpz,gllwz
real(kind=kreal),dimension(nipx) :: igllpx,igllwx
real(kind=kreal),dimension(nipx) :: igllpy,igllwy
real(kind=kreal),dimension(nipx) :: igllpz,igllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx
integer :: iINF

ddir=one
if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

! interpolation points
! compute everything in indexed order
! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz, need to call only once
call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

! integration points
! gauss-jacobi or gauss-legendre points and weights
call zwgjd(igllpx,igllwx,nipx,jacobi_alpha,jacobi_beta)
if(nip.ne.8.and.nip.ne.27.and.nip.ne.64)then
  print*,'ERROR: illegal number of Gauss points:',nip,'!'
  stop
endif
igllpy=igllpx; igllpz=igllpx;
igllwy=igllwx; igllwz=igllwx;

ii=0
do k1=1,nipx
  do j1=1,nipx
    do i1=1,nipx
      ii=ii+1

      ! integration points
      xi(1)=igllpx(i1) !xi,   gll_points(1,ii)
      xi(2)=igllpy(j1) !eta,  gll_points(2,ii)
      xi(3)=igllpz(k1) !zeta, gll_points(3,ii)

      ! integration weights
      GLw(ii)=igllwx(i1)*igllwy(j1)*igllwz(k1)      
      
      call lagrange1dGLL(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLL(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLL(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! shape functions for HEX8
      ! compute 1d lagrange polynomials
      call lagrange1dGEN(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGEN(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGEN(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! consider 3 nodes but compute only at 2 nodes
      call lagrange1d_infiniteZW(3,xi(iINF),lagrangeINF_x(iINF,:),         &
      lagrangeINF_dx(iINF,:))      
      n=0
      do k=1,2
        do j=1,2
          do i=1,2
            n=n+1
            shape_infinite(ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)=lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo

return
end subroutine shape_function_infiniteGLHEX8ZW_GQ
!===============================================================================

! This subroutine computes GLL quadrature points and weights for 3D
subroutine shape_function_infiniteGLHEX8MO(ngllx,nglly,ngllz,ngll,nipx,nip,   &
iface,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL,lagrange1dGEN,zwgjd
implicit none
integer,intent(in) :: ngllx,nglly,ngllz,ngll,nipx,nip,iface
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,8),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,8),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1
real(kind=kreal) :: ddir,xi(ndim)
! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpx,gllwx
real(kind=kreal),dimension(nipx) :: gllpy,gllwy
real(kind=kreal),dimension(nipx) :: gllpz,gllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx
integer :: iINF !,ngllxINF(ndim)

ddir=one
if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

! compute everything in indexed order
! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz, need to call only once
! gauss-jacobi or gauss-legendre points and weights
call zwgjd(gllpx,gllwx,nipx,jacobi_alpha,jacobi_beta)
! for an infinite element we use Gauss-Legendre quadrature
!if(nip==8)then
!  gllpx(1)=-one/sqrt(3.0_kreal); gllpx(2)=-gllpx(1)
!  gllwx(1)=one; gllwx(2)=one;  
!elseif(nip==27)then
!  gllpx(1)=-sqrt(3.0_kreal/five); gllpx(2)=0.0_kreal; gllpx(3)=-gllpx(1)
!  gllwx(1)=five/9.0_kreal; gllwx(2)=8.0_kreal/9.0_kreal; gllwx(3)=gllwx(1);
!else
if(nip.ne.8.and.nip.ne.27.and.nip.ne.64)then
  print*,'ERROR: illegal number of Gauss points:',nip,'!'
  stop
endif
gllpy=gllpx; gllpz=gllpx;
gllwy=gllwx; gllwz=gllwx;

ii=0
do k1=1,nipx
  do j1=1,nipx
    do i1=1,nipx
      ii=ii+1
      xi(1)=gllpx(i1) !xi,   gll_points(1,ii)
      xi(2)=gllpy(j1) !eta,  gll_points(2,ii)
      xi(3)=gllpz(k1) !zeta, gll_points(3,ii)

      GLw(ii)=gllwx(i1)*gllwy(j1)*gllwz(k1)      
      
      call lagrange1dGLL(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLL(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLL(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))     

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! shape functions for HEX8
      ! compute 1d lagrange polynomials
      call lagrange1dGEN(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGEN(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGEN(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! consider 3 nodes but compute only at 2 nodes
      call lagrange1d_infiniteMO(3,nd,xi(iINF),lagrangeINF_x(iINF,:),lagrangeINF_dx(iINF,:))      
      n=0
      do k=1,2
        do j=1,2
          do i=1,2
            n=n+1
            shape_infinite(ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)=lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo

return
end subroutine shape_function_infiniteGLHEX8MO
!===============================================================================

! This subroutine extracts the nodes for HEX8 of the finite region of an       
! infinite element
subroutine get_gnodinfHEX8M(ndim,ngllx,nglly,ngllz,ngllinf,isfaces,gnodinf)
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngllinf
logical,intent(in) :: isfaces(6)
integer,intent(out) :: gnodinf(ngllinf)
integer :: i,j,k,inum
integer :: inc(ndim),ngllxINF0(ndim),ngllxINF(ndim),iINF
real(kind=kreal) :: ddir
real(kind=kreal) :: sdir(ndim)
real(kind=kreal),parameter :: one=1.0_kreal

! x:1, y:2, z:3 direction for each face
integer,parameter :: idir_face(6)=(/ 2, 1, 2, 1, 3, 3 /)

! Sign of the direction for each face
integer,parameter :: sdir_face(6)=(/ -ONE, ONE, ONE, -ONE, -ONE, ONE /)
integer :: i_dir,i_face,idir
logical :: isface
logical :: isINF(ndim)

! Element nodes must be reordered such that the infinite faces are always in the
! far ends
isINF=.false.
sdir=one
do i_face=1,6
  isface=isfaces(i_face)
  if(isface)then
    idir=idir_face(i_face)
    isINF(idir)=.true.
    sdir(idir)=sdir_face(i_face)
  endif
enddo

! Initialize ngllINF indices
ngllxINF0=1
ngllxINF(1)=ngllx
ngllxINF(2)=nglly
ngllxINF(3)=ngllz

! Modify and overwrite for infinite direction/s
do i_dir=1,NDIM
  if(isINF(i_dir))then
    if(sdir(i_dir)<0)then
      write(*,*)'ERROR: infinite face must be in far ends!'
      ngllxINF0(i_dir)=2
    else
      ngllxINF(i_dir)=ngllxINF(i_dir)-1
    endif
  endif
enddo

!if(iface==1)then
!  iINF=2; ddir=-one
!elseif(iface==2)then
!  iINF=1; ddir=one
!elseif(iface==3)then
!  iINF=2; ddir=one
!elseif(iface==4)then
!  iINF=1; ddir=-one
!elseif(iface==5)then
!  iINF=3; ddir=-one
!elseif(iface==6)then
!  iINF=3; ddir=one
!endif
!
!if(ddir<0)then
!  ngllxINF0(iINF)=2
!else
!  ngllxINF(iINF)=ngllxINF(iINF)-1
!endif

! extract only the corner nodes
inc=ngllxINF-ngllxINF0
inum=0
do k=ngllxINF0(3),ngllxINF(3),inc(3)
  do j=ngllxINF0(2),ngllxINF(2),inc(2)
    do i=ngllxINF0(1),ngllxINF(1),inc(1)
      inum=inum+1
      gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
    enddo
  enddo
enddo
end subroutine get_gnodinfHEX8M
!===============================================================================

! This subroutine extracts the nodes for HEX8 of the finite region of an       
! infinite element
subroutine get_gnodinfHEX8(ndim,ngllx,nglly,ngllz,ngllinf,iface,gnodinf)
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngllinf,iface
integer,intent(out) :: gnodinf(ngllinf)
integer :: i,j,k,inum
integer :: inc(ndim),ngllxINF0(ndim),ngllxINF(ndim),iINF
real(kind=kreal) :: ddir
real(kind=kreal),parameter :: one=1.0_kreal

if(iface.lt.1.or.iface.gt.6)then
  write(*,*)'ERROR: illegal outer face ID:',iface
  stop
endif

! initialize ngllINF indices
ngllxINF0=1
ngllxINF(1)=ngllx
ngllxINF(2)=nglly
ngllxINF(3)=ngllz

if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

if(ddir<0)then
  ngllxINF0(iINF)=2
else
  ngllxINF(iINF)=ngllxINF(iINF)-1
endif

! extract only the corner nodes
inc=ngllxINF-ngllxINF0
inum=0
do k=ngllxINF0(3),ngllxINF(3),inc(3)
  do j=ngllxINF0(2),ngllxINF(2),inc(2)
    do i=ngllxINF0(1),ngllxINF(1),inc(1)
      inum=inum+1
      gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
    enddo
  enddo
enddo
end subroutine get_gnodinfHEX8
!===============================================================================

! This subroutine computes GLL quadrature points and weights for 3D
subroutine shape_function_infiniteGLHEXN2DF(ndim,ngllx,nglly,ngllz,ngll,      &
ngllinf,nipx,nip,iface,nd,gam,a,rgll,coordinf,shape_infinite,dshape_infinite, &
lagrange_gl,dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL,zwgjd
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,ngllinf,nipx,nip,iface
real(kind=kreal),intent(in) :: gam,a,rgll(ngll),coordinf(ngllinf),nd !nd: order
real(kind=kreal),dimension(nip,ngllinf),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,ngllinf),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,    &
one=1.0_kreal,five=5.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1
real(kind=kreal) :: ddir,xi(ndim)
! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpx,gllwx
real(kind=kreal),dimension(nipx) :: gllpy,gllwy
real(kind=kreal),dimension(nipx) :: gllpz,gllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,ngllx) :: lagrangeINF_x,lagrangeINF_dx
integer :: inc(ndim),ngllxINF0(ndim),ngllxINF(ndim),iINF,ind0(3),ind1(3)

if(nip.ne.8.and.nip.ne.27.and.nip.ne.64)then
  print*,'ERROR: illegal number of Gauss points:',nip,'!'
  stop
endif

ngllxINF0=1
ngllxINF(1)=ngllx; ngllxINF(2)=nglly; ngllxINF(3)=ngllz
inc=1

ddir=one
if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

if(ddir<0)then
  ngllxINF0(iINF)=2
else
  ngllxINF(iINF)=ngllxINF(iINF)-1
endif
inc(iINF)=ngllxINF(iINF)-ngllxINF0(iINF) ! extract only two nodes in infinity
ind0=1
ind1(1)=ngllx; ind1(2)=nglly; ind1(3)=ngllz
ind1(iINF)=2
!direction

! compute everything in indexed order
! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz, need to call only once
! gauss-jacobi or gauss-legendre points and weights
call zwgjd(gllpx,gllwx,nipx,jacobi_alpha,jacobi_beta)
! for an infinite element we use Gauss-Legendre quadrature
!if(nip==8)then
!  gllpx(1)=-one/sqrt(3.0_kreal); gllpx(2)=-gllpx(1)
!  gllwx(1)=one; gllwx(2)=one;  
!elseif(nip==27)then
!  gllpx(1)=-sqrt(3.0_kreal/five); gllpx(2)=0.0_kreal; gllpx(3)=-gllpx(1)
!  gllwx(1)=five/9.0_kreal; gllwx(2)=8.0_kreal/9.0_kreal; gllwx(3)=gllwx(1);
!else

gllpy=gllpx; gllpz=gllpx;
gllwy=gllwx; gllwz=gllwx;

ii=0
do k1=1,nipx
  do j1=1,nipx
    do i1=1,nipx
      ii=ii+1
      xi(1)=gllpx(i1) !xi,   gll_points(1,ii)
      xi(2)=gllpy(j1) !eta,  gll_points(2,ii)
      xi(3)=gllpz(k1) !zeta, gll_points(3,ii)

      GLw(ii)=gllwx(i1)*gllwy(j1)*gllwz(k1)      
      
      call lagrange1dGLL(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLL(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLL(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))     

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! shape functions for HEXN2
      ! compute 1d lagrange polynomials
      call lagrange1dGLL(ngllx,gllpx,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGLL(nglly,gllpy,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGLL(ngllz,gllpz,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! consider 3 nodes but compute only at 2 nodes
      call lagrange_infinite(3,nd,xi(iINF),ddir,gam,lagrangeINF_x(iINF,1:3),lagrangeINF_dx(iINF,1:3))      
      n=0
      do k=ind0(3),ind1(3) !ngllxINF0(3),ngllxINF(3),inc(3)
        do j=ind0(2),ind1(2) !ngllxINF0(2),ngllxINF(2),inc(2)
          do i=ind0(1),ind1(1) !ngllxINF0(1),ngllxINF(1),inc(1)
            n=n+1
            shape_infinite(ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)=lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo

!r=matmul(matmul(shape_infinite,coordinf),matmul(shape_infinite,coordinf))
!dr_dxi=matmul(shape_infinite,coordinf)*matmul(dshape_infinite,xoordinf)/r
!dM_dxi=dlagrange_gl*ri/r-(ri/r^3)*dr_dxi
return
end subroutine shape_function_infiniteGLHEXN2DF
!===============================================================================

! This subroutine computes GLL quadrature points and weights for 3D
subroutine shape_function_infiniteGLHEXN2(ndim,ngllx,nglly,ngllz,ngll,ngllinf,&
nipx,nip,iface,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,           &
dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL,lagrange1dGLLAS,zwgjd
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,ngllinf,nipx,nip,iface
!of decay
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,ngllinf),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,ngllinf),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal,five=5.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1
real(kind=kreal) :: ddir,xi(ndim)
! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpx,gllwx
real(kind=kreal),dimension(nipx) :: gllpy,gllwy
real(kind=kreal),dimension(nipx) :: gllpz,gllwz
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(ndim,ngllx) :: lagrangeINF_x,lagrangeINF_dx
integer :: inc(ndim),ngllxINF0(ndim),ngllxINF(ndim),iINF,ind0(3),ind1(3)
if(nip.ne.8.and.nip.ne.27.and.nip.ne.64)then
  print*,'ERROR: illegal number of Gauss points:',nip,'!'
  stop
endif

ngllxINF0=1
ngllxINF(1)=ngllx;  ngllxINF(2)=nglly;  ngllxINF(3)=ngllz
inc=1

ddir=one
if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

if(ddir<0)then
  ngllxINF0(iINF)=2
else
  ngllxINF(iINF)=ngllxINF(iINF)-1
endif
!ngllxINF(iINF)=ngllxINF(iINF)-1
inc(iINF)=ngllxINF(iINF)-ngllxINF0(iINF) ! extract only two nodes in infinity
ind0=1
ind1(1)=ngllx; ind1(2)=nglly; ind1(3)=ngllz
ind1(iINF)=2
!direction

! compute everything in indexed order
! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz, need to call only once
! gauss-jacobi or gauss-legendre points and weights
call zwgjd(gllpx,gllwx,nipx,jacobi_alpha,jacobi_beta)
! for an infinite element we use Gauss-Legendre quadrature
!if(nip==8)then
!  gllpx(1)=-one/sqrt(3.0_kreal); gllpx(2)=-gllpx(1)
!  gllwx(1)=one; gllwx(2)=one;  
!elseif(nip==27)then
!  gllpx(1)=-sqrt(3.0_kreal/five); gllpx(2)=0.0_kreal; gllpx(3)=-gllpx(1)
!  gllwx(1)=five/9.0_kreal; gllwx(2)=8.0_kreal/9.0_kreal; gllwx(3)=gllwx(1);
!else

gllpy=gllpx; gllpz=gllpx;
gllwy=gllwx; gllwz=gllwx;

ii=0
do k1=1,nipx
  do j1=1,nipx
    do i1=1,nipx
      ii=ii+1
      !do ii=1,ngll ! ngllx*nglly*ngllz
      xi(1)=gllpx(i1) !xi,   gll_points(1,ii)
      xi(2)=gllpy(j1) !eta,  gll_points(2,ii)
      xi(3)=gllpz(k1) !zeta, gll_points(3,ii)

      !xi(iINF)=ddir*xi(iINF)

      GLw(ii)=gllwx(i1)*gllwy(j1)*gllwz(k1)      
      
      call lagrange1dGLLAS(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLLAS(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLLAS(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))     

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! shape functions for HEXN2
      ! compute 1d lagrange polynomials
      call lagrange1dGLL(ngllx,gllpx,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
      call lagrange1dGLL(nglly,gllpy,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
      call lagrange1dGLL(ngllz,gllpz,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
      ! consider 3 nodes but compute only at 2 nodes
      call lagrange_infinite(3,nd,xi(iINF),ddir,gam,lagrangeINF_x(iINF,1:3),lagrangeINF_dx(iINF,1:3))      
      n=0
      do k=ind0(3),ind1(3) !ngllxINF0(3),ngllxINF(3),inc(3)
        do j=ind0(2),ind1(2) !ngllxINF0(2),ngllxINF(2),inc(2)
          do i=ind0(1),ind1(1) !ngllxINF0(1),ngllxINF(1),inc(1)
            n=n+1
            shape_infinite(ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(1,ii,n)=lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(2,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
            dshape_infinite(3,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
          enddo
        enddo
      enddo

    enddo
  enddo
enddo

return
end subroutine shape_function_infiniteGLHEXN2
!===============================================================================

! This subroutine extracts the nodes for HEXN2 of the finite region of an
! infinite element
subroutine get_gnodinfHEXN2(ndim,ngllx,nglly,ngllz,ngllinf,iface,gnodinf)
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngllinf,iface
integer,intent(out) :: gnodinf(ngllinf)
integer :: i,j,k,inum
integer :: inc(ndim),ngllxINF0(ndim),ngllxINF(ndim),iINF
real(kind=kreal) :: ddir
real(kind=kreal),parameter :: one=1.0_kreal

if(iface.lt.1.or.iface.gt.6)then
  write(*,*)'ERROR: illegal outer face ID:',iface
  stop
endif

! initialize ngllINF indices
ngllxINF0=1
ngllxINF(1)=ngllx; ngllxINF(2)=nglly; ngllxINF(3)=ngllz
inc=1

if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

if(ddir<0)then
  ngllxINF0(iINF)=2
else
  ngllxINF(iINF)=ngllxINF(iINF)-1
endif

! extract mapping nodes
! extract only two nodes along infinite direction first and last nodes of the
! finite region
inc(iINF)=ngllxINF(iINF)-ngllxINF0(iINF)
inum=0
do k=ngllxINF0(3),ngllxINF(3),inc(3)
  do j=ngllxINF0(2),ngllxINF(2),inc(2)
    do i=ngllxINF0(1),ngllxINF(1),inc(1)
      inum=inum+1
      gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
    enddo
  enddo
enddo
end subroutine get_gnodinfHEXN2
!===============================================================================

! This subroutine computes GLL quadrature points and weights for 3D
subroutine shape_function_infiniteGL1(ndim,ngllx,nglly,ngllz,ngll,ngllinf,    &
nipx,nip,iface,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,           &
dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,ngllinf,nipx,nip,iface
!integer,parameter :: ngllinf=ngll-nglly*ngllz
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,ngllinf),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,ngllinf),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)
real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal,five=5.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1
real(kind=kreal) :: ddir,xi(ndim) !,eta,zeta
real(kind=kreal),dimension(nipx) :: gllpx,gllwx ! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpy,gllwy ! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpz,gllwz ! gll points and weights
real(kind=kreal),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
integer :: ngllxINF(ndim),iINF

ngllxINF(1)=ngllx
ngllxINF(2)=nglly
ngllxINF(3)=ngllz

ddir=one
if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

ngllxINF(iINF)=ngllxINF(iINF)-1

! compute everything in indexed order
! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz, need to call only once

! for an infinite element we use Gauss-Legendre quadrature
if(nip==8)then
  gllpx(1)=-one/sqrt(3.0_kreal); gllpx(2)=-gllpx(1)
  gllwx(1)=one; gllwx(2)=one;  
elseif(nip==27)then
  gllpx(1)=-sqrt(3.0_kreal/five); gllpx(2)=0.0_kreal; gllpx(3)=-gllpx(1)
  gllwx(1)=five/9.0_kreal; gllwx(2)=8.0_kreal/9.0_kreal; gllwx(3)=gllwx(1);
else
  print*,'ERROR: illegal number of Gauss points:',nip,'!'
  stop
endif
gllpy=gllpx; gllpz=gllpx;

gllwy=gllwx; gllwz=gllwx;

ii=0
do k1=1,nipx
  do j1=1,nipx
    do i1=1,nipx
      ii=ii+1
      !do ii=1,ngll ! ngllx*nglly*ngllz
      xi(1)=gllpx(i1) !xi,   gll_points(1,ii)
      xi(2)=gllpy(j1) !eta,  gll_points(2,ii)
      xi(3)=gllpz(k1) !zeta, gll_points(3,ii)

      xi(iINF)=ddir*xi(iINF)

      GLw(ii)=gllwx(i1)*gllwy(j1)*gllwz(k1)
      
      
      call lagrange1dGLL(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
      call lagrange1dGLL(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
      call lagrange1dGLL(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))     

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      ! shape functions

      ! compute 1d lagrange polynomials
      call lagrange_infinite(ngllx,nd,xi(iINF),ddir,gam,lagrange_x(iINF,:),lagrange_dx(iINF,:))
      n=0
      do k=1,ngllxINF(3)
        do j=1,ngllxINF(2)
          do i=1,ngllxINF(1)
            n=n+1
            shape_infinite(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dshape_infinite(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
            dshape_infinite(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
            dshape_infinite(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
          enddo
        enddo
      enddo

      !enddo
    enddo
  enddo
enddo

return
end subroutine shape_function_infiniteGL1
!===============================================================================

! This subroutine extracts the nodes of the finite region of an infinite
! element
subroutine get_gnodinf(ndim,ngllx,nglly,ngllz,ngllinf,iface,gnodinf)
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngllinf,iface
integer,intent(out) :: gnodinf(ngllinf)
integer :: i,j,k,inum
integer :: ngllxINF0(ndim),ngllxINF(ndim),iINF
real(kind=kreal) :: ddir
real(kind=kreal),parameter :: one=1.0_kreal

! initialize ngllINF indices
ngllxINF0=1
ngllxINF(1)=ngllx
ngllxINF(2)=nglly
ngllxINF(3)=ngllz

if(iface==1)then
  iINF=2; ddir=-one
elseif(iface==2)then
  iINF=1; ddir=one
elseif(iface==3)then
  iINF=2; ddir=one
elseif(iface==4)then
  iINF=1; ddir=-one
elseif(iface==5)then
  iINF=3; ddir=-one
elseif(iface==6)then
  iINF=3; ddir=one
endif

if(ddir<0)then
  ngllxINF0(iINF)=2
else
  ngllxINF(iINF)=ngllxINF(iINF)-1
endif

if(iface.lt.1.or.iface.gt.6)then
  write(*,*)'ERROR: illegal outer face ID:',iface
  stop
endif

inum=0
do k=ngllxINF0(3),ngllxINF(3)
    do j=ngllxINF0(2),ngllxINF(2)
      do i=ngllxINF0(1),ngllxINF(1)
        inum=inum+1
        gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
      enddo
    enddo
  enddo
end subroutine get_gnodinf
!===============================================================================

! This subroutine computes GLL quadrature points and weights for 3D
subroutine shape_function_infiniteGLCOR(ndim,ngllx,nglly,ngllz,ngll,ngllinf,  &
nipx,nip,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,ngllinf,nipx,nip
!integer,parameter :: ngllinf=ngll-nglly*ngllz
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,ngllinf),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,ngllinf),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal,   &
one=1.0_kreal,five=5.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1
real(kind=kreal) :: xi,eta,zeta
real(kind=kreal),dimension(nipx) :: gllpx,gllwx ! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpy,gllwy ! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpz,gllwz ! gll points and weights
real(kind=kreal),dimension(ngllx-1) :: lagrangeINF_y,lagrangeINF_dy
real(kind=kreal),dimension(ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(nglly) :: lagrange_y,lagrange_dy
real(kind=kreal),dimension(ngllz) :: lagrange_z,lagrange_dz

! compute everything in indexed order

! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz=ngll, need to call only once
if(nip==8)then
  gllpx(1)=-one/sqrt(3.0_kreal); gllpx(2)=-gllpx(1)
  gllwx(1)=one; gllwx(2)=one;  
elseif(nip==27)then
  gllpx(1)=-sqrt(3.0_kreal/five); gllpx(2)=0.0_kreal; gllpx(3)=-gllpx(1)
  gllwx(1)=five/9.0_kreal; gllwx(2)=8.0_kreal/9.0_kreal; gllwx(3)=gllwx(1);
else
  print*,'ERROR: illegal number of Gauss points:',nip,'!'
  stop
endif
gllpy=gllpx; gllpz=gllpx;

gllwy=gllwx; gllwz=gllwx;

ii=0
do k1=1,nipx
  do j1=1,nipx
    do i1=1,nipx
      ii=ii+1
      !do ii=1,ngll ! ngllx*nglly*ngllz
      xi=gllpx(i1) !gll_points(1,ii)
      eta=gllpy(j1) !gll_points(2,ii)
      zeta=gllpz(k1) !gll_points(3,ii)

      GLw(ii)=gllwx(i1)*gllwy(j1)*gllwz(k1)

      ! compute 1d lagrange polynomials
      call lagrange_infinite(ngllz,nd,eta,one,gam,lagrangeINF_y,lagrangeINF_dy)
      call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)
      call lagrange1dGLL(nglly,gllpy,eta,lagrange_y,lagrange_dy)
      call lagrange1dGLL(ngllz,gllpz,zeta,lagrange_z,lagrange_dz)

      ! shape functions
      n=0
      do k=1,ngllz
        do j=1,nglly-1
          do i=1,ngllx
            n=n+1
            shape_infinite(ii,n)=lagrange_x(i)*lagrangeINF_y(j)*lagrange_z(k)
            dshape_infinite(1,ii,n)=lagrange_dx(i)*lagrangeINF_y(j)*lagrange_z(k)
            dshape_infinite(2,ii,n)=lagrange_x(i)*lagrangeINF_dy(j)*lagrange_z(k)
            dshape_infinite(3,ii,n)=lagrange_x(i)*lagrangeINF_y(j)*lagrange_dz(k)
          enddo
        enddo
      enddo

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(i)*lagrange_y(j)*lagrange_z(k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(i)*lagrange_y(j)*lagrange_z(k)
            dlagrange_gl(2,ii,n)=lagrange_x(i)*lagrange_dy(j)*lagrange_z(k)
            dlagrange_gl(3,ii,n)=lagrange_x(i)*lagrange_y(j)*lagrange_dz(k)
          enddo
        enddo
      enddo

      !enddo
    enddo
  enddo
enddo

return
end subroutine shape_function_infiniteGLCOR
!===============================================================================

! This subroutine computes GLL quadrature points and weights for 3D
subroutine shape_function_infiniteGL(ndim,ngllx,nglly,ngllz,ngll,ngllinf,nipx,&
nip,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
use gll_library,only:lagrange1dGLL
implicit none
integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,ngllinf,nipx,nip
!integer,parameter :: ngllinf=ngll-nglly*ngllz
real(kind=kreal),intent(in) :: gam,a,nd !nd: order
real(kind=kreal),dimension(nip,ngllinf),intent(out) :: shape_infinite
real(kind=kreal),dimension(ndim,nip,ngllinf),intent(out) :: dshape_infinite
real(kind=kreal),dimension(nip,ngll),intent(out) :: lagrange_gl
real(kind=kreal),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
real(kind=kreal),intent(out) :: GLw(nip)

real(kind=kreal),parameter :: jacobi_alpha=0.0_kreal,jacobi_beta=0.0_kreal, &
one=1.0_kreal,five=5.0_kreal
integer :: i,ii,j,k,n,i1,j1,k1
real(kind=kreal) :: xi,eta,zeta
real(kind=kreal),dimension(nipx) :: gllpx,gllwx ! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpy,gllwy ! gll points and weights
real(kind=kreal),dimension(nipx) :: gllpz,gllwz ! gll points and weights
real(kind=kreal),dimension(ngllx-1) :: lagrangeINF_x,lagrangeINF_dx
real(kind=kreal),dimension(ngllx) :: lagrange_x,lagrange_dx
real(kind=kreal),dimension(nglly) :: lagrange_y,lagrange_dy
real(kind=kreal),dimension(ngllz) :: lagrange_z,lagrange_dz

! compute everything in indexed order

! get gll points
! for alpha=beta=0, jacobi polynomial is legendre polynomial
! for ngllx=nglly=ngllz=ngll, need to call only once
!print*,nip,nipx; stop
if(nip==8)then
  gllpx(1)=-one/sqrt(3.0_kreal); gllpx(2)=-gllpx(1)
  gllwx(1)=one; gllwx(2)=one;  
elseif(nip==27)then
  gllpx(1)=-sqrt(3.0_kreal/five); gllpx(2)=0.0_kreal; gllpx(3)=-gllpx(1)
  gllwx(1)=five/9.0_kreal; gllwx(2)=8.0_kreal/9.0_kreal; gllwx(3)=gllwx(1);
else
  print*,'ERROR: illegal number of Gauss points:',nip,'!'
  stop
endif
gllpy=gllpx; gllpz=gllpx;

gllwy=gllwx; gllwz=gllwx;

ii=0
do k1=1,nipx
  do j1=1,nipx
    do i1=1,nipx
      ii=ii+1
      !do ii=1,ngll ! ngllx*nglly*ngllz
      xi=gllpx(i1) !gll_points(1,ii)
      eta=gllpy(j1) !gll_points(2,ii)
      zeta=gllpz(k1) !gll_points(3,ii)

      GLw(ii)=gllwx(i1)*gllwy(j1)*gllwz(k1)

      ! compute 1d lagrange polynomials
      call lagrange_infinite(ngllx,nd,xi,one,gam,lagrangeINF_x,lagrangeINF_dx)
      call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)
      call lagrange1dGLL(nglly,gllpy,eta,lagrange_y,lagrange_dy)
      call lagrange1dGLL(ngllz,gllpz,zeta,lagrange_z,lagrange_dz)

      ! shape functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx-1
            n=n+1
            shape_infinite(ii,n)=lagrangeINF_x(i)*lagrange_y(j)*lagrange_z(k)
            dshape_infinite(1,ii,n)=lagrangeINF_dx(i)*lagrange_y(j)*lagrange_z(k)
            dshape_infinite(2,ii,n)=lagrangeINF_x(i)*lagrange_dy(j)*lagrange_z(k)
            dshape_infinite(3,ii,n)=lagrangeINF_x(i)*lagrange_y(j)*lagrange_dz(k)
          enddo
        enddo
      enddo

      ! interpolation functions
      n=0
      do k=1,ngllz
        do j=1,nglly
          do i=1,ngllx
            n=n+1
            lagrange_gl(ii,n)=lagrange_x(i)*lagrange_y(j)*lagrange_z(k)        
            dlagrange_gl(1,ii,n)=lagrange_dx(i)*lagrange_y(j)*lagrange_z(k)
            dlagrange_gl(2,ii,n)=lagrange_x(i)*lagrange_dy(j)*lagrange_z(k)
            dlagrange_gl(3,ii,n)=lagrange_x(i)*lagrange_y(j)*lagrange_dz(k)
          enddo
        enddo
      enddo

      !enddo
    enddo
  enddo
enddo

return
end subroutine shape_function_infiniteGL
!===============================================================================

! This subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.
subroutine lagrange1d_infiniteMO(nenod,nd,xi,phi,dphi_dxi)
implicit none
! number of nodes in an 1d element
integer,intent(in) :: nenod
! xi: point where to calculate lagrange function and its derivative
real(kind=kreal),intent(in) :: nd,xi
real(kind=kreal),dimension(nenod-1),intent(out) :: phi,dphi_dxi
real(kind=kreal) :: fac !dx
real(kind=kreal),parameter :: one=1.0_kreal,two=2.0_kreal

if(nenod.ne.3)then
  write(*,*)'ERROR: infinite element is currently implemented only for 3 nodes!'
  stop
endif
!! compute natural coordnates
!dx=2.0_kreal/real((nenod-1),kreal)! length = 2.0 as xi is taken -1 to +1
!do i=1,nenod
!  ! coordinates when origin is in the left
!  xii(i)=real((i-1),kreal)*dx
!enddo

!! origin is tranformed to mid point
!xii=xii-1.0_kreal

fac=one/(one-xi)

phi(1)=-two*xi*fac
phi(2)=one-phi(1)

dphi_dxi(1)=-two*fac*fac
dphi_dxi(2)=two*fac*fac

return
end subroutine lagrange1d_infiniteMO
!===============================================================================

! This subroutine computes the 1d infinite element shape functions and their
! derivatives at a given point xi.
! X3=\infty when \xi=1.0
subroutine lagrange1d_infiniteZW(nenod,xi,phi,dphi_dxi)
implicit none
integer,intent(in) :: nenod ! number of nodes in an 1d element
! xi: point where to calculate lagrange function and
real(kind=kreal),intent(in) :: xi
!its derivative
real(kind=kreal),dimension(nenod-1),intent(out) :: phi,dphi_dxi
real(kind=kreal) :: fac !dx
real(kind=kreal),parameter :: one=1.0_kreal

if(nenod.ne.3)then
  write(*,*)'ERROR: infinite element is currently implemented only for 3 nodes!'
  stop
endif

fac=one/(one-xi)

phi(1)=-xi*fac
phi(2)=fac

dphi_dxi(1)=-fac*fac
dphi_dxi(2)=fac*fac

return
end subroutine lagrange1d_infiniteZW
!===============================================================================

! This subroutine computes the 1d infinite element shape functions and their
! derivatives at a given point xi.
! X3=\infty when \xi=1.0
subroutine lagrange1d_infiniteZWAS(nenod,xi,phi,dphi_dxi)
implicit none
integer,intent(in) :: nenod ! number of nodes in an 1d element
! xi: point where to calculate lagrange function and
real(kind=kreal),intent(in) :: xi
!its derivative
real(kind=kreal),dimension(:),intent(out) :: phi,dphi_dxi

real(kind=kreal) :: fac !dx
real(kind=kreal),parameter :: one=1.0_kreal

! Check ubound
if( ubound(phi,1).ne.nenod-1 .or. &
    ubound(dphi_dxi,1).ne.nenod-1)then
  write(*,*)'ERROR: array sizes do not match!'
  stop
endif

if(nenod.ne.3)then
  write(*,*)'ERROR: infinite element is currently implemented only for 3 nodes!'
  stop
endif

fac=one/(one-xi)

phi(1)=-xi*fac ! M0
phi(2)=fac     ! M2

dphi_dxi(1)=-fac*fac
dphi_dxi(2)=fac*fac

return
end subroutine lagrange1d_infiniteZWAS
!===============================================================================

! This subroutine computes the 1d infinite element shape functions and their
! derivatives at a given point xi for the reverse direction, i.e.,
! X3=\infty when \xi=-1.0
subroutine lagrange1d_infiniteZWRAS(nenod,xi,phi,dphi_dxi)
implicit none
integer,intent(in) :: nenod ! number of nodes in an 1d element
! xi: point where to calculate lagrange function and
real(kind=kreal),intent(in) :: xi
!its derivative
real(kind=kreal),dimension(:),intent(out) :: phi,dphi_dxi

real(kind=kreal) :: fac !dx
real(kind=kreal),parameter :: one=1.0_kreal

! Check ubound
if( ubound(phi,1).ne.nenod-1 .or. &
    ubound(dphi_dxi,1).ne.nenod-1)then
  write(*,*)'ERROR: array sizes do not match!'
  stop
endif

if(nenod.ne.3)then
  write(*,*)'ERROR: infinite element is currently implemented only for 3 nodes!'
  stop
endif

fac=one/(one+xi)

phi(1)=fac     ! M2
phi(2)=xi*fac  ! M0

dphi_dxi(1)=-fac*fac
dphi_dxi(2)=fac*fac

return
end subroutine lagrange1d_infiniteZWRAS
!===============================================================================

! This subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.
subroutine lagrange_infinite(nenod,nd,xi,ddir,gam,phi,dphi_dxi)
implicit none
! number of nodes in an 1d element
integer,intent(in) :: nenod
integer :: i
! point where to calculate lagrange function and its derivative
real(kind=kreal),intent(in) :: nd,xi
real(kind=kreal),intent(in) :: ddir,gam
! will be computed only for nenod-1
real(kind=kreal),dimension(nenod-1),intent(out) :: phi,dphi_dxi
real(kind=kreal),dimension(nenod) :: xii!,term,dterm,sum_term
real(kind=kreal) :: A0,A1,B0,B1,dx,fac,rxi1,rxi2

if(nenod.ne.3)then
  write(*,*)'ERROR: infinite element is currently implemented only for 3 nodes!'
  stop
endif

phi=0.0_kreal; dphi_dxi=0.0_kreal

! compute natural coordnates
dx=2.0_kreal/real((nenod-1),kreal)! length = 2.0 as xi is taken -1 to +1
do i=1,nenod
  ! coordinates when origin is in the left
  xii(i)=real((i-1),kreal)*dx
enddo

! origin is transformed to mid point
xii=xii-1.0_kreal
if(ddir>0)then
  rxi1=decay_function(nd,ddir,gam,xii(1))
  rxi2=decay_function(nd,ddir,gam,xii(2))
else
  rxi1=decay_function(nd,ddir,gam,xii(2))
  rxi2=decay_function(nd,ddir,gam,xii(3))
endif

fac=1.0_kreal/(rxi1-rxi2)

A0=-rxi2*fac
A1=fac
B0=rxi1*fac
B1=-fac

phi(1)=A0+A1*decay_function(nd,ddir,gam,xi)
phi(2)=B0+B1*decay_function(nd,ddir,gam,xi)

dphi_dxi(1)=A1*ddecay_function(nd,ddir,gam,xi)
dphi_dxi(2)=B1*ddecay_function(nd,ddir,gam,xi)

return
end subroutine lagrange_infinite
!===============================================================================

! decay function can be defined independent of the pole distance (a) for
! 1/r^n decay type
! this function computes the decay function r(xi) for 1/r decay type
function decay_function(nd,ddir,gam,xi) result(r)
implicit none
!integer,intent(in) :: n
real(kind=kreal),intent(in) :: ddir,gam,nd,xi
real(kind=kreal) :: recn,r
real(kind=kreal),parameter :: two=2.0_kreal
!print*,gam; stop
recn=1.0_kreal/nd
r=(two**recn)*gam/(xi*xi*(gam**nd-two)-ddir*xi*gam**nd+two)**recn
end function decay_function
!=i=============================================================================

! this function computes the decay function r(xi) for 1/r decay type
function ddecay_function(nd,ddir,gam,xi) result(dr)
implicit none
!integer,intent(in) :: n
real(kind=kreal),intent(in) :: ddir,gam,nd,xi
real(kind=kreal) :: dr
real(kind=kreal) :: fac,gamn,recn
real(kind=kreal),parameter :: two=2.0_kreal
gamn=gam**nd
recn=1.0_kreal/nd
fac=1.0_kreal/(xi*xi*(gamn-two)-ddir*xi*gamn+two)**(1.0_kreal+recn)
dr=-recn*(two**recn)*gam*fac*(two*xi*(gamn-two)-ddir*gamn)
end function ddecay_function
!===============================================================================

! this subroutine extract and reorder the nodes of the finite region
! of an element such that the decaying direction always match with 
! the mapping element
subroutine get_gnodinf1(iface,ngllx,nglly,ngllz,ngllinf,gnodinf,poleid)
implicit none
integer,intent(in) :: iface,ngllx,nglly,ngllz,ngllinf
integer,intent(out) :: gnodinf(ngllinf),poleid(9)
integer :: i,j,k,inum,ip

if(iface.lt.1.or.iface.gt.6)then
  write(*,*)'ERROR: illegal outer face ID:',iface
  stop
endif
inum=0; ip=0
if(iface==1)then
  do k=1,ngllz
    do i=1,ngllx
      do j=nglly,2,-1
        inum=inum+1
        gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
        if(j==nglly)then
          ip=ip+1
          poleid(ip)=inum
        endif
      enddo
    enddo
  enddo
elseif(iface==2)then
  do k=1,ngllz
    do j=1,nglly
      do i=1,ngllx-1
        inum=inum+1
        gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
        if(i==1)then
          ip=ip+1
          poleid(ip)=inum
        endif
      enddo
    enddo
  enddo
elseif(iface==3)then
  do i=1,ngllx
    do k=1,ngllz
      do j=1,nglly-1
        inum=inum+1
        gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
        if(j==1)then
          ip=ip+1
          poleid(ip)=inum
        endif
      enddo
    enddo
  enddo
elseif(iface==4)then
  do j=1,nglly
    do k=1,ngllz
      do i=ngllx,2,-1
        inum=inum+1
        gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
        if(i==ngllx)then
          ip=ip+1
          poleid(ip)=inum
        endif
      enddo
    enddo
  enddo
elseif(iface==5)then
  do i=1,ngllx
    do j=1,nglly
      do k=ngllz,2,-1
        inum=inum+1
        gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
        if(k==ngllz)then
          ip=ip+1
          poleid(ip)=inum
        endif
      enddo
    enddo
  enddo
elseif(iface==6)then
   do j=1,nglly
    do i=1,ngllx
      do k=1,ngllz-1
        inum=inum+1
        gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
        if(k==1)then
          ip=ip+1
          poleid(ip)=inum
        endif
      enddo
    enddo
  enddo
endif
if(inum.ne.ngllinf)then
  !print*,iface,inum,ngllinf
  !stop
endif
return
end subroutine get_gnodinf1
!===============================================================================

! Revision:
!   HNG, Apr 19,2012
! The subroutine radau_quadrature computes a Radau quadrature rule.
! the Radau rule is distinguished by the fact that the left endpoint
! (-1) is always an abscissa.
!
! the integral:
! integral ( -1 <= x <= 1 ) f(x) dx
!
! the quadrature rule:
! sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
! the quadrature rule will integrate exactly all polynomials up to
! X**(2*N-2).
!
! Licensing:!
! this code is distributed under the GNU LGPL license.
!
! Modified:
! 06 February 2007
!
! Author:
!    Original MATLAB code by Greg von Winckel.
!    This MATLAB version by John Burkardt.
!
! References:
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
!    Spectral Methods in Fluid Dynamics,
!    Springer, 1993,
!    ISNB13: 978-3540522058,
!    LC: QA377.S676.
!
!    Francis Hildebrand, 
!    Section 8.11,
!    Introduction to Numerical Analysis,
!    Dover, 1987,
!    ISBN13: 978-0486653631,
!    LC: QA300.H5.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
! Input:
! N: the order or the number of integration points (>0, integer)
! Output:
! X(N): the abscissas
! W(N): the weights
subroutine radau_quadrature(n,x,w)
use math_constants
implicit none
integer,intent(in) :: n
integer :: i,j
real(kind=kreal),intent(out) :: x(n),w(n)
real(kind=kreal) :: rj,rn,fac,p(n,n+1),xold(n)

if(n.lt.1)then
  write(*,*)'ERROR: number of quadrature points must be > 1!'
  stop
endif

x=zero; w=zero
rn=real(n,kreal)

! initial estimate for the abscissas is the Chebyshev-Gauss-Radau nodes.
fac=two*pi/(two*rn-one)

! initialize the Legendre Vandermonde matrix.
p=zero
p(2:n,1) = one;
do i=1,n
  x(i)=-cos(fac*real(i-1,kreal))
  p(1,i) = (-one)**(i-1)
enddo
p(1,n+1)=(-one)**(n)

! compute P using the recursion relation.
! compute its first and second derivatives and 
! update X using the Newton-Raphson method.
xold=two
do i=1,100
  if(maxval(abs(x-xold)).le.zerotol)exit
  if(i.ge.100)then
    write(*,*)'ERROR: Legendre Vandermonde matrix does not converge!'
    stop
  endif  
  xold = x; 
  p(2:n,2) = x(2:n);       
  do j=2,n
    rj=real(j,kreal)
    p(2:n,j+1) = ((two*rj-one)*x(2:n)*p(2:n,j)+(-rj+one)*p(2:n,j-1))/rj
  enddo     
  x(2:n) = xold(2:n)- &
           ((one-xold(2:n))/rn)*(p(2:n,n)+p(2:n,n+1))/(p(2:n,n)-p(2:n,n+1))
enddo

! compute the weights.
w = zero
w(1) = two/(rn*rn)
w(2:n)=(one-x(2:n))/(rn*p(2:n,n)*rn*p(2:n,n))
return
end subroutine radau_quadrature
!===============================================================================

end module infinite_element
!===============================================================================
