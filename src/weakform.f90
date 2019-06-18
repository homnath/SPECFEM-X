! TODO: lot of optimization is possible using orthogonality
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
! compute weighting terms of the weak form
module weakform
use set_precision
use global,only:ndim
use math_constants,only:ZERO
implicit none
character(len=250),private :: myfname=" => weakform.f90"
character(len=500),private :: errsrc
contains
!-------------------------------------------------------------------------------

! This subroutine forms the strain-displacement matrix (bmat)
! REFERENCE:
!  modified from
!  Smith and Griffiths (2004): Programming the finite element method
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
! Engineering shear strain is used for concise representation. Hence, no 1/2
! factor below.
! Engineering shear strain = 2 * Ordinary shear strain
! NOTE: strain is symmetric and its components are in the order:
! 1: exx
! 2: eyy
! 3: ezz
! 4: exy
! 5: eyz
! 6: ezx
subroutine compute_bmat_stress(deriv,bmat_stress)
implicit none
real(kind=kreal),intent(in) :: deriv(:,:)
real(kind=kreal),intent(out) :: bmat_stress(:,:)
integer :: n,nod,nst,j1,j2,j3
real(kind=kreal) :: dx,dy,dz

errsrc=trim(myfname)//' => compute_bmat_stress'

! check size
nst=ubound(bmat_stress,1)
if(nst.ne.6)then
  write(*,*)'ERROR: wrong size of the stress tensor!'//trim(errsrc)
  stop
endif
nod=ubound(deriv,2)

!bmat=ZERO 
!DO m=1,nod
!  n=3*m !j3
!  k=n-1 !j2
!  l=k-1 !j1
!  dx=deriv(1,m)
!  dy=deriv(2,m)
!  dz=deriv(3,m)
!  bmat(1,l)=dx
!  bmat(4,k)=dx
!  bmat(6,n)=dx
!  bmat(2,k)=dy
!  bmat(4,l)=dy
!  bmat(5,n)=dy
!  bmat(3,n)=dz
!  bmat(5,k)=dz
!  bmat(6,l)=dz
!end DO

bmat_stress=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  dx=deriv(1,n)
  dy=deriv(2,n)
  dz=deriv(3,n)

  bmat_stress(1,j1)=dx
  bmat_stress(2,j2)=dy
  bmat_stress(3,j3)=dz
  bmat_stress(4,j1)=dy
  bmat_stress(4,j2)=dx
  bmat_stress(5,j2)=dz 
  bmat_stress(5,j3)=dy
  bmat_stress(6,j1)=dz
  bmat_stress(6,j3)=dx
enddo
return
end subroutine compute_bmat_stress
!===============================================================================

subroutine compute_rmat_term1(rho,interpf,rmat_term1)
implicit none
real(kind=kreal),intent(in) :: rho
real(kind=kreal),intent(in) :: interpf(:)
real(kind=kreal),intent(out) :: rmat_term1(:,:)
integer :: n,nod,nst,j1,j2,j3

errsrc=trim(myfname)//' => compute_rmat_term1'

! check size
nst=ubound(rmat_term1,1)
if(nst.ne.3)then
  write(*,*)'ERROR: wrong size of the term1!'//trim(errsrc)
  print*,'hello!',nst
  stop
endif
nod=ubound(interpf,1)

rmat_term1=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  rmat_term1(1,j1)=interpf(n)*rho
  rmat_term1(2,j2)=interpf(n)*rho
  rmat_term1(3,j3)=interpf(n)*rho
enddo
return
end subroutine compute_rmat_term1
!===============================================================================

subroutine compute_rmat_term2(g,dg,interpf,dinterpf,rmat_term2)
implicit none
!g: (gx,gy,gz)
!dg: dg(1)=dgmat(1,1)
!    dg(2)=dgmat(2,2)
!    dg(3)=dgmat(3,3)
!    dg(4)=dgmat(1,2)
!    dg(5)=dgmat(1,3)
!    dg(6)=dgmat(2,3)
!dgmat is symmetric by way of second derivative of the gravitational potential
real(kind=kreal),intent(in) :: g(:),dg(:),interpf(:),dinterpf(:,:)
real(kind=kreal),intent(out) :: rmat_term2(:,:)
integer :: n,nod,nst,j1,j2,j3

errsrc=trim(myfname)//' => compute_rmat_term2'

! check size
nst=ubound(rmat_term2,1)
if(nst.ne.3)then
  write(*,*)'ERROR: wrong size of the term2!'//trim(errsrc)
  stop
endif
nod=ubound(dinterpf,2)

rmat_term2=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  rmat_term2(1,j1)=interpf(n)*dg(1)+dinterpf(1,n)*g(1)
  rmat_term2(2,j1)=interpf(n)*dg(4)+dinterpf(1,n)*g(1)
  rmat_term2(3,j1)=interpf(n)*dg(5)+dinterpf(1,n)*g(1)
               
  rmat_term2(1,j2)=interpf(n)*dg(4)+dinterpf(2,n)*g(2)
  rmat_term2(2,j2)=interpf(n)*dg(2)+dinterpf(2,n)*g(2)
  rmat_term2(3,j2)=interpf(n)*dg(6)+dinterpf(2,n)*g(2)
               
  rmat_term2(1,j3)=interpf(n)*dg(5)+dinterpf(3,n)*g(3)
  rmat_term2(2,j3)=interpf(n)*dg(6)+dinterpf(3,n)*g(3)
  rmat_term2(3,j3)=interpf(n)*dg(3)+dinterpf(3,n)*g(3)
enddo
return
end subroutine compute_rmat_term2
!===============================================================================

subroutine compute_rmat_term3(rho,grav,interpf,rmat_term3)
implicit none
real(kind=kreal),intent(in) :: rho,grav(:)
real(kind=kreal),intent(in) :: interpf(:)
real(kind=kreal),intent(out) :: rmat_term3(:,:)
integer :: n,nod,nst,j1,j2,j3

errsrc=trim(myfname)//' => compute_rmat_term3'

! check size
nst=ubound(rmat_term3,1)
if(nst.ne.1)then
  write(*,*)'ERROR: wrong size of the term4!'//trim(errsrc)
  print*,'hello!',nst
  stop
endif
nod=ubound(interpf,1)

rmat_term3=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  !rmat_term3(1,j1)=interpf(n)*rho(n)*grav(1,n)
  !rmat_term3(1,j2)=interpf(n)*rho(n)*grav(2,n)
  !rmat_term3(1,j3)=interpf(n)*rho(n)*grav(3,n)
  rmat_term3(1,j1)=interpf(n)*rho*grav(1)
  rmat_term3(1,j2)=interpf(n)*rho*grav(2)
  rmat_term3(1,j3)=interpf(n)*rho*grav(3)
enddo
return
end subroutine compute_rmat_term3
!===============================================================================

subroutine compute_rmat_term4(rho,dinterpf,rmat_term4)
implicit none
real(kind=kreal),intent(in) :: rho
real(kind=kreal),intent(in) :: dinterpf(:,:)
real(kind=kreal),intent(out) :: rmat_term4(:,:)
integer :: n,nod,nst,j1,j2,j3

errsrc=trim(myfname)//' => compute_rmat_term4'

! check size
nst=ubound(rmat_term4,1)
if(nst.ne.1)then
  write(*,*)'ERROR: wrong size of the term4!'//trim(errsrc)
  print*,'hello!',nst
  stop
endif
nod=ubound(dinterpf,2)

rmat_term4=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  !rmat_term4(1,j1)=dinterpf(1,n)*rho(n)
  !rmat_term4(1,j2)=dinterpf(2,n)*rho(n)
  !rmat_term4(1,j3)=dinterpf(3,n)*rho(n)
  rmat_term4(1,j1)=dinterpf(1,n)*rho
  rmat_term4(1,j2)=dinterpf(2,n)*rho
  rmat_term4(1,j3)=dinterpf(3,n)*rho
enddo
return
end subroutine compute_rmat_term4
!===============================================================================

!subroutine compute_wmat_term1(grav,hmat,interpf,dinterpf,wmat_term1)
!implicit none
!real(kind=kreal),intent(in) :: grav(:),hmat(:),interpf(:),dinterpf(:,:)
!real(kind=kreal),intent(out) :: wmat_term1(:,:)
!integer :: n,nod,nst,j1,j2,j3
!
!errsrc=trim(myfname)//' => compute_wmat_term1'
!
!! check size
!nst=ubound(wmat_term1,2)
!if(nst.ne.3)then
!  write(*,*)'ERROR: wrong size of the term1!'//trim(errsrc)
!  stop
!endif
!nod=ubound(dinterpf,2)
!
!wmat_term1=ZERO
!j3=0
!do n=1,nod
!  j1=j3+1; j2=j1+1; j3=j2+1
!
!  wmat_term1(j1,1)=interpf(n)*hmat(1)+dinterpf(1,n)*grav(1)
!  wmat_term1(j1,2)=interpf(n)*hmat(4)+dinterpf(1,n)*grav(2)
!  wmat_term1(j1,3)=interpf(n)*hmat(5)+dinterpf(1,n)*grav(3)
!
!  wmat_term1(j2,1)=interpf(n)*hmat(4)+dinterpf(2,n)*grav(1)
!  wmat_term1(j2,2)=interpf(n)*hmat(2)+dinterpf(2,n)*grav(2)
!  wmat_term1(j2,3)=interpf(n)*hmat(6)+dinterpf(2,n)*grav(3)
!  
!  wmat_term1(j3,1)=interpf(n)*hmat(5)+dinterpf(3,n)*grav(1)
!  wmat_term1(j3,2)=interpf(n)*hmat(6)+dinterpf(3,n)*grav(2)
!  wmat_term1(j3,3)=interpf(n)*hmat(3)+dinterpf(3,n)*grav(3)
!enddo
!return
!end subroutine compute_wmat_term1
!!===============================================================================

subroutine compute_wmat_term1(g,dg,interpf,dinterpf,wmat_term1)
implicit none
!g: (gx,gy,gz)
!dg: dg(1)=dgmat(1,1)
!    dg(2)=dgmat(2,2)
!    dg(3)=dgmat(3,3)
!    dg(4)=dgmat(1,2)
!    dg(5)=dgmat(1,3)
!    dg(6)=dgmat(2,3)
!dgmat is symmetric by way of second derivative of the gravitational potential
real(kind=kreal),intent(in) :: g(:),dg(:),interpf(:),dinterpf(:,:)
real(kind=kreal),intent(out) :: wmat_term1(:,:)

integer :: n,nod,ncol,j1,j2,j3

errsrc=trim(myfname)//' => compute_wmat_term1'

! check size
ncol=ubound(wmat_term1,2)
if(ncol.ne.ndim)then
  write(*,*)'ERROR: wrong size of the term1!'//trim(errsrc)
  stop
endif
nod=ubound(dinterpf,2)

wmat_term1=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  wmat_term1(j1,1)=interpf(n)*dg(1)+dinterpf(1,n)*g(1)
  wmat_term1(j1,2)=interpf(n)*dg(4)+dinterpf(1,n)*g(1)
  wmat_term1(j1,3)=interpf(n)*dg(5)+dinterpf(1,n)*g(1)

  wmat_term1(j2,1)=interpf(n)*dg(4)+dinterpf(2,n)*g(2)
  wmat_term1(j2,2)=interpf(n)*dg(2)+dinterpf(2,n)*g(2)
  wmat_term1(j2,3)=interpf(n)*dg(6)+dinterpf(2,n)*g(2)
  
  wmat_term1(j3,1)=interpf(n)*dg(5)+dinterpf(3,n)*g(3)
  wmat_term1(j3,2)=interpf(n)*dg(6)+dinterpf(3,n)*g(3)
  wmat_term1(j3,3)=interpf(n)*dg(3)+dinterpf(3,n)*g(3)
enddo
return
end subroutine compute_wmat_term1
!===============================================================================

subroutine compute_wmat_term2(rho,interpf,wmat_term2)
implicit none
real(kind=kreal),intent(in) :: rho
real(kind=kreal),intent(in) :: interpf(:)
real(kind=kreal),intent(out) :: wmat_term2(:,:)
integer :: n,nod,nst,j1,j2,j3

errsrc=trim(myfname)//' => compute_wmat_term2'

! check size
nst=ubound(wmat_term2,2)
if(nst.ne.3)then
  write(*,*)'ERROR: wrong size of the term2!'//trim(errsrc)
  print*,'hello!',nst
  stop
endif
nod=ubound(interpf,1)

wmat_term2=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  wmat_term2(j1,1)=interpf(n)*rho
  wmat_term2(j2,2)=interpf(n)*rho
  wmat_term2(j3,3)=interpf(n)*rho
enddo
return
end subroutine compute_wmat_term2
!===============================================================================

subroutine compute_wmat_term3(dinterpf,wmat_term3)
implicit none
real(kind=kreal),intent(in) :: dinterpf(:,:)
real(kind=kreal),intent(out) :: wmat_term3(:,:)
integer :: n,nod,nst,j1,j2,j3

errsrc=trim(myfname)//' => compute_wmat_term3'

! check size
nst=ubound(wmat_term3,2)
if(nst.ne.1)then
  write(*,*)'ERROR: wrong size of the term3!'//trim(errsrc)
  print*,'hello!',nst
  stop
endif
nod=ubound(dinterpf,2)

wmat_term3=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  wmat_term3(j1,1)=dinterpf(1,n)
  wmat_term3(j2,1)=dinterpf(2,n)
  wmat_term3(j3,1)=dinterpf(3,n)
enddo
return
end subroutine compute_wmat_term3
!===============================================================================

subroutine compute_wmat_term4(grav,interpf,wmat_term4)
implicit none
real(kind=kreal),intent(in) :: grav(:)
real(kind=kreal),intent(in) :: interpf(:)
real(kind=kreal),intent(out) :: wmat_term4(:,:)
integer :: n,nod,nst,j1,j2,j3

errsrc=trim(myfname)//' => compute_wmat_term4'

! check size
nst=ubound(wmat_term4,2)
if(nst.ne.1)then
  write(*,*)'ERROR: wrong size of the term4!'//trim(errsrc)
  print*,'hello!',nst
  stop
endif
nod=ubound(interpf,1)

wmat_term4=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1

  wmat_term4(j1,1)=interpf(n)*grav(1)
  wmat_term4(j2,1)=interpf(n)*grav(2)
  wmat_term4(j3,1)=interpf(n)*grav(3)
enddo
return
end subroutine compute_wmat_term4
!===============================================================================

subroutine compute_wmat_lf(interpf,wmat_lf)
implicit none
real(kind=kreal),intent(in) :: interpf(:)
real(kind=kreal),intent(out) :: wmat_lf(:)
integer :: n,nod,j1,j2,j3

errsrc=trim(myfname)//' => compute_wmat_lf'

nod=ubound(interpf,1)

wmat_lf=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1
  wmat_lf(j1)=interpf(n)
  wmat_lf(j2)=interpf(n)
  wmat_lf(j3)=interpf(n)
enddo
return
end subroutine compute_wmat_lf
!===============================================================================

subroutine compute_wmat_gradphi(interpf,wmat_gradphi)
implicit none
real(kind=kreal),intent(in) :: interpf(:)
real(kind=kreal),intent(out) :: wmat_gradphi(:,:)
integer :: n,nod,nst,j1,j2,j3

errsrc=trim(myfname)//' => compute_wmat_gradphi'

! check size
nst=ubound(wmat_gradphi,2)
if(nst.ne.3)then
  write(*,*)'ERROR: wrong size of the gradphi!'//trim(errsrc)
  stop
endif
nod=ubound(interpf,1)

wmat_gradphi=ZERO
j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1
  wmat_gradphi(j1,1)=interpf(n)
  wmat_gradphi(j2,2)=interpf(n)
  wmat_gradphi(j3,3)=interpf(n)
enddo
return
end subroutine compute_wmat_gradphi
!===============================================================================

subroutine compute_rmat_gradphi(rho,dinterpf,rmat_gradphi)
use global,only:NDIM,NGLL
implicit none
real(kind=kreal),intent(in) :: rho
real(kind=kreal),intent(in) :: dinterpf(ndim,NGLL)
real(kind=kreal),intent(out) :: rmat_gradphi(:,:)
integer :: n,nod,nst
real(kind=kreal) :: rhodN(ndim,NGLL)

errsrc=trim(myfname)//' => compute_rmat_gradphi'

! check size
nst=ubound(rmat_gradphi,1)
if(nst.ne.3)then
  write(*,*)'ERROR: wrong size of the (gradphi)!'//trim(errsrc)
  stop
endif
nod=ubound(dinterpf,2)

rhodN=rho*dinterpf

rmat_gradphi=ZERO

do n=1,nod
  rmat_gradphi(1,n)=rhodN(1,n)
  rmat_gradphi(2,n)=rhodN(2,n)
  rmat_gradphi(3,n)=rhodN(3,n)
enddo
return
end subroutine compute_rmat_gradphi
!===============================================================================

subroutine compute_rmat_sPE(rho,interpf,rmat_sPE)
use global,only:NDIM,NGLL
implicit none
real(kind=kreal),intent(in) :: rho,interpf(:)
real(kind=kreal),intent(out) :: rmat_sPE(:,:)
real(kind=kreal) :: facN(NGLL)
integer :: n,nod,nst,j1,j2,j3
!real(kind=kreal),parameter :: FOUR_PI_G=4.0_kreal*PI*GRAV

errsrc=trim(myfname)//' => compute_rmat_sPE'

! check size
nst=ubound(rmat_sPE,1)
if(nst.ne.3)then
  write(*,*)'ERROR: wrong size of the sPE!'//trim(errsrc)
  print*,'hello!',nst
  stop
endif
nod=ubound(interpf,1)

! factored interpolation functions
facN=rho*interpf
rmat_sPE=ZERO

j3=0
do n=1,nod
  j1=j3+1; j2=j1+1; j3=j2+1
  rmat_sPE(1,j1)=facN(n)
  rmat_sPE(2,j2)=facN(n)
  rmat_sPE(3,j3)=facN(n)
enddo
return
end subroutine compute_rmat_sPE
!===============================================================================

subroutine compute_wmat_sPE(dinterpf,wmat_sPE)
implicit none
real(kind=kreal),intent(in) :: dinterpf(:,:)
real(kind=kreal),intent(out) :: wmat_sPE(:,:)
integer :: n,nod,nst

errsrc=trim(myfname)//' => compute_wmat_sPE'

! check size
nst=ubound(wmat_sPE,2)
if(nst.ne.3)then
  write(*,*)'ERROR: wrong size of the wmat_sPE!'//trim(errsrc)
  stop
endif
nod=ubound(dinterpf,2)

wmat_sPE=ZERO

do n=1,nod
  wmat_sPE(n,1)=dinterpf(1,n)
  wmat_sPE(n,2)=dinterpf(2,n)
  wmat_sPE(n,3)=dinterpf(3,n)
enddo
return
end subroutine compute_wmat_sPE
!===============================================================================

subroutine compute_rmat_sn(nfgll,uvect,interpf,rmat_sn)
use global,only:nndofu
implicit none
integer,intent(in) :: nfgll
real(kind=kreal),intent(in) :: uvect(:),interpf(:)
real(kind=kreal),intent(out) :: rmat_sn(:,:)
integer :: np,nrow,ncol,ncomp
integer :: ip,j1,j2,j3

errsrc=trim(myfname)//' => compute_rmat_sn'

! check size
np=ubound(interpf,1)
if(np.ne.nfgll)then
  write(*,'(a,1x,i0)')'ERROR: wrong dimension for "interpf" in '// &
  trim(errsrc),np
  stop
endif
ncomp=1
nrow=ubound(rmat_sn,1)
ncol=ubound(rmat_sn,2)
if(nrow.ne.ncomp .or. &
   ncol.ne.nfgll*nndofu)then
  write(*,'(a,1x,i0,1x,i0)')'ERROR: wrong dimensions for "rmat_sn" in '// &
  trim(errsrc),nrow,ncol
  stop
endif

rmat_sn=ZERO
j3=0
do ip=1,np
  j1=j3+1; j2=j1+1; j3=j2+1
  rmat_sn(1,j1)=interpf(ip)*uvect(1)
  rmat_sn(1,j2)=interpf(ip)*uvect(2)
  rmat_sn(1,j3)=interpf(ip)*uvect(3)
enddo
return
end subroutine compute_rmat_sn
!===============================================================================

end module weakform
!===============================================================================
