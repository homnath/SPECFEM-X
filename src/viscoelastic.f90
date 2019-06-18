! DESCRIPTION
!  This module conains the routines for viscoelastic rheology.
! DEVELOPER
!  Hom Nath Gharti, Princeton University
! REVISION
!  HNG, Feb 19, 2016; HNG, Jul 12,2011; HNG, Apr 09,2010; HNG, Dec 08,2010
! TODO
!  - implement power law
module viscoelastic
use math_constants
use global,only:myrank
contains
! ----------------------------------------------------------------------

! TODO modify the below function according to original one
real(kind=kreal) function hvisc(x,expx)
!     Purpose: Compute integration factor for viscoelastic material
!              H_visc(x) = [ 1 - exp(-x)]/x
implicit  none
real(kind=kreal) :: x, expx
if(x.lt.1.d-04) then
  hvisc = 1.d0 - 0.5d0*x*(1.d0 - x/3.d0*(1.d0 &
          - 0.25d0*x*(1.d0 - 0.2d0*x)))
else
  hvisc = (1.d0 - expx)/x
endif
end function hvisc
!===============================================================================

! Compute prony series (1-exp(-x))/x which is unstable/undefined for very 
! small x.
! Purpose: Compute integration factor for viscoelastic material
! visco_dh(x) = [ 1 - exp(-x)]/x
! depend only on the time step & relaxation time
! applicable to both Maxwell and general Maxwell rheology
real(kind=kreal) function visco_dq(x)
implicit  none
real(kind=kreal),intent(in) :: x
integer,parameter :: nterm=10
real(kind=kreal),parameter :: tol=1e-12_kreal
integer :: i
real(kind=kreal) :: factorial,frac,fsign,term

if(x.eq.zero)then
  visco_dq=one
  return
elseif(x.lt.1.e-04_kreal)then
  fsign = one;
  factorial = one;
  frac = one;
  visco_dq = one;

  term=one
  do i=2,nterm
    factorial=factorial*real(i,kreal)
    fsign=-one*fsign
    frac=frac*x
    term=fsign*frac/factorial
    visco_dq=visco_dq+term
    if(term<=tol)return
  enddo ! do
  write(*,*)'ERROR: nonconvergence of prony series!'
  stop
  !visco_dq = 1.d0 - 0.5d0*x*(1.d0 - x/3.d0*(1.d0 &
  !        - 0.25d0*x*(1.d0 - 0.2d0*x)))
else
  visco_dq = (one-exp(-x))/x
  return
  !visco_dq = (1.d0 - expx)/x
endif
end function visco_dq
!===============================================================================

! PURPOSE:
!  Compute stress for Maxwell element
! INPUT:
!  K: Bulk modulus
!  G: Shear modulus
subroutine compute_stress_maxwell(K,G,tratio,eps,e0,q0,sig,vsig)
use global,only:NST
implicit  none
real(kind=kreal),intent(in) :: G,K,tratio
real(kind=kreal),intent(in) :: eps(NST)
real(kind=kreal),intent(inout) :: e0(NST),q0(NST)
real(kind=kreal),intent(out) :: sig(NST),vsig(NST)

integer :: i
real(kind=kreal) :: Kth,expx,dq_n,theta
real(kind=kreal) :: e(NST)

! Compute volumetric strain and deviatoric components
theta = ONE_THIRD*(eps(1) + eps(2) + eps(3)) ! multiplied by a factor only to
!reduce computation steps later

do i = 1,3
  e(i) = eps(i) - theta
enddo ! i

! Original shear strain was a engineering shear strain, i.e., 2 X ordinary
! strain. Therefore we need to multiply by factor 0.5
do i = 4,NST
  e(i) = eps(i)*HALF
enddo ! i

! Set properties for integrating the q terms
sig=ZERO

expx = exp(-tratio)

dq_n = visco_dq(tratio)
! Update history and compute viscoelastic deviatoric stress
do i = 1,NST
  q0(i) = expx*q0(i) + dq_n*(e(i) - e0(i))
  sig(i)  = sig(i) + q0(i)
enddo ! i

do i = 1,NST
  sig(i) = TWO*G*sig(i)
  e0(i)  = e(i)
enddo ! i
vsig=sig
! Add elastic bulk term
Kth = K*theta*three ! cancel out the factor multiplied above

do i = 1,3
  sig(i) = sig(i) + Kth
enddo ! i
end subroutine compute_stress_maxwell
!===============================================================================

! compute tangent modulus for Maxwell element
subroutine compute_cmat_maxwell(K,G,tratio,cmat)
use global,only:NST
implicit  none
real(kind=kreal),intent(in) :: G,K,tratio
real(kind=kreal),intent(out) :: cmat(NST,NST)

integer :: i,j
real(kind=kreal) :: Gv,KGv,gfac

! vicous factor
gfac = visco_dq(tratio)


cmat=ZERO
Gv = G*gfac
KGv = K - TWO_THIRD*Gv
do j =1,3
  do i = 1,3
    cmat(i,j) = KGv
  enddo
  ! diagonal elements
  cmat(j,j) = cmat(j,j) + TWO*Gv
enddo
! diagonal elements of last three rows
do i = 4,NST
  cmat(i,i) = Gv
enddo ! i
end subroutine compute_cmat_maxwell
!===============================================================================

! TODO: it has to be decomped into two functions similar to above
subroutine compute_stress_genmaxwell(K,G,tratio,muratio,eps,e0,q0,sig,vsig)
!,cmatve,cmat)
use global,only:nmaxwell,NST
implicit  none
real(kind=kreal),intent(in) :: G,K,tratio(nmaxwell)
real(kind=kreal),intent(in) :: muratio(nmaxwell)
real(kind=kreal),intent(in) :: eps(NST)
real(kind=kreal),intent(inout) :: e0(NST),q0(NST,nmaxwell)
real(kind=kreal),intent(out) :: sig(NST),vsig(NST)
integer :: i,j,n
real(kind=kreal) :: Gg,Kg,Kth, gfac,exp_n,mu_0,mu_n,dq_n,theta
!real(kind=kreal) :: cmatve(NST,NST),cmat(NST,NST)
real(kind=kreal) :: e(NST)
save

! Compute volumetric strain and deviatoric components
theta = ONE_THIRD* (eps(1) + eps(2) + eps(3))

do i = 1,3
  e(i) = eps(i) - theta
enddo ! i

do i = 4,NST
  e(i) = HALF*eps(i)
enddo ! i

! Set properties for integrating the q_i terms
sig = ZERO

mu_0 = ZERO
gfac = ZERO

do n = 1,nmaxwell
  mu_n  = muratio(n)
  exp_n = exp(-tratio(n))
  dq_n = mu_n * visco_dq(tratio(n)) !hvisc(tratio(n),exp_n) ! WARNING0: need to check this
  gfac = gfac + dq_n
  mu_0 = mu_0 + mu_n

  ! Update history and compute viscoelastic deviatoric stress
  do i = 1,NST
    q0(i,n) = exp_n*q0(i,n) + dq_n*(e(i) - e0(i)) ! WARNING: need to check this
    sig(i)  = sig(i) + q0(i,n)
  enddo ! i
enddo ! n

! Finish updates and save the strains
mu_0 = ONE - mu_0
gfac = gfac + mu_0

vsig=ZERO
do i = 1,NST
  vsig(i) = TWO*G*sig(i)
  sig(i) = TWO*G*(mu_0*e(i) + sig(i))
  e0(i)  = e(i)
enddo ! i

! Add elastic bulk term
Kth = K*theta*THREE

do i = 1,3
  sig(i) = sig(i) + Kth
enddo ! i

!! Set tangent parameters
!Gg = G*gfac
!Kg = K - TWO_THIRD*Gg
!K  = K - TWO_THIRD*G
!do j =1,3
!  do i = 1,3
!    cmatve(i,j) = Kg
!    cmat(i,j) = K
!  enddo ! i
!  cmatve(j,j) = cmatve(j,j) + TWO*Gg
!  cmat(j,j) = cmat(j,j) + TWO*G
!enddo ! i
!
!do i = 4,NST
!  cmatve(i,i) = Gg
!  cmat(i,i) = G
!enddo ! i
end subroutine compute_stress_genmaxwell
!===============================================================================

! Compute tangent modulus for General Maxwell element
subroutine compute_cmat_genmaxwell(K,G,tratio,muratio,cmat)
use global,only:nmaxwell,NST
implicit  none
real(kind=kreal),intent(in) :: G,K
real(kind=kreal),intent(in) :: tratio(nmaxwell),muratio(nmaxwell)
real(kind=kreal),intent(out) :: cmat(NST,NST)

integer :: i,j,n
real(kind=kreal) :: Gv,KGv,gfac
real(kind=kreal) :: exp_n,muratio_n,tratio_n
real(kind=kreal) :: dq_n,mu_0,mu_n

! vicous factor
mu_0 = ZERO
gfac = ZERO
do n = 1,nmaxwell
  mu_n  = muratio(n)
  tratio_n=tratio(n)
  exp_n = exp(-tratio_n)
  dq_n = mu_n * visco_dq(tratio_n) !hvisc(tratio_n,exp_n) ! WARNING0: need to check this
  gfac = gfac + dq_n
  mu_0 = mu_0 + mu_n
enddo ! n

! Finish updates and save the strains
mu_0 = ONE - mu_0
gfac = gfac + mu_0

cmat=ZERO
Gv = G*gfac
KGv = K - TWO_THIRD*Gv
do j =1,3
  do i = 1,3
    cmat(i,j) = KGv
  enddo
  ! diagonal elements
  cmat(j,j) = cmat(j,j) + TWO*Gv
enddo
! diagonal elements of last three rows
do i = 4,NST
  cmat(i,i) = Gv
enddo ! i
end subroutine compute_cmat_genmaxwell
!===============================================================================

! Compute stress tensor for Standard Linear Solid (SLS) or Zener element
subroutine compute_stress_zener(K,G,tratio,muratio,eps,e0,q0,sig,vsig)
!,cmatve,cmat)
use global,only:NST
implicit  none
real(kind=kreal),intent(in) :: G,K,tratio
real(kind=kreal),intent(in) :: muratio
real(kind=kreal),intent(in) :: eps(NST)
real(kind=kreal),intent(inout) :: e0(NST),q0(NST)
real(kind=kreal),intent(out) :: sig(NST),vsig(NST)
integer :: i,j,n
real(kind=kreal) :: Gg,Kg,Kth, gfac,exp_n,mu_0,mu_n,dq_n,theta
!real(kind=kreal) :: cmatve(NST,NST),cmat(NST,NST)
real(kind=kreal) :: e(NST)
save

! Compute volumetric strain and deviatoric components
theta = ONE_THIRD* (eps(1) + eps(2) + eps(3))

do i = 1,3
  e(i) = eps(i) - theta
enddo ! i

do i = 4,NST
  e(i) = HALF*eps(i)
enddo ! i

! Set properties for integrating the q_i terms
sig = ZERO

mu_0 = ZERO
gfac = ZERO

mu_n  = muratio
exp_n = exp(-tratio)
dq_n = mu_n * visco_dq(tratio) !hvisc(tratio,exp_n) ! WARNING0: need to check this
gfac = gfac + dq_n
mu_0 = mu_0 + mu_n

! Update history and compute viscoelastic deviatoric stress
do i = 1,NST
  q0(i) = exp_n*q0(i) + dq_n*(e(i) - e0(i)) ! WARNING: need to check this
  sig(i)  = sig(i) + q0(i)
enddo ! i

! Finish updates and save the strains
mu_0 = ONE - mu_0
gfac = gfac + mu_0

vsig=ZERO
do i = 1,NST
  vsig(i) = TWO*G*sig(i)
  sig(i) = TWO*G*(mu_0*e(i) + sig(i))
  e0(i)  = e(i)
enddo ! i

! Add elastic bulk term
Kth = K*theta*THREE

do i = 1,3
  sig(i) = sig(i) + Kth
enddo ! i

end subroutine compute_stress_zener
!===============================================================================

! Compute tangent modulus for Standard Linear Solid (SLS) or Zener element
subroutine compute_cmat_zener(K,G,tratio,muratio,cmat)
use global,only:nmaxwell,NST
implicit  none
real(kind=kreal),intent(in) :: G,K
real(kind=kreal),intent(in) :: tratio,muratio
real(kind=kreal),intent(out) :: cmat(NST,NST)

integer :: i,j
real(kind=kreal) :: Gv,KGv,gfac
real(kind=kreal) :: exp_n,muratio_n,tratio_n
real(kind=kreal) :: dq_n,mu_0,mu_n

! vicous factor
mu_0 = ZERO
gfac = ZERO

mu_n  = muratio
tratio_n=tratio
exp_n = exp(-tratio_n)
dq_n = mu_n * visco_dq(tratio_n) !hvisc(tratio_n,exp_n) ! WARNING0: need to check this
gfac = gfac + dq_n
mu_0 = mu_0 + mu_n

! Finish updates and save the strains
mu_0 = ONE - mu_0
gfac = gfac + mu_0

cmat=ZERO
Gv = G*gfac
KGv = K - TWO_THIRD*Gv
do j =1,3
  do i = 1,3
    cmat(i,j) = KGv
  enddo
  ! diagonal elements
  cmat(j,j) = cmat(j,j) + TWO*Gv
enddo
! diagonal elements of last three rows
do i = 4,NST
  cmat(i,i) = Gv
enddo ! i
end subroutine compute_cmat_zener
!===============================================================================

end module viscoelastic
!===============================================================================
