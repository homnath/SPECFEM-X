module source_function
use set_precision
use math_constants
contains
!-------------------------------------------------------------------------------
real(kind=kreal) function source_frequency_function(freq,hdur)

implicit none
integer,parameter :: SFTYPE=0

real(kind=kreal),intent(in) :: freq,hdur
real(kind=kreal) :: omegath

if(SFTYPE==0)then
  ! Heaviside function
  source_frequency_function = ONE
elseif(SFTYPE==1)then
  omegath=TWO*freq*hdur
  source_frequency_function = sin(omegath)/omegath 
else
  write(*,*)'ERROR: invalid DFTYPE for source frequency function!'
  stop
endif

end function source_frequency_function
!===============================================================================

real(kind=kreal) function source_time_function(t,hdur)

implicit none

real(kind=kreal),intent(in) :: t,hdur


! quasi Heaviside
source_time_function = 0.5d0*(1.0d0 + erf(t/hdur))
! comp_source_time_function = dexp(-(t/hdur)**2)/(dsqrt(PI)*hdur)

end function source_time_function
!===============================================================================

real(kind=kreal) function source_time_function_rickr(t,f0)

implicit none
real(kind=kreal),intent(in) :: t,f0

! ricker
source_time_function_rickr = (ONE-TWO*PI*PI*f0*f0*t*t)*exp(-PI*PI*f0*f0*t*t)

!!! another source time function they have called 'ricker' in some old papers,
!!! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006)
!!! in order to benchmark those simulations, the following formula is needed.
! comp_source_time_function_rickr = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

end function source_time_function_rickr
!===============================================================================
end module source_function
!===============================================================================
