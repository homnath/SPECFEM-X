!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  double precision function comp_source_frequency_function(freq,hdur)

  implicit none
  integer,parameter :: SFTYPE=0
  include "constants.h"

  double precision,intent(in) :: freq,hdur
  double precision :: omegath

 if(SFTYPE==0)then
   ! Heaviside function
   comp_source_frequency_function = 1.0d0
 elseif(SFTYPE==1)then
   omegath=TWO*freq*hdur
   comp_source_frequency_function = sin(omegath)/omegath 
 else
   write(*,*)'ERROR: invalid DFTYPE for source frequency function!'
   stop
 endif

  end function comp_source_frequency_function
