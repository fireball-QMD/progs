! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Motorola, Physical Sciences Research Labs - Alex Demkov
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio University - Dave Drabold

!
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! factorial.f90
! Program Description
! ===========================================================================
!       This computes the factorial and returns it as a real*4
!
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
 
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real function factorial (ifac)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: ifac
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer jj
        integer countit
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        if (ifac .lt. 0) then
         write (*,*) ' Error in factorial  ----  ifac < 0 '
         stop
        else
         countit = 1
! Note: if ifac=0,1 then this do loop is skipped, and factorial=1
         do jj = 2, ifac
          countit = countit*jj
         end do
        end if
        factorial = real(countit)
 
! Format Statements
! ===========================================================================
        return
        end
