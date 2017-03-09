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

 
! corrector.f90
! Program Description
! ===========================================================================
!       Corrects positions and velocities based on fifth order gear
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
!       subroutine corrector (natoms, xdot, dt, xmass, ftot)
        subroutine corrector (itime, dt, ftot)
        use configuration
        use constants_fireball
        use dimensions
        use fragments
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: itime

        real, intent(in) :: dt

        real, intent(in), dimension (3, natoms) :: ftot

! Output
! xdot(0,:,:) = positions
! xdot(1,:,:) = velocities
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iorder

        real dtfactor
        real, external :: factorial
        real cfactor

        real, dimension (3) :: acceleration
        real, dimension (3, natoms) :: difference
 
! Procedure
! ===========================================================================
! First calculate the accelerations at the predicted points and the 
! differences between the predicted and actual accelerations
        if (itime .eq. 1) then
         do iatom = 1, natoms
          xdot(2,:,iatom) = fovermp*ftot(:,iatom)/xmass(iatom)
         end do 
        end if 
        do iatom = 1, natoms
         acceleration = fovermp*ftot(:,iatom)/xmass(iatom)
         difference(:,iatom) = acceleration - xdot(2,:,iatom)
        end do

! Gear (often fifth-order)
        do iatom = 1, natoms
         do iorder = 0, gear_order
          cfactor = cfac(iorder)

! Use lower order gear for fragments
          if (numfrags .ne. 0) then
           if (fraggots(iatom) .ne. 0) then
            if (iorder .eq. 1 .or. iorder .eq. 2) then
             cfactor = 1.0d0
            else
             cfactor = 0.0d0
            end if
           end if
          end if

          if (iorder .eq. 2) then
           dtfactor = 1.0d0
          else
           dtfactor = dt**(2-iorder)
          end if
          xdot(iorder,:,iatom) = xdot(iorder,:,iatom)                       &
     &     + cfactor*difference(:,iatom)*(factorial(iorder)/2.0d0)*dtfactor
         end do
        end do
 
! Format Statements
! ===========================================================================
 
        return
        end
 
