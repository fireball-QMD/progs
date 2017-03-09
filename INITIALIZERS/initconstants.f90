! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
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


! initconstants.f90
! Program Description
! ===========================================================================
!       This routine initializes some constants_fireball and other parameters at
! the beginning of the main program.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initconstants (sigma, sigmaold, scf_achieved)
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Output
        real, intent (out) :: sigma
        real, intent (out) :: sigmaold
 
        logical, intent (out) :: scf_achieved

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
        sigma = 0.0d0
        sigmaold = 0.0d0
 
! Now set up Kronecker delta and levi-civita symbols.
        delk = 0.0d0
        delk(1,1) = 1.0d0
        delk(2,2) = 1.0d0
        delk(3,3) = 1.0d0
 
        xlevi = 0.0d0
        xlevi(1,2,3) = 1.0d0
        xlevi(1,3,2) = -1.0d0
        xlevi(3,1,2) = 1.0d0
        xlevi(3,2,1) = -1.0d0
        xlevi(2,3,1) = 1.0d0
        xlevi(2,1,3) = -1.0d0

        scf_achieved = .true.
 
! Format Statements
! ===========================================================================
 
        return
        end
