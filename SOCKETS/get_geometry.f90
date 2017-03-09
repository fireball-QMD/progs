! copyright info:
!
!                             @Copyright 2003
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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
! University of Regensburg - Juergen Fritsch

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

 
! get_geometry.f90
! Program Description
! ===========================================================================
! Gets a new geometry from another code that is driving fireball.
! The other code is using fireball to generate the energy and forces.
! Real work is done in C code.
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! LLNL
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine get_geometry(natoms)
        use configuration
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: natoms
 
! Local Variable Declaration and Description
! ===========================================================================
        integer realsize
 
! Procedure
! ===========================================================================
        realsize = 4 ! Real(4)
        if (precision(ratom) .ge. 10) realsize = 8 ! Real(8)
        call soc_recv(ratom,3*natoms*realsize)
        xdot(0,:,1:natoms) = ratom(:,1:natoms)
 
        return
        end
