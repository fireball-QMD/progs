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


! initmasses.f90
! Program Description
! ===========================================================================
!       This routine initializes the masses according to the mass of the
! species given in the info.dat file. If the user wishes to make the
! masses of particular atoms any different, then a MASS file can be
! created and masses read from there.
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
        subroutine initmasses (natoms, symbol, smass, xmass)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: natoms
 
        real, dimension(nspec_max), intent(in) :: smass
 
        character (len=2), dimension(natoms), intent(in) :: symbol
 
! Output
        real, intent(out), dimension (natoms) :: xmass
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
        integer num_atoms
 
        logical massfile
 
! Procedure
! ===========================================================================
! By default the input masses are initialized to the masses in info.dat
        xmass = 0.0d0
        do iatom = 1, natoms
         in1 = imass(iatom)
         xmass(iatom) = smass(in1)
        end do
 
! Option
! In some cases the user may wish to make the masses very large or
! different than the elemental masses. For example, in surface calculations.
        inquire (file = 'MASSES', exist = massfile)
        if (massfile) then
         write (*,*) '  '
         write (*,*) ' We are reading from a mass file. '
         open (unit = 12, file = 'MASSES', status = 'old')
         read (12,*) num_atoms
         if (num_atoms .ne. natoms) then
          write (*,*) ' The mass file that you are using must not '
          write (*,*) ' belong to the basis file that you are now '
          write (*,*) ' calculating.  The number of atoms differs '
          write (*,*) ' between the two. '
         end if
         do iatom = 1, natoms
          read (12,*) xmass(iatom)
         end do
         close (unit = 12)
        end if
 
! Write out the input masses.
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' The Input Masses Are: '
        write (*,100)
        write (*,101)
        write (*,100)
        do iatom = 1, natoms
         write (*,102) iatom, symbol(iatom), xmass(iatom)
        end do
        write (*,100)
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, ' Atom # ', 2x, ' Type ', 2x, ' Mass ')
102     format (3x, i5, 7x, a2, 4x, f6.2)
 
        return
        end
