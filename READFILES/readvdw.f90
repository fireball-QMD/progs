! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
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

 
! initvdw.f90
! Program Description
! ===========================================================================
!       This routine will read the vdw parameters based on the input file
! given by the user - vdw.optional
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
        subroutine readvdw (nspecies, symbolA, ivdw)
        use charges
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input 
        integer, intent (in) :: ivdw
        integer, intent (in) :: nspecies

        character(len=2), intent (in), dimension (nspec_max) :: symbolA 
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer ispec
        integer nzx_temp
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        if (ivdw .eq. 1) then
         open (unit = 27, file = 'vdw.optional', status = 'unknown')
         
         write (*,*) ' Read in the van der Waals parameters.  If you have '
         write (*,*) ' chosen to add in the van der Waals interactions, then '
         write (*,*) ' the parameters will be stored in the file vdw.optional. '
         write (*,*) ' In this file give the cutoff radius for the interaction '
         write (*,*) ' range, and the parameters C6, p_alpha for each atom, '
         write (*,*) ' where the C6''s are the atomic van der Waals terms and '
         write (*,*) ' alpha are the polarizabilities. '

         allocate (C6 (nspec_max))
         allocate (p_alpha (nspec_max))
         allocate (R0 (nspec_max))

         read (27,*) range_vdw
         write (*,*) ' Interaction radius for van der Waals = ', range_vdw 
         write (*,*) '  '
         write (*,*) ' Parameters for van der Waals: '
         write (*,*) '  '
         write (*,100) 
         write (*,101) 
         write (*,100) 
         do ispec = 1, nspecies
          read (27,*) nzx_temp, C6(ispec), p_alpha(ispec), R0(ispec)
          if (nzx_temp .ne. nzx(ispec)) then
           write (*,*) ' This species in the van der Waals parameter file '
           write (*,*) ' does not match the one in the info.dat file. '
           write (*,*) ' nzx from vdW file = ', nzx_temp
           write (*,*) ' nzx from info.dat = ', nzx(ispec)
           write (*,*) ' Fix the vdw.optional file and start again. '
           stop
          end if
          write (*,102) symbolA(ispec), C6(ispec), p_alpha(ispec), R0(ispec)
         end do
         close (unit = 27)
        end if
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, ' Type ', 2x, ' C6 ', 2x, ' alpha ', 2x, ' R0 ')
102     format (5x, a2, 4x, f6.2, 4x, f6.2, 4x, f6.2)

        return
        end
