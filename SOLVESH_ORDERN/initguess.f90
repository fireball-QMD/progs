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

 
! initguess.f90
! Program Description
! ===========================================================================
!        Routine to assign an initial guess for an atomic orbital in a 
! localized wave function centered on iatom.  Assigns a random guess if the 
! orbital belongs to the first 'zeta' of the atom (lowest energy shell of its 
! angular momentum), and if the angular momentum is populated in the free 
! atom. Otherwise, sets coefficient to zero.

! ===========================================================================
! This code was originally provided by Pablo Ordejon.
!
! Code rewritten by:
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
        subroutine initguess (natoms, nspecies, iatom, imu, in1, cvalue)
        use charges
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iatom
        integer, intent (in) :: imu
        integer, intent (in) :: in1
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies

! Output
        real, intent (out) :: cvalue

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer index
        integer inu
        integer issh
        integer jssh
        integer Lmax, Lmaxp
        integer Lvalue_1, Lvalue_2

        real xrandom

        logical first
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize cvalue to zero
        cvalue = 0.0d0

! Find out angular momentum of orbital iorb
        index = 0
        do issh = 1, nssh(in1)
         Lvalue_1 = lssh(issh,in1)
         first = .true.
         if (issh .gt. 1) then
          do jssh = 1, issh - 1
           Lvalue_2 = lssh(jssh,in1) 
           if (Lvalue_1 .eq. Lvalue_2) first = .false.
          end do
         end if
         Lmax = 2*Lvalue_1 + 1 
         do inu = 1, Lmax
          index = index + 1
         end do
         if (index .ge. imu) exit
        end do

! Return if orbital is not the first zeta of its shell.
        if (.not. first) return

! Assign initial guess.
! If 2 or less electrons, populate lowest s orbital
! If 8 or less electrons, populate lowest s and p  orbitals
! If 18 or less electrons, populate lowest s, p and d orbitals
! If 32 or less electrons, populate lowest s, p, d and f orbitals
        if (nelectron(iatom) .le. 32) Lmaxp = 3
        if (nelectron(iatom) .le. 18) Lmaxp = 2
        if (nelectron(iatom) .le. 8) Lmaxp = 1
        if (nelectron(iatom) .le. 2) Lmaxp = 0

        if (Lmaxp .gt. Lmax) then
         write (*,*) ' Cannot build initial guess in initguess.f90 '
         write (*,*) ' Reason: Maximum angular moment for atom ', iatom,     &
                     '         is not large enough. '
         stop
        end if 

        if (nelectron(iatom) .gt. 32) then
         write (*,*) ' Cannot build initial guess in initguess.f90 '
         write (*,*) ' Too many valence electrons in atom ', iatom
         stop
        end if 

        if (Lvalue_1 .le. Lmaxp) then
         index = 0
         do 
          index = index + 1          
          call random_number (xrandom) 
          cvalue = (xrandom - 0.5d0)*2.0d0
          if (abs(cvalue) .gt. 1.0d-5) exit 
         end do
        end if

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
