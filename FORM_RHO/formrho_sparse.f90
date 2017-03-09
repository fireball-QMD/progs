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

 
! formrho_sparse.f90
! Program Description
! ===========================================================================
!       This routine creates the density matrices (rho and cape) which are in 
! the neighbor- mapped form (sparse form) from the compact form of the 
! matrices.
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
        subroutine formrho_sparse (natoms, iprows, isendstart)  
        use density
        use interactions
        use neighbor_map
        use ordern
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        integer, intent (in) :: iprows
        integer, intent (in) :: isendstart
    
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu, inu
        integer in1, in2
        integer index
        integer ineigh
        integer iwrtdensity
        integer jatom
        integer mmu, mmup, nnu

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Loop over all atoms
        do iatom = 1, natoms
         in1 = imass(iatom)
	
! Loop over neighbors of iatom
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)

! Loop over orbitals
 	  do imu = 1, num_orb(in1)   
           mmu = imu + degelec(iatom)
	   do inu = 1, num_orb(in2)   
            nnu = inu + degelec(jatom)

! Overlap matrix, Hamiltonian and control vector
            if (mmu .ge. isendstart .and. mmu .le. iprows) then
             mmup = mmu - isendstart + 1
             do index = 1, numrho_local(mmup)
              if (listrho_local(index,mmup) .eq. nnu) then
               rho(imu,inu,ineigh,iatom) = rho_compact_local(index,mmup)
               exit
              end if
             end do
             do index = 1, numcape_local(mmup)
              if (listcape_local(index,mmup) .eq. nnu) then
               cape(imu,inu,ineigh,iatom) = cape_compact_local(index,mmup)
               exit
              end if
             end do
            end if
           end do
          end do
         end do
        end do

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
400     format (9f9.4)
 
        return
        end
