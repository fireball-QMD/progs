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


! formsh_compact.f90
! Program Description
! ===========================================================================
!       This routine creates a compact form of the Hamiltonian and Overlap
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
        subroutine formsh_compact (natoms, iprows, isendstart)
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
        integer index1, index2
        integer ineigh
        integer iplace
        integer jatom
        integer mmu, nnu

        integer, dimension (:), allocatable :: index
 
! Allocate Arrays
! ===========================================================================
        allocate (index (norbitals))

! Procedure
! ===========================================================================
! Initialize the compact matrix array to zero. 
        numh = 0.0d0
        listh = 0.0d0
 	s_compact = 0.0d0
	h_compact = 0.0d0
        index = 0

! Loop over all atoms
        do iatom = 1, natoms
         in1 = imass(iatom)
	
! Loop over neighbors of iatom
	 index1 = 0
         index2 = 0
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)

          iplace = 1 + degelec(jatom)
          if (index(iplace) .eq. 0) then
           index(iplace) = index1 + 1

! Loop over orbitals of jatom
	   do inu = 1, num_orb(in2)   
	    index1 = index1 + 1

! Loop over orbitals of iatom
 	    do imu = 1, num_orb(in1)   
             mmu = imu + degelec(iatom)

! Overlap matrix, Hamiltonian and control vector
             if (mmu .ge. isendstart .and. mmu .le. iprows) then
              nnu = mmu - isendstart + 1
              s_compact(index1,nnu) = s_mat(imu,inu,ineigh,iatom)
              h_compact(index1,nnu) = h_mat(imu,inu,ineigh,iatom)
              listh(index1,nnu) = inu + degelec(jatom)
             end if
	    end do
           end do
	  else
           index2 = index(iplace) - 1

! Loop over orbitals of atom j
	   do inu = 1, num_orb(in2)   
            index2 = index2 + 1

! Loop over orbitals of atom i
 	    do imu = 1, num_orb(in1)   
             mmu = imu + degelec(iatom)

! Overlap matrix, Hamiltonian and control vector
             if (mmu .ge. isendstart .and. mmu .le. iprows) then
              nnu = mmu - isendstart + 1
              s_compact(index2,nnu) =                                        &
     &         s_compact(index2,nnu) + s_mat(imu,inu,ineigh,iatom)
              h_compact(index2,nnu) =                                        &
     &         h_compact(index2,nnu) + h_mat(imu,inu,ineigh,iatom)
              listh(index2,nnu) = inu + degelec(jatom)
             end if
	    end do
           end do
	  end if
	 end do
         if (index2 .gt. index1) index1 = index2

! Set counter of number of non-zero elements
	 do imu = 1, num_orb(in1)
          mmu = imu + degelec(iatom)
          if (mmu .ge. isendstart .and. mmu .le. iprows)                     & 
     &     numh(mmu - isendstart + 1) = index1
	 end do

! Initialize control vector index
         index = 0
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (index)
 
! Format Statements
! ===========================================================================
 
        return
        end
