! copyright info:
!
! @Copyright 2006
! Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

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

 
! allocate_neigh.f90
! Program Description
! ===========================================================================
! This routine allocates the neighbor arrays based on the maximum number
! of neighbors. 
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine allocate_neigh (nprocs, my_proc, iordern, &
 &icluster, ivdw, ifixneigh,iwrthampiece, iwrtatom)
	use configuration
	use interactions
	use neighbor_map
	use module_dos
!CHROM
	use options, only: iclassicMD
	use classicMD, only: find_neigh_max_class
!END CHROM
	implicit none
	
! Argument Declaration and Description
! ===========================================================================
! Input
	integer, intent (in) :: icluster
	integer, intent (in) :: ifixneigh
	integer, intent (in) :: iordern
	integer, intent (in) :: ivdw
	integer, intent (in) :: my_proc
	integer, intent (in) :: nprocs
	integer, intent (in) :: iwrthampiece
	integer, intent (in) :: iwrtatom


! Parameters and Data Declaration
! ===========================================================================
 
! Variable Declaration and Description
! ===========================================================================
	integer num_atoms

	character (len=40) fromwhichfile
	logical neighborsfile
 
! Procedure
! ===========================================================================
! If we are fixing the number of neighbors, then get the dimensions for
! the maximum number of neighbors from the NEIGHBORS file.

	if (ifixneigh .eq. 1) then
 		inquire (file = 'NEIGHBORS', exist = neighborsfile)
 		if (.not. neighborsfile) then
 			write (*,*) ' ifixneigh = 1, but there is no NEIGHBORS file! '
 			stop
 		end if 
 		open (unit = 20, file = 'NEIGHBORS', status = 'old')
 		read (20,*) num_atoms, neigh_max, fromwhichfile
 		write (*,101) fromwhichfile
 		if (num_atoms .ne. natoms) then
 			write (*,*) ' The neighbors file that you are using must not belong '
 			write (*,*) ' to the basis file that you are now calculating. '
	 		write (*,*) ' The number of atoms differs between the two. ' 
			write (*,*) ' Check the NEIGHBORS file and start over! '
			stop
	 	end if
 		close (unit = 20)
! For vdw interactions
 		if (ivdw .eq. 1) then
 			inquire (file = 'NEIGHBORS_VDW', exist = neighborsfile)
			if (.not. neighborsfile) then
 				write (*,*) ' ifixneigh = 1, but there is no NEIGHBORS file! '
 				stop
			end if
 		end if
	 	open (unit = 21, file = 'NEIGHBORS', status = 'old')
 		read (21,*) num_atoms, neigh_max_vdw, fromwhichfile
	 	write (*,101) fromwhichfile
 		if (num_atoms .ne. natoms) then
			write (*,*) ' The neighbors file that you are using must not belong '
			write (*,*) ' to the basis file that you are now calculating. '
			write (*,*) ' The number of atoms differs between the two. ' 
			write (*,*) ' Check the NEIGHBORS file and start over! '
			stop
 		end if
	 	close (unit = 21)
! For PP interactions
		inquire (file = 'NEIGHBORS_PP', exist = neighborsfile)
 		if (.not. neighborsfile) then
			write (*,*) ' ifixneigh = 1, but there is no NEIGHBORS_PP file! '
			stop
 		end if 
 		open (unit = 22, file = 'NEIGHBORS_PP', status = 'old')
 		read (22,*) num_atoms, neighPP_max, fromwhichfile
 		write (*,101) fromwhichfile
 		if (num_atoms .ne. natoms) then
			write (*,*) ' The neighbors file that you are using must not belong '
			write (*,*) ' to the basis file that you are now calculating. '
			write (*,*) ' The number of atoms differs between the two. ' 
			write (*,*) ' Check the NEIGHBORS_PP file and start over! '
			stop
	 	end if
 		close (unit = 22)
	else 
		call find_neigh_max (nprocs, my_proc, iordern, icluster, ivdw)
 		call find_neighPP_max (nprocs, my_proc, iordern, icluster)
	end if

	allocate (neigh_b (neigh_max, natoms))
	allocate (neigh_j (neigh_max, natoms))
	allocate (neighn (natoms))
	allocate (neigh_comb (2, neigh_max**2, natoms))
	allocate (neigh_comj (2, neigh_max**2, natoms))
	allocate (neigh_comm (neigh_max**2, natoms))
	allocate (neigh_comn (natoms))
	allocate (neigh_back (natoms, neigh_max))
	allocate (neigh_self (natoms))
 	allocate (nPP_b (neighPP_max, natoms))
	allocate (nPP_j (neighPP_max, natoms))
	allocate (nPP_map (neighPP_max, natoms))
	allocate (nPPn (natoms))
	allocate (nPP_self (natoms))
 	allocate (nPPx_b (neighPP_max, natoms))
	allocate (nPPx_j (neighPP_max, natoms))
	allocate (nPPx_map (neighPP_max, natoms))
	allocate (nPPx_point (neighPP_max, natoms))
	allocate (nPPxn (natoms))
	allocate (nPPx_self (natoms))
! neighPP
 	allocate (neighPP_b (neighPP_max**2, natoms))
 	allocate (neighPP_j (neighPP_max**2, natoms))
 	allocate (neighPPn (natoms))

! 3. party PPcommon pairs
 	allocate (neighPP_comb (2, neighPP_max**2, natoms))
 	allocate (neighPP_comj (2, neighPP_max**2, natoms))
 	allocate (neighPP_comm (neighPP_max**2, natoms))
 	allocate (neighPP_comn (natoms))
 	allocate (neighPP_self (natoms))

! Total neighbor list (mapping together neigh and neighPP part)
 	allocate (neighj_tot (neigh_max+neighPP_max, natoms))
 	allocate (neighb_tot (neigh_max+neighPP_max, natoms))
 	allocate (neighn_tot (natoms))
 	if (ivdw .eq. 1) then
 		allocate (neigh_b_vdw (neigh_max_vdw, natoms))
 		allocate (neigh_j_vdw (neigh_max_vdw, natoms))
 		allocate (neighn_vdw (natoms))
 		allocate (neigh_vdw_self (natoms))
	end if

 	if (iwrtatom .ge. 1) then 
 		allocate (hr_box (numorb_max, numorb_max, natoms,0:(neigh_max+neighPP_max)))
 	end if
!CHROM
	if (iclassicMD>0)then
		neigh_max_class = find_neigh_max_class(icluster)
		allocate(neigh_classic(neigh_max_class,natoms))
		allocate(neighn_classic(natoms))
		allocate(neigh_b_classic(neigh_max_class,natoms))
	endif
!END CHROM

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
	101 format (2x, ' Neighbor file corresponds to basisfile = ', a40)
	return
end subroutine allocate_neigh
