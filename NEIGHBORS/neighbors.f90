! copyright info:
!
!                             @Copyright 2002
!                           Fireball Committee
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


! neighbors.f90
! Program Description
! ===========================================================================
!       Find all neighbors to atoms in the central cell.  An atom is at 
! lattice vector xl(mbeta), and basis ratom(iatom).  We refer to this as 
! (mbeta,iatom). An atom in the central cell is at (0,iatom).
! Find all neighbors to atom (0,iatom).
!
! neighn(i)=# of neighbors of atom i
! neighj(i,m)= j-sub-m, the j value of the m'th neighbor.
! neighb(i,m)= beta-sub-m, the beta value for the m'th neighbor.
!
!       The important quantity here is neighj, which indicates which basis
! vector we have: This identifies the species a through the array 
! imass (neighj(iatom,ineigh)) is the type (1 or 2) of the ineigh'th neighbor u
! of atom iatom.  Furthermore, this atom is located at 
! xl(:,neighb(iatom,ineigh)) + ratom(:,neighj(iatom,ineigh)).
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
! ==========================================================================
! This file was changed at 03 Decv 2009 - there were added neigbours according to 
! classic potentials for MD simulations. The main part of the old routine was moved 
! to subroutine fillneig.
! The arrays neigh_XX_classic are alloacted and filled out only in the case that iclassicMD > 0. 
! Zdenka Chromcova, FZU AVCR, Czech rep., kuntova @ fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================

subroutine neighbors (nprocs, my_proc, iordern, icluster,    &
     &                        iwrtneigh, ivdw)
    use configuration, only: nspecies, mbeta_max, natoms
    use options, only: iclassicMD
    use dimensions
    use interactions
    use neighbor_map
    use classicMD,only: potential
    implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
    integer, intent (in) :: icluster
    integer, intent (in) :: iordern
    integer, intent (in) :: ivdw
    integer, intent (in) :: iwrtneigh
    integer, intent (in) :: my_proc
    integer, intent (in) :: nprocs

!$ volatile rcutoff
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
    integer iatomstart
    integer natomsp
!CHROM classic empirical potential(Zdenka Chromcova, chrom@fzu.cz)
	interface
	subroutine fillneigh(iatomstart,natomsp,nprocs, my_proc)
		implicit none
		integer, intent(in) :: nprocs, my_proc,iatomstart,natomsp
	end subroutine fillneigh

	subroutine fillneigh_class(iatomstart,natomsp,nprocs, my_proc)
		integer, intent(in) :: nprocs, my_proc,iatomstart,natomsp
	end subroutine
	end interface
!END CHROM
! Procedure
! ===========================================================================



!    if (my_proc .eq. 0)                                                  &
!     &   write (*,*) ' Welcome to neighbors - determine mapping of neighbors. Npoc=',nprocs,'iordern=',iordern
        
    if (icluster .eq. 1) mbeta_max = 0

! Determine which atoms are assigned to this processor.
    if (iordern .eq. 1) then
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          	natomsp = natomsp + 1
          	iatomstart = natomsp*my_proc + 1
         else
          	iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
     &                + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
    else
         iatomstart = 1
         natomsp = natoms
    end if
    
!CHROM
	if (iclassicMD > 0)then
		call fillneigh_class(iatomstart,natomsp,nprocs, my_proc)
		return
	endif

	call fillneigh(iatomstart,natomsp,nprocs,my_proc)
	
    if (iordern .eq. 1) call neighbors_ordern_final (natoms, nprocs, my_proc, ivdw)
	call writeout_neighbors (nprocs, ivdw, iwrtneigh)
!END CHROM
	return
end subroutine neighbors


!CHROM
subroutine fillneigh(iatomstart,natomsp, nprocs,my_proc)

	use interactions, only: nssh,nsh_max,wrtout,imass
	use neighbor_map, only: neigh_j_vdw,neigh_b_vdw,neigh_b,range_vdw,neighn_vdw,neigh_b_classic,neighn,neigh_j,neigh_max
	use configuration, only: ratom,xl,mbeta_max,rcutoff,natoms
	use options, only: iordern,icluster,ivdw
	use dimensions, only: nspec_max
	implicit none

	integer, intent(in) :: nprocs, my_proc,natomsp,iatomstart

!local variables
	integer :: iatom,jatom,mbeta,num_neigh,num_neigh_vdw,in1,imu,in2,neighcount
	real :: distance2,rcutoff_j, rcutoff_i,distance,range2,rc_max
	integer :: mbeta_max2
	integer :: neigh_max_old, ii, jj

! Loop over all atoms.
!$omp parallel do private (num_neigh, num_neigh_vdw, rcutoff_i, rcutoff_j) &
!$omp&private (in1, in2, distance, distance2, range2, jatom, mbeta)

	rc_max = 0.00
	do iatom = iatomstart, iatomstart - 1 + natomsp
 		num_neigh = 0
		num_neigh_vdw = 0
 		rcutoff_i = 0.0d0
		in1 = imass(iatom)
 		do imu = 1, nssh(in1)
 			if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
		end do

! Loop over all possible neighbors
 		do mbeta = 0, mbeta_max
 			do jatom = 1, natoms
				rcutoff_j = 0.0d0
 				in2 = imass(jatom)
				do imu = 1, nssh(in2)
					if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
 				end do
 
! Find the distance from (mbeta,jatom) to (0,iatom)
				distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2&
					&+ (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
					&+ (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
				distance = sqrt(distance2)
 
! Add a small displacement to the sum of cutoffs.
				range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
 				rc_max = max(rc_max,(rcutoff_i + rcutoff_j))
 
				if (distance2 .le. range2) then
					if (distance2 .lt. 0.7d0 .and. distance .gt. 1.0d-4 .and.&
						 &iatom .ne. jatom .and. wrtout) then
						write (*,*) ' WARNING - atoms dangerously close! '
						write (*,*) ' WARNING - atoms dangerously close! '
						write (*,*) ' WARNING - atoms dangerously close! '
						write (*,*) ' iatom, jatom, distance = ', iatom, jatom, distance
					end if
 					num_neigh = num_neigh + 1
 
! The num_neigh'th neighbor to (0,iatom) at (mbeta,jatom)
 					neigh_j(num_neigh,iatom) = jatom
 					neigh_b(num_neigh,iatom) = mbeta
				end if

! This is a neighbor mapping for including the van der Waals interactions.
 				if (ivdw .eq. 1 .and. distance2 .le. range_vdw**2) then
 					num_neigh_vdw = num_neigh_vdw + 1

! The num_neigh'th neighbor to (0,iatom) at (mbeta,jatom)
 					neigh_j_vdw(num_neigh_vdw,iatom) = jatom
 					neigh_b_vdw(num_neigh_vdw,iatom) = mbeta
				end if
			end do
 		end do
 
! The number of neighbors to atom iatom is num_neigh.
! Remember, that in the total count, the atom itself is a neighbor!
		neighn(iatom) = num_neigh
 		if (ivdw .eq. 1) neighn_vdw(iatom) = num_neigh_vdw
 	end do
 
 	if (iordern .eq. 1)&
		& call neighbors_ordern_final (natoms, nprocs, my_proc, ivdw)

end subroutine fillneigh


subroutine fillneigh_class(iatomstart,natomsp,nprocs, my_proc)

	use neighbor_map, only: neigh_b_classic, neigh_max_class, neighn_classic, neigh_classic
	use configuration, only: ratom, xl, mbeta_max, natoms
	use interactions, only: imass
	use classicMD, only: potential
	use options, only: icluster
	implicit none

	interface
		subroutine reallocate_2D_iarray(x,y,newx,newy,array)
		    integer, intent (in) :: x,y
	        integer, intent(in) :: newx, newy
		    integer, dimension(:,:), allocatable :: array
		end subroutine
	end interface                                                                                                    

	integer, intent(in) :: nprocs, my_proc, iatomstart, natomsp

	integer :: iatom,jatom,mbeta,num_neigh,in1,in2, beta_max
	real :: distance2,range2

	if( icluster == 1 )then  
		beta_max = 0
	else 
		beta_max = mbeta_max	
	endif

	do iatom = iatomstart, iatomstart - 1 + natomsp
 		num_neigh = 0
 		in1 = imass(iatom)
 		do mbeta = 0, beta_max
 			do jatom = 1, natoms
 				in2 = imass(jatom)
				distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2&
					&+ (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
					&+ (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
 
				range2 = (potential(in1,in2)%cutoff)**2
				if (distance2 <= range2) then
					if ( iatom /= jatom  .and. distance2 < 0.7d0 .and. distance2 > 1.0d-8 ) then
						write (*,*) ' WARNING - atoms dangerously close! '
						write (*,*) ' iatom, jatom, distance = ', iatom, jatom, sqrt(distance2)
					end if
 					num_neigh = num_neigh + 1
 					if( num_neigh == neigh_max_class ) then
						neigh_max_class = floor(1.5*num_neigh)
						call reallocate_2D_iarray(num_neigh,natoms,neigh_max_class,natoms,neigh_classic)
						call reallocate_2D_iarray(num_neigh,natoms,neigh_max_class,natoms,neigh_b_classic)
					endif
					if(jatom==0)then
						write(*,*)'ERROR - id of neighbor is 0!'
						stop
					endif
 					neigh_classic(num_neigh,iatom) = jatom
 					neigh_b_classic(num_neigh,iatom) = mbeta
				end if
			end do
 		end do
 		neighn_classic(iatom) = num_neigh
! 		write(*,*)iatom,num_neigh
 	end do
! 	write(*,*)neigh_classic(4,1)
! 	write(*,*)neighn_classic(1)
! 	stop
end subroutine fillneigh_class

subroutine reallocate_2D_iarray(x,y,newx,newy,array)
	implicit none
	integer, intent (in) :: x,y
	integer, intent(in) :: newx, newy
	integer, dimension(:,:), allocatable :: array
	integer, dimension(x,y) :: tmp_array
	integer :: i, j

	do i = 1, x
		do j = 1, y
			tmp_array(i,j)=array(i,j)
		enddo
	enddo

	deallocate(array)
	allocate(array(newx,newy))

	do i = 1, x
		do j = 1, y
			array(i,j)=tmp_array(i,j)
		enddo
	enddo
end subroutine
!END CHROM
