!
!                             @Copyright 2010
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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! getforces_classic.f90
! Program Description
! ===========================================================================
! This is used in ALLOCATIONS/allocate_neigh.f90 only one times. It find the
! maximum number of neighbours (size of arrays) in allocation. In the next run
! the arrays are realloced in case that number of neighbors is bigger that array size.
! ===========================================================================
! Code written by:
! Zdenka Chromcova
! Institute of Physics of  the  AS CR,  v. v. i.
! Cukrovarnicka 10
! CZ-162 00 Praha 6
! Czech Republic
! email: chrom@fzu.cz
! webpage with description: http://nanosurf.fzu.cz/wiki/doku.php?id=classical_md
! ===========================================================================

integer function find_neigh_max_class(icluster)
	use configuration, only: natoms, mbeta_max, ratom,xl,nspecies
	use interactions, only: imass,wrtout
	use neighbor_map, only: neigh_max_class
	use classicMD, only: potential
	use options, only: iclassicMD
	implicit none
 
	integer, intent (in) :: icluster
 
	integer iatom
	integer in1, in2
	integer jatom
	integer mbeta
	integer nn
	real distance2

	if (icluster .eq. 1) mbeta_max = 0

	do iatom = 1, natoms
 		nn = 0
 		in1 = imass(iatom)
 		do mbeta = 0, mbeta_max
			do jatom = 1, natoms
 				in2 = imass(jatom)
 				distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2&
 					&+ (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
 					&+ (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
 				if (distance2 < ( potential(in1,in2)%cutoff)**2) then
					nn = nn + 1
 				end if
			end do
 		end do
		if(neigh_max_class < nn)  find_neigh_max_class = nn
	end do	
!write(*,*)'==================', find_neigh_max_class,'================'
end