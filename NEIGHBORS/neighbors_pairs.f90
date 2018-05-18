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
!
! Code written by:
! ==========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================

subroutine neighbors_pairs (icluster)
    use configuration, only: natoms
    use dimensions
    use interactions
    use neighbor_map
    implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
    integer, intent (in) :: icluster

!$ volatile rcutoff
 
! ===========================================================================
! Local Variable Declaration and Description
! ===========================================================================
      integer :: num_pairs
      integer :: iatom
      integer :: ineigh
      integer :: jatom
      integer :: jneigh

!CAREFUL!!   For the time being, this subroutine only works with icluster = 1.
!It must be extended to work on periodic systems!

!Create list of pairs of neighbors without repetitions: SYMMETRIC FIREBALL,
!APRIL 2018

num_pairs = 0

if (icluster .eq. 1) then
do iatom = 1,natoms

  do ineigh = 1,neighn(iatom)
    jatom = neigh_j(ineigh,iatom)
    if (jatom .gt. iatom) then   !Doubt.. .gt. or .ge.? what's best?

      num_pairs = num_pairs+1
      jneigh = neigh_back(iatom,ineigh)

      neigh_pair_a1(num_pairs) = iatom
      neigh_pair_a2(num_pairs) = jatom
      neigh_pair_n1(num_pairs) = ineigh
      neigh_pair_n2(num_pairs) = jneigh

    end if !if jatom .gt. iatom

  end do !end do ineigh = 1,neighn(iatom) 

end do !end do iatom = 1,natoms

else  !if icluster .eq. 1

!Here the same stuff for periodic systems



end if !end if icluster .eq. 1

tot_pairs = num_pairs    !tot_pairs stores the total number of non-repeated
!pairs of neighbors

!End of creating list of pairs of neighbors without repetitions: SYMMETRIC
!FIREBALL,
!APRIL 2018

return 
end
