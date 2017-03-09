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


! projectrho.f90
! Program Description
! ===========================================================================
!       Uses projection coefficients to multiply each element of rho, projecting
! out contributions from particular classes of states or atoms.
!
! ===========================================================================
! Code written by:
! Bret Hess
! Department of Physics and Astronomy
! Brigham Young University
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine projectrho (rhomat, norb)
        use configuration
        use interactions
        use neighbor_map
        use dimensions
        use density
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: norb
		complex, intent(inout), dimension (norb,norb) :: rhomat

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer :: iatom
        integer :: imu
        integer :: info
        integer :: inu
        integer :: in1, in2
        integer :: ineigh
        integer :: ishort
        integer :: jatom
        integer :: jmu
        integer :: jnu
        integer :: mbeta
        integer :: mineig

        integer iplace, jplace, Li, Lj, mi, mj, issh, jssh, indexi, indexj
        real ui, uj, rmagi, rmagj



! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Loop through all orbital pairs
 ! We now loop over all neighbors jatom of iatom.
        do iatom = 1, natoms
         in1 = imass(iatom)

! Now loop over all neighbors jatom of iatom.
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)

! So this matrix element goes in the iplace,jplace slot in flat matrix. Need to know type (l,m) of orbital as well.
		 indexi = 0
		 indexj = 0
         do issh = 1, nssh(in1)
           Li = lssh(issh,in1)
           do mi = -Li,Li
             iplace =  degelec(iatom) + indexi
             indexi = indexi + 1
		 	 do jssh = 1, nssh(in2)
              Lj = lssh(jssh,in2)
              do mj = -Lj,Lj
               jplace =  degelec(jatom) + indexj

               !== projection ==
               !pi orbital for c60, radially outward
	      	    rmagi = sqrt((ratom(1,iatom))**2 + (ratom(2,iatom))**2 + (ratom(3,iatom))**2)
	      	    rmagj = sqrt((ratom(1,jatom))**2 + (ratom(2,jatom))**2 + (ratom(3,jatom))**2)
			    if (Li .eq. 1) then
			     if (mi .eq. -1) then !y state
			      ui = (ratom(2,iatom))/rmagi
			     else if (mi .eq. 0) then !z state
			      ui = (ratom(3,iatom))/rmagi
			     else if (mi .eq. 1) then !x state
			      ui = (ratom(1,iatom))/rmagi
			     end if !mi
			    else !s or d state
			     ui = 0.0d0
			    end if !Li

			    if (Lj .eq. 1) then
			     if (mj .eq. -1) then !y state
			      uj = (ratom(2,jatom))/rmagj
			     else if (mj .eq. 0) then !z state
			      uj = (ratom(3,jatom))/rmagj
			     else if (mi .eq. 1) then !x state
			      uj = (ratom(1,jatom))/rmagj
			     end if !mj
			    else !s or d state
			     uj = 0.0d0
			    end if !Lj
  			    rhomat(iplace, jplace) = rhomat(iplace, jplace) * ui * uj
               !== end projection ==

                indexj = indexj + 1
              end do !mj
             end do !jssh
            end do ! mi
           end do ! issh
         end do ! do ineigh
        end do ! do iatom

! Format Statements
! ===========================================================================
 50     format (4i4, 2f12.6, 3f12.6)

        return
      end subroutine projectrho
