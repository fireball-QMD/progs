! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
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

 
! reallocate_neigh.f90
! Program Description
! ===========================================================================
!       This routine reallocates the neighbor arrays if the maximum number
! of neighbors has changed. 
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
        subroutine reallocate_neigh (nprocs, my_proc, iordern,       &
     &                               itheory, itheory_xc, igauss, icluster,  &
     &                               ivdw, iwrthampiece,     &
     &                               iwrtatom, igrid)
        use configuration
        use interactions
        use neighbor_map
        use module_dos
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: igauss
        integer, intent (in) :: iordern
        integer, intent (in) :: itheory
        integer, intent (in) :: itheory_xc
        integer, intent (in) :: ivdw
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs
        integer, intent (in) :: iwrthampiece
        integer, intent (in) :: iwrtatom
        integer, intent (in) :: igrid



! Parameters and Data Declaration
! ===========================================================================
 
! Variable Declaration and Description
! ===========================================================================
        integer neigh_max_old
        integer neigh_max_vdw_old
        integer neighPP_max_old

! Procedure
! ===========================================================================
        neigh_max_old = neigh_max
        neighPP_max_old = neighPP_max
        if(ivdw .eq. 1) neigh_max_vdw_old = neigh_max_vdw

        call find_neigh_max (nprocs, my_proc, iordern, icluster, ivdw)
        if (neigh_max_old .ne. neigh_max) then 
         if (my_proc .eq. 0) then
          write (*,*) ' The maximum number of neighbors has changed. '
          write (*,*) ' Reallocating neighbor_map arrays! '
          write (*,*) ' neigh_max = ', neigh_max
          write (*,*) ' neigh_max_old = ', neigh_max_old
         end if

! deallocation
         deallocate (neigh_b)
         deallocate (neigh_j)
         deallocate (neigh_comb)
         deallocate (neigh_comj)
         deallocate (neigh_comm)
         deallocate (neigh_back)
         deallocate (neigh_pair_a1)
         deallocate (neigh_pair_a2)
         deallocate (neigh_pair_n1)
         deallocate (neigh_pair_n2)


! new allocation
         allocate (neigh_b (neigh_max, natoms))
         allocate (neigh_j (neigh_max, natoms))
         allocate (neigh_comb (2, neigh_max**2, natoms))
         allocate (neigh_comj (2, neigh_max**2, natoms))
         allocate (neigh_comm (neigh_max**2, natoms))
         allocate (neigh_back (natoms, neigh_max))
         allocate (neigh_pair_a1 (neigh_max*natoms))
         allocate (neigh_pair_a2 (neigh_max*natoms))
         allocate (neigh_pair_n1 (neigh_max*natoms))
         allocate (neigh_pair_n2 (neigh_max*natoms))
         

         call reallocate_h (natoms, neigh_max, neighPP_max, itheory,   &
     &                     itheory_xc, igauss)
         call reallocate_f (natoms, neigh_max, neighPP_max, itheory,   &
     &                     itheory_xc, igauss)
! jel-grid
         call reallocate_rho (natoms, neigh_max, neighPP_max,          &
     &                        itheory_xc, igrid)
! end jel-grid
        end if

! VdW part
        if (ivdw .eq. 1) then 
         if ( neigh_max_vdw_old .ne. neigh_max_vdw) then
          deallocate (neigh_b_vdw)
          deallocate (neigh_j_vdw)

          allocate (neigh_b_vdw (neigh_max_vdw, natoms))
          allocate (neigh_j_vdw (neigh_max_vdw, natoms))
         endif 
        end if

! PP part
        call find_neighPP_max (nprocs, my_proc, iordern, icluster)

        if (neighPP_max_old .ne. neighPP_max) then 
         if (my_proc .eq. 0) then
          write (*,*) ' The maximum number of PP-neighbors has changed. '
          write (*,*) ' Reallocating PP-neighbor_map arrays! '
          write (*,*) ' neighPP_max = ', neighPP_max
          write (*,*) ' neighPP_max_old = ', neighPP_max_old
         end if

! reallocation of PP neighbor arrays
         deallocate (nPP_b)
         deallocate (nPP_j)
         deallocate (nPP_map)

         deallocate (nPPx_b)
         deallocate (nPPx_j)
         deallocate (nPPx_map)
         deallocate (nPPx_point)
! neighPP
         deallocate (neighPP_b)
         deallocate (neighPP_j)
! 3. party PP common pairs
         deallocate (neighPP_comb)
         deallocate (neighPP_comj)
         deallocate (neighPP_comm)

! new allocation
         allocate (nPP_b (neighPP_max, natoms))
         allocate (nPP_j (neighPP_max, natoms))
         allocate (nPP_map (neighPP_max, natoms))

         allocate (nPPx_b (neighPP_max, natoms))
         allocate (nPPx_j (neighPP_max, natoms))
         allocate (nPPx_map (neighPP_max, natoms))
         allocate (nPPx_point (neighPP_max, natoms))
! neighPP
         allocate (neighPP_b (neighPP_max**2, natoms))
         allocate (neighPP_j (neighPP_max**2, natoms))
! 3. party PP  common pairs
         allocate (neighPP_comb (2, neighPP_max**2, natoms))
         allocate (neighPP_comj (2, neighPP_max**2, natoms))
         allocate (neighPP_comm (neighPP_max**2, natoms))

! total neighbor list (mapping together neigh and  neighPP part)
         deallocate (neighj_tot)  
         deallocate (neighb_tot)  
         allocate (neighj_tot (neigh_max+neighPP_max, natoms))
         allocate (neighb_tot (neigh_max+neighPP_max, natoms))
         if (iwrtatom .ge. 1) then                   
          deallocate (hr_box)  
          allocate (hr_box (numorb_max, numorb_max, natoms,                  &
     &              0:(neigh_max+neighPP_max)))
         end if

! FIX ME
! In the future we must to optimize this part. This means we separate neighPP 
! and neigh into two different subroutines. To avoid double reallocation in 
! the case as neighPP_max and neigh_max are changed
         call reallocate_h (natoms, neigh_max, neighPP_max, itheory,   &
     &                      itheory_xc, igauss)
         call reallocate_f (natoms, neigh_max, neighPP_max, itheory,   &
     &                      itheory_xc, igauss)
         call reallocate_rho (natoms, neigh_max, neighPP_max, itheory_xc, &
     &                        igrid)

        end if 

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end subroutine reallocate_neigh
