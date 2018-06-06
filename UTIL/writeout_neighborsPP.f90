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

 
! writeout_neighbors.f90
! Program Description
! ===========================================================================
!       Write out the charges for restart capabilities.
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
        subroutine writeout_neighborsPP (nprocs, iwrtneigh)
        use configuration  
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iwrtneigh
        integer, intent (in) :: nprocs

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer ineigh
        integer jatom
        integer jneigh
        integer jjneigh
        integer mbeta
        integer num_neigh

        real distance

        real, dimension (3) :: dvec
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        if (iwrtneigh .eq. 1) then
! nPP list
         write (*,*) '  '
         write (*,*) ' nPP list : '
         do iatom = 1, natoms
          num_neigh = nPPn(iatom)
          write (*,*) '  ' 
          write (*,*) ' num_neigh = ', num_neigh 
          write (*,*) ' iatom ineigh mbeta jatom   distance            vector '
          write (*,100) 
          do ineigh = 1, num_neigh 
           mbeta = nPP_b(ineigh,iatom) 
           jatom = nPP_j(ineigh,iatom) 
           distance = sqrt((xl(1,mbeta) + ratom(1,jatom) - ratom(1,iatom))**2& 
     &                   + (xl(2,mbeta) + ratom(2,jatom) - ratom(2,iatom))**2& 
     &                   + (xl(3,mbeta) + ratom(3,jatom) - ratom(3,iatom))**2) 
           dvec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom) 
           write (*,101) iatom, ineigh, mbeta, jatom, distance, dvec 
          end do 
         end do 
! nPPx list 
        write (*,*) '  '
         write (*,*) ' nPPx list : '
         do iatom = 1, natoms
          num_neigh = nPPxn(iatom)
          write (*,*) '  ' 
          write (*,*) ' num_neigh = ', num_neigh 
          write (*,*) ' iatom ineigh mbeta jatom   distance            vector '
          write (*,100) 
          do ineigh = 1, num_neigh 
           mbeta = nPPx_b(ineigh,iatom) 
           jatom = nPPx_j(ineigh,iatom) 
           distance = sqrt((xl(1,mbeta) + ratom(1,jatom) - ratom(1,iatom))**2& 
     &                   + (xl(2,mbeta) + ratom(2,jatom) - ratom(2,iatom))**2& 
     &                   + (xl(3,mbeta) + ratom(3,jatom) - ratom(3,iatom))**2) 
           dvec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom) 
           write (*,101) iatom, ineigh, mbeta, jatom, distance, dvec 
          end do 
         end do 
! neighPP list
         write (*,*) '  '
         write (*,*) ' neighPP list : '
         do iatom = 1, natoms
          num_neigh = neighPPn(iatom)
          write (*,*) '  ' 
          write (*,*) ' num_neigh = ', num_neigh 
          write (*,*) ' iatom ineigh mbeta jatom   distance            vector '
          write (*,100) 
          do ineigh = 1, num_neigh 
           mbeta = neighPP_b(ineigh,iatom) 
           jatom = neighPP_j(ineigh,iatom) 
           distance = sqrt((xl(1,mbeta) + ratom(1,jatom) - ratom(1,iatom))**2& 
     &                   + (xl(2,mbeta) + ratom(2,jatom) - ratom(2,iatom))**2& 
     &                   + (xl(3,mbeta) + ratom(3,jatom) - ratom(3,iatom))**2) 
           dvec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom) 
           write (*,101) iatom, ineigh, mbeta, jatom, distance, dvec 
          end do 
         end do 
        end if 

! ****************************************************************************
!
! W R I T E    O U T    N E I G H B O R S    F I L E
! ****************************************************************************
! Open the file NEIGHBORS which contain the neighbor map information for 
! restart purposes.

        !open (unit = 20, file = 'NEIGHBORS_PP', status = 'unknown')
        !write (20,200) natoms, neighPP_max, basisfile
! neighPP list
        do iatom = 1, natoms
         num_neigh = neighPPn(iatom)
         !write (20,*) num_neigh
         do ineigh = 1, num_neigh 
          mbeta = neighPP_b(ineigh,iatom) 
          jatom = neighPP_j(ineigh,iatom) 
          !write (20,*) iatom, ineigh, mbeta, jatom
         end do 
        end do
        write(*,*) '  '
! nPP list
        do iatom = 1, natoms
         num_neigh = nPPn(iatom)
         !write (20,*) num_neigh
         do ineigh = 1, num_neigh 
          mbeta = nPP_b(ineigh,iatom) 
          jatom = nPP_j(ineigh,iatom) 
          jneigh = nPP_map(ineigh,iatom)
          !write (20,*) iatom, ineigh, mbeta, jatom, jneigh
         end do 
        end do
        write(*,*) '  ' 
! nPPx list
        do iatom = 1, natoms
         num_neigh = nPPxn(iatom)
         !write (20,*) num_neigh
         do ineigh = 1, num_neigh 
          mbeta = nPPx_b(ineigh,iatom) 
          jatom = nPPx_j(ineigh,iatom) 
          jneigh = nPPx_map(ineigh,iatom)
          jjneigh = nPPx_point(ineigh,iatom)
          !write (20,*) iatom, ineigh, mbeta, jatom, jneigh, jjneigh
         end do 
        end do
!
! nPPx_self
! nPP_map
! nPPx_map
! nPPx_point
        close (unit = 20)
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (70('='))
101     format (2x, i5, 4x, i2, 3x, i3, 3x, i4, 2x, f9.4, 3x, 3f9.4)
200     format (2x, i5, 2x, i4, 2x, a40)
 
        return
      end subroutine writeout_neighborsPP
