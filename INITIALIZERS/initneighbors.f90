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

 
! initneighbors.f90
! Program Description
! ===========================================================================
!       This routine initializes the neighbors for the simulation. If the
! old neighbor file exists, then this will be used to determine the neighbors
! If not, then the neighbors will be calculated as normal. 
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
        subroutine initneighbors (natoms, ivdw, nstepi)
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ivdw
        integer, intent (in) :: natoms
        integer, intent (in) :: nstepi
 
! Output
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ineigh
        integer jatom
        integer jneigh
        integer jjneigh
        integer katom
        integer kneigh
        integer mbeta
        integer num_neigh
        integer num_neigh_vdw


! Allocate arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! If there is a neighbor file from a previous run, or the user added one,
! then initialize the input neighbor map accordingly.
        if (nstepi .gt. 1) then
            
        open (unit = 20, file = 'NEIGHBORS', status = 'unknown')
        read (20,*) 
        do iatom = 1, natoms
         read (20,*) num_neigh 
         neighn(iatom) = num_neigh
         do ineigh = 1, num_neigh 
          read (20,*) katom, kneigh, mbeta, jatom 
          if (katom .ne. iatom .or. kneigh .ne. ineigh) then
           write (*,*) ' iatom, katom = ', iatom, katom
           write (*,*) ' ineigh, kneigh = ', ineigh, kneigh
           write (*,*) ' Problem in NEIGHBORS, atoms not lined up correctly. '
           write (*,*) ' Fix and restart! '
           stop
          end if 
          neigh_b(ineigh,iatom) = mbeta
          neigh_j(ineigh,iatom) = jatom
         end do 
        end do 
        close (unit = 20)

! Initialize PP neighbor map 
        open (unit = 22, file = 'NEIGHBORS_PP', status = 'unknown')
        read (22,*)
! neighPP list
        do iatom = 1, natoms
         read (22,*) num_neigh 
         neighPPn(iatom) = num_neigh
         do ineigh = 1, num_neigh 
          read (22,*) katom, kneigh, mbeta, jatom 
          if (katom .ne. iatom .or. kneigh .ne. ineigh) then
           write (*,*) ' iatom, katom = ', iatom, katom
           write (*,*) ' ineigh, kneigh = ', ineigh, kneigh
           write (*,*) ' Problem in NEIGHBORS_PP (neighPPn)',           &
     &                 ' atoms not lined up correctly. '
           write (*,*) ' Fix and restart! '
           stop
          end if 
          neighPP_b(ineigh,iatom) = mbeta
          neighPP_j(ineigh,iatom) = jatom
         end do 
        end do 
! nPP list
        do iatom = 1, natoms
         read (22,*) num_neigh 
         nPPn(iatom) = num_neigh
         do ineigh = 1, num_neigh 
          read (22,*) katom, kneigh, mbeta, jatom, jneigh 
          if (katom .ne. iatom .or. kneigh .ne. ineigh) then
           write (*,*) ' iatom, katom = ', iatom, katom
           write (*,*) ' ineigh, kneigh = ', ineigh, kneigh
           write (*,*) ' Problem in NEIGHBORS_PP (nPPn)',               &
     &                 '  atoms not lined up correctly. '
           write (*,*) ' Fix and restart! '
           stop
          end if 
          nPP_b(ineigh,iatom) = mbeta
          nPP_j(ineigh,iatom) = jatom
          nPP_map(ineigh,iatom) = jneigh
         end do 
        end do 
! nPPx list
        do iatom = 1, natoms
         read (22,*) num_neigh 
         nPPxn(iatom) = num_neigh
         do ineigh = 1, num_neigh 
          read (22,*) katom, kneigh, mbeta, jatom, jneigh, jjneigh 
          if (katom .ne. iatom .or. kneigh .ne. ineigh) then
           write (*,*) ' iatom, katom = ', iatom, katom
           write (*,*) ' ineigh, kneigh = ', ineigh, kneigh
           write (*,*) ' Problem in NEIGHBORS, atoms not lined up correctly. '
           write (*,*) ' Fix and restart! '
           stop
          end if 
          nPPx_b(ineigh,iatom) = mbeta
          nPPx_j(ineigh,iatom) = jatom
          nPPx_map(ineigh,iatom) = jneigh
          nPPx_point(ineigh,iatom) = jjneigh
         end do 
        end do 
        close (unit = 22)
        end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   SELF nPP_self  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! What's the neighbor of the atom itself?
        nPP_self = -999
        do iatom = 1, natoms
          do ineigh = 1, nPPn(iatom)
            mbeta = nPP_b(ineigh,iatom)
            jatom = nPP_j(ineigh,iatom)
            if (iatom .eq. jatom .and. mbeta .eq. 0) nPP_self(iatom) = ineigh
          end do
        end do

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   SELF nPPx_self  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        nPPx_self = -999
        do iatom = 1, natoms
          do ineigh = 1, nPPxn(iatom)
            mbeta = nPPx_b(ineigh,iatom)
            jatom = nPPx_j(ineigh,iatom)
            if (iatom .eq. jatom .and. mbeta .eq. 0) nPPx_self(iatom) = ineigh
          end do
        end do


! For van der Waals interactions.        
        if (ivdw .eq. 1) then
         open (unit = 21, file = 'NEIGHBORS_VDW', status = 'unknown')
         read (21,*) 
         do iatom = 1, natoms
          read (21,*) num_neigh_vdw 
          neighn_vdw(iatom) = num_neigh_vdw
          do ineigh = 1, num_neigh_vdw
           read (21,*) katom, kneigh, mbeta, jatom 
           if (katom .ne. iatom .or. kneigh .ne. ineigh) then
            write (*,*) ' iatom, katom = ', iatom, katom
            write (*,*) ' ineigh, kneigh = ', ineigh, kneigh
            write (*,*) ' Problem in NEIGHBORS, atoms not lined up correctly. '
            write (*,*) ' Fix and restart! '
            stop
           end if 
           neigh_b_vdw(ineigh,iatom) = mbeta
           neigh_j_vdw(ineigh,iatom) = jatom
          end do 
         end do 
         close (unit = 21)
        end if        
 
! Format Statements
! ===========================================================================
 
        return
        end
