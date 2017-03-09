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


! dohermit.f90
! Program Description
! ===========================================================================
! make all matrices hermitian, it means we force those matrices to be symetric
!
! ===========================================================================
! Code written by:

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine dohermit (nprocs, my_proc, iordern)
        use configuration
        use dimensions
        use interactions
        use neighbor_map
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iordern
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer inu
        integer in1
        integer in2
        integer ineigh
        integer jatom
        integer jbeta
        integer jneigh
        integer katom
        integer kbeta
        integer kneigh
        integer iatomstart, natomsp

        real dr
        real rmax
        real mat

        real, dimension (3) :: ri
        real, dimension (3) :: rj
        real, dimension (3) :: rk
        real, dimension (3) :: rij
        real, dimension (3) :: rkj

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
        if (my_proc .eq. 0 .and. wrtout) then
         write (*,*) ' Welcome to dohermit subroutine '
        end if

! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
                      + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if

! Loop over all atoms
        do iatom = iatomstart, iatomstart - 1 + natomsp

         ri(:) = ratom(:,iatom)

! Loop over all known neighbors of iatom. Call these atoms ineigh.
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          jbeta = neigh_b(ineigh,iatom)


! Keep only third party common neighbors.
          if (.not. (iatom .eq. jatom .and. jbeta .eq. 0)) then
           in1 = imass(iatom)
           rj (:) = ratom(:,jatom)
! distance iatom-jatom
           rij(:) = ri(:) - rj(:)

! Loop over all known neighbors of atom jatom, and check to see if atoms
! iatom and jatom are neighbors.
           do kneigh = 1, neighn(jatom)
            katom = neigh_j(jneigh,jatom)
            kbeta = neigh_b(jneigh,jatom)

! Keep only third party common neighbors.
            if (.not. (jatom .eq. katom .and. kbeta .eq. 0)                   &
     &          .and. (ineigh .ne. kneigh)) then
             in2 = imass(jatom)
             rk(:) = ratom(:,katom)
             rkj(:) = rk(:) - rj(:)

             rmax = 0.0d0
             do i = 1,3
              dr = abs(rij(i) - rkj(i))
              if (dr .gt. rmax) rmax = dr
             enddo

! check if rij & rkj are equivalent
             if (rmax .le. range2) then

! Found one!, symetrize matrix
              do imu = 1, num_orb(in1)
               do inu = 1, num_orb(in2)
                mat = 0.5d0*(h_mat(imu,inu,ineigh,iatom)             &
      &           + h_mat(inu,imu,kneigh,jatom))
                h_mat(imu,inu,ineigh,iatom) = mat
                h_mat(inu,imu,kneigh,jatom) = mat
               enddo ! do inu
              enddo ! do imu
             end if ! if(rmax)
            end if ! if
           end do ! do kneigh
          end if ! if
         end do ! do ineigh
        end do ! do iatom


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end
