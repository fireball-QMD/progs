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


! find_neigh_max.f90
! Program Description
! ===========================================================================
!       Finds all the maximum number of neighbors to atoms in the central cell.
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
        subroutine find_neigh_max (nprocs, my_proc, iordern, icluster, ivdw)
        use configuration  
        use dimensions
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: iordern
        integer, intent (in) :: ivdw
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs



!$ volatile rcutoff, icluster, ivdw, natoms, nprocs, my_proc, iordern
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer imu
        integer in1, in2
        integer jatom
        integer mbeta
        integer natomsp
        integer num_neigh
        integer num_neigh_vdw

        real distance2
        real range2
        real rcutoff_i
        real rcutoff_j

! Procedure
! ===========================================================================
!        if (my_proc .eq. 0 .and. wrtout)                                                  &
!     &   write (*,*) ' Determine maximum number of neighbors. ' 

        if (icluster .eq. 1) mbeta_max = 0

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

! Loop over all atoms.
        neigh_max = -99
        if (ivdw .eq. 1) neigh_max_vdw = -99
!$omp parallel do private (num_neigh, num_neigh_vdw, rcutoff_i, rcutoff_j)   &
!$omp&            private (in1, in2, imu, distance2, range2)
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
           distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
     &                + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
     &                + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2

! Add a small displacement to the sum of cutoffs. 
           range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2

           if (distance2 .le. range2) then
            num_neigh = num_neigh + 1
           end if
           
           if (ivdw .eq. 1 .and. distance2 .le. range_vdw**2) then 
            num_neigh_vdw = num_neigh_vdw + 1
           end if
          end do
         end do

! Maximum number of neighbors thus far.
! FIXME: I tried using omp atomic here and it crashed.
!$omp critical
         neigh_max = max(neigh_max, num_neigh)
         if (ivdw .eq. 1) neigh_max_vdw = max(neigh_max_vdw, num_neigh_vdw)
!$omp end critical
        end do
        if (iordern .eq. 1) call find_neigh_max_ordern_final()

! Format Statements
! ===========================================================================
!write(*,*)'find_neigh_max(): neigh_max=', neigh_max 
        return
      end subroutine find_neigh_max
