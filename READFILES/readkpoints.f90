! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

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


! readkpoints.f90
! Program Description
! ===========================================================================
!       This routine will read in the special k-points from the *.kpts file.
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
        subroutine readkpoints (kptsfile, nkpoints, icluster)
        use dimensions
        use kpoints
        use options, only : verbosity
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        character (len=40), intent (in) :: kptsfile
        integer, intent(in) :: icluster
 
! Output
        integer, intent (out) :: nkpoints
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
 
        real sum_weight
 
! Procedure
! ===========================================================================
! Open the file containing the special k_points.
        if (icluster .ne. 0) then
          write (*,*) '  '
          write (*,*) ' You are doing a cluster calculation '
          write (*,*) ' k-points are not read in '
          write (*,*) '  '
          nkpoints = 1
          allocate (special_k(3, nkpoints))
          allocate (special_k_orig(3, nkpoints))
          allocate (scale_k(3, nkpoints))
          allocate (weight_k(nkpoints))
          allocate (weight_k_orig(nkpoints))
          special_k_orig(:,1) = 0
          weight_k_orig(1) = 1
          return
        end if

        open (unit = 54, file = kptsfile, status = 'old')
        read (54,*) nkpoints
 
        if (verbosity .ge. 3) write (*,*) '  '
        write (*,*) '  Reading k-points '
        if (verbosity .ge. 3) write (*,100)
        if (verbosity .ge. 3) write (*,*) ' nkpoints = ', nkpoints
 
        if (nkpoints .le. 0) then
         write (*,*) ' nkpoints .le. 0 Huh! Fix XXX.kpts file! '
         stop
        end if

! Allocate the kpoints module memory
        allocate (special_k(3, nkpoints))
        allocate (special_k_orig(3, nkpoints))
        allocate (scale_k(3, nkpoints))
        allocate (weight_k(nkpoints))
        allocate (weight_k_orig(nkpoints))

        sum_weight = 0.0d0
        do ikpoint = 1, nkpoints
         read (54,*) special_k_orig(:,ikpoint), weight_k_orig(ikpoint)
         if (verbosity .ge. 3) write (*,300) ikpoint, special_k_orig(:,ikpoint),                   &
     &                 weight_k_orig(ikpoint)
         sum_weight = sum_weight + weight_k_orig(ikpoint)
        end do
        if (verbosity .ge. 3) write (*,100)
 
        if (abs(sum_weight - 1.0d0) .gt. 1.0d-3) then
         write (*,*) ' Sum of k-point weights = ', sum_weight
         write (*,*) ' They do not add to one (within 1.0d-3) '
         write (*,*) ' Sorry -- Stop --- fix XXX.kpts'
         write (*,*) ' The k-points sum of weights does not add to 1! '
         stop
        end if
 
        close (unit = 54)
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
300     format (2x, i4, ' - ikpoint = ', 3f11.6, 2x, ' weight = ', f9.4)
 
        return
        end
