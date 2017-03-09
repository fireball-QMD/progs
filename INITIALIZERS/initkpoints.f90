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


! initkpts.f90
! Program Description
! ===========================================================================
!       This routine initializes the k-points to those read in from the
! k-point file. However, if this is a restarted run, then the k-points
! used are from the file called KPOINTS.
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
        subroutine initkpoints ( )
        use dimensions
        use kpoints
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer ix
        integer num_kpoints
 
        logical kpointfile
 
! Procedure
! ===========================================================================
! If this is a RESTART, then read in the k-points from a file called KPOINT.
        inquire (file = 'KPOINTS', exist = kpointfile)
        if (kpointfile) then
         write (*,*) '  '
         write (*,*) ' We are reading from a KPOINT file. '
         open (unit = 12, file = 'KPOINTS', status = 'old')
         read (12,*) num_kpoints
         if (num_kpoints .ne. nkpoints) then
          write (*,*) ' The k-point file that you are using must not '
          write (*,*) ' belong to the k-point file that you are now '
          write (*,*) ' calculating.  The number of k-points differs '
          write (*,*) ' between the two. '
         end if
         do ikpoint = 1, nkpoints
          read (12,*) (special_k(ix,ikpoint), ix = 1, 3), weight_k(ikpoint)
          if (weight_k(ikpoint) .ne. weight_k_orig(ikpoint)) then
           write (*,*) ' The weights in your KPOINT file do not match '
           write (*,*) ' the weights from the special k-point file. '
           write (*,*) ' Make sure that you are not mixing k-points! '
          end if
         end do
         close (unit = 12)
        else
         do ikpoint = 1, nkpoints
          special_k(:,ikpoint) = special_k_orig(:,ikpoint)
          weight_k(ikpoint) = weight_k_orig(ikpoint)
         end do
        end if
 
 
! Format Statements
! ===========================================================================
 
        return
        end
