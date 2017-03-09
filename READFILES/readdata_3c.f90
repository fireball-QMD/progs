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


! readdata_3c.f90
! Program Description
! ===========================================================================
!       This routine reads the data from the 3-center integral files. When
! read, the information is stored in the array threecint.  This array
! is the field that stores all non-vanishing matrix elements for a general
! 3-center integral.  There are maximal ME3c_max non-vanishing matrix
! elements given on a grid of maximal nfofx data points.  The exact dimensions
! for a given interaction, and a given pair of atoms are numz and num_nonzero.
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
        subroutine readdata_3c (iounit, numx, numy, num_nonzero, isorp,    &
     &                          maxtype, index, ioff3c, nspecies, xintegral)
        use dimensions
        use interactions
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: index
        integer, intent (in) :: ioff3c
        integer, intent (in) :: iounit
        integer, intent (in) :: isorp
        integer, intent (in) :: maxtype
        integer, intent (in) :: nspecies
        integer, intent (in) :: num_nonzero
        integer, intent (in) :: numx, numy 
        
! Output
        real, intent (out), dimension (numXmax, numYmax, ME3c_max,           &
     &                                 0:maxtype, nspecies**3) :: xintegral
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint
        integer integral
        integer jpoint

        real, dimension (ME3c_max, numXmax, numYmax) :: gstore
        real, allocatable, save, dimension(:,:) :: binomial
        integer maxmax
        real, external :: factorial
 
! Procedure
! ===========================================================================
        do jpoint = 1, numy
         do ipoint = 1, numx
          read (iounit,*)                                                    &
     &     (gstore(integral,ipoint,jpoint), integral = 1, num_nonzero)
         end do
        end do
 
        do jpoint = 1, numy
         do ipoint = 1, numx
          do integral = 1, num_nonzero
           xintegral(ipoint,jpoint,integral,isorp,index)                     &
     &      = gstore(integral,ipoint,jpoint) * ioff3c
          end do
         end do
        end do

        return
        end subroutine readdata_3c

