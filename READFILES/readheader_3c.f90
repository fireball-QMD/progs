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


! readheader_3c.f90
! Program Description
! ===========================================================================
!       Read the header of the 3-center data files.
!
! ===========================================================================
! Code rewritten by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readheader_3c (iounit, numx, numy, xmax, ymax)
        use dimensions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iounit
 
! Output
! ymax, numy: bond charge distances grid
! xmax, numx: neutral atom distances grid
        integer, intent (out) :: numx, numy
        real, intent (out) :: xmax, ymax
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iline
        integer nucZ1, nucZ2, nucZ3, nr, ntheta_in, nphi2
        real rc1a, rc2a, rc3a
 
        character (len = 70) message
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        do iline = 1, 10
         read (iounit,100) message
        end do
 
        read (iounit,*) nphi2, nr, ntheta_in
 
        read (iounit,*) ymax, numy
        read (iounit,*) xmax, numx
 
        read (iounit,100) message
 
        read (iounit,*) nucZ1, rc1a
        read (iounit,*) nucZ2, rc2a
        read (iounit,*) nucZ3, rc3a
 
        read (iounit,100) message
 
        if (numx .gt. numXmax .or. numy .gt. numYmax) then
         write (*,*) ' Courseness too fine in 3c data files. '
         write (*,*) ' numx = ', numx, ' numXmax = ', numXmax
         write (*,*) ' numy = ', numy, ' numYmax = ', numYmax
         write (*,*) ' Change numXmax and numYmax in MODULES/dimensions.f90! '
         stop
        end if
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (a70)
 
        return
        end
