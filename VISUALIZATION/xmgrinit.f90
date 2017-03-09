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


! xmgrinit.f90
! Program Description
! ===========================================================================
!       This routine opens and initializes the xmgr program so that 
! data can be plotted 'on the fly'.
!
! ===========================================================================
! Original code written by Aaron Lefohn.

! Code rewritten by:
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
        subroutine xmgrinit (time, etotper, getotper)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: etotper
        real, intent (in) :: getotper
        real, intent (in) :: time
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ACEgrOpenf
        integer isOpen

        character (len=64) buffer

        external ACEgrOpenf
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Start xmgr with a buffer size of 512 and open the pipe
        isOpen = ACEgrOpenf (1024)

! Initialize the graph legend in xmgr
        if (isOpen .eq. 0) then
         call ACEgrCommandf ('xaxis label "Time (fs)"')
         call ACEgrCommandf ('yaxis label "Energy (eV/atom)"')
         call ACEgrCommandf ('legend on')
         call ACEgrCommandf ('legend string 0 "etotper"')
         call ACEgrCommandf ('legend string 1 "getotper"')

         write (buffer, '("g0.s0 point ", f10.4, ",", f12.4)') time, etotper
         call ACEgrCommandf (buffer)

         write (buffer, '("g0.s1 point ", f10.4, ",", f12.4)') time, getotper
         call ACEgrCommandf (buffer)

         call ACEgrCommandf ('autoscale')
        else
         write (*,*) ' Cannot run xmgr! Mission aborted! '
         write (*,*)
         stop 
        end if

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
