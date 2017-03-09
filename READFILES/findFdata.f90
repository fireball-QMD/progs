! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
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


! fireball.f90
! Program Description
! ===========================================================================
!       This finds the Fdata if it does not exist in cwd.
!
! ===========================================================================
! Code written by:
! J.B. Keith
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine findFdata (fdataLocation)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        character(len = 200), intent(out) :: fdataLocation
        logical qexist

! Procedure
! ===========================================================================

        inquire(file='Fdata/info.dat',exist=qexist)
        if(qexist) then
          fdataLocation='Fdata'
        else
          open(unit=412,file='Fdata.optional',status='unknown',action='read')
          read(412,'(a)')fdataLocation
          close(412)
        end if  

        end subroutine
