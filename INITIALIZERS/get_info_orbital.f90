! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young of Utah - James P. Lewis, Chair
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


! make_munu.f90
! Program Description
! ===========================================================================
! ===========================================================================
! Original code written by Jose Ortega & Daniel Gonzalez.
 
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
        subroutine get_info_orbital (natoms)
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer issh
        integer iorb

! Allocate Arrays
! ===========================================================================
        allocate (getiatom(norbitals))
        allocate (getissh(norbitals))
        allocate (getlssh(norbitals))
        allocate (getmssh(norbitals))
 
! Procedure
! ===========================================================================
        imu=0
        do iatom=1,natoms                              !example iatom = 12 with Z=14 ans s,p,d
         do issh=1,nssh(imass(iatom))                  ! nssh = 3  (s s+ p) 
          do iorb = 1, 2*lssh(issh,imass(iatom))+1     ! ilssh = 1, [0->1,1->3,2->5] 
           imu=imu+1
           getmssh(imu)=iorb   
           getissh(imu)=issh
           getlssh(imu)=lssh(issh,imass(iatom))    
           getiatom(imu)=iatom
          end do
         end do
        end do 

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine get_info_orbital
