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


! sparse_getdimension.f90
! Program Description
! ===========================================================================
!       Before performing a matrix multiplication, this routine determines the
! dimension of the resultant product.  The dimension obtained is then used to
! allocate the array with the appropriate dimensions. 
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
        subroutine sparse_getdimension (nprows, nAmax, nBmax, norbitals,     &
     &                                  ncolsAmax, ncolsB, ncolsBmax,        &
     &                                  numA, listA, istartrowA, numB, listB,&
     &                                  istartrowB, append, numC, listC, nCmax)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: istartrowA
        integer, intent (in) :: istartrowB
        integer, intent (in) :: nprows
        integer, intent (in) :: nAmax
        integer, intent (in) :: nBmax
        integer, intent (in) :: ncolsAmax
        integer, intent (in) :: ncolsB
        integer, intent (in) :: ncolsBmax
        integer, intent (in) :: norbitals
 
        integer, intent (in), dimension (ncolsAmax) :: numA
        integer, intent (in), dimension (ncolsBmax) :: numB
        integer, intent (inout), dimension (nprows) :: numC
        integer, intent (in), dimension (nAmax, ncolsAmax) :: listA
        integer, intent (in), dimension (nBmax, ncolsBmax) :: listB
        integer, intent (inout), dimension (norbitals, nprows) :: listC

        logical, intent (in) :: append

! Output
        integer, intent (inout) :: nCmax

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iA, iB
        integer imu
        integer index1, index2
        integer inu
        integer mmu

        integer, dimension (:), allocatable :: indexloc
        integer, dimension (:), allocatable :: indexv

! Allocate Arrays
! ===========================================================================
        allocate (indexloc (norbitals))
        allocate (indexv (norbitals))
 
! Procedure
! ===========================================================================
! Initialize the resultant matrix to zero.
        if (.not. append) then
         numC = 0 
         listC = 0
         nCmax = -99
        end if

! Set up indices for C sparse matrix. Determine maximum number of indices.
        indexloc = 0
        indexv = 0
        do iA = 1, nprows
         index2 = numC(iA)
         if (append) then              ! If we are appending, then do not
          do imu = 1, index2           ! count this index in listC again.
           index1 = listC(imu,iA)
           indexloc(index1) = 1
           indexv(imu) = index1
          end do
         end if
         do imu = 1, numA(iA)
          iB = listA(imu,iA)
          if (iB .ge. istartrowB .and. iB .le. ncolsB + istartrowB - 1) then
           iB = iB - istartrowB + 1
           do inu = 1, numB(iB)
            index1 = listB(inu,iB)
            if (indexloc(index1) .eq. 0) then
             indexloc(index1) = 1
             index2 = index2 + 1
             indexv(index2) = index1
            end if
           end do
          end if
         end do
         numC(iA) = index2 
         if (index2 .gt. nCmax) nCmax = index2 
         do imu = 1, index2
          index1 = indexv(imu)
          indexv(imu) = 0
          indexloc(index1) = 0
          listC(imu,iA) = index1
         end do
        end do


! Deallocate Arrays
! ===========================================================================
        deallocate (indexloc)
        deallocate (indexv)
 
! Format Statements
! ===========================================================================
 
        return
        end
