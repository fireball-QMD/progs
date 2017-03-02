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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! sparse_copy.f90
! Program Description
! ===========================================================================
!       This routine copies a segment of size 'rows' x 'cols', starting at row 
! index 'startrow', from sparse matrix A into sparse matrix B, multiplying each
! element by 'mult' before copying.  (NOTE: The routine will NOT warn you if 
! you pass in mult = 0.)  Copying from a sparse matrix to itself should be 
! safe as long as 'append' is true.
!
! ===========================================================================
! Original code written by Spencer Shellman.

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
        subroutine sparse_copy (nAmax, nBmax, ncolsA, ncolsAmax, ncolsB,     &
     &                          ncolsBmax, numA, listA, elementsA,           &
     &                          istartrowA, istartrowB, isparsesize, append, &
     &                          transpose, mult, numB, listB, elementsB)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: isparsesize
        integer, intent (in) :: istartrowA
        integer, intent (in) :: istartrowB
        integer, intent (in) :: nAmax
        integer, intent (in) :: nBmax
        integer, intent (in) :: ncolsA
        integer, intent (in) :: ncolsAmax
        integer, intent (in) :: ncolsB
        integer, intent (in) :: ncolsBmax

        integer, intent (in), dimension (ncolsAmax) :: numA
        integer, intent (in), dimension (nAmax, ncolsAmax) :: listA

        real, intent (in) :: mult

        real, intent (in), dimension (nAmax, ncolsAmax) :: elementsA
 
! append = true if the source matrix should be appended to the sparse matrix
! already contained in 'numB'/'listB'/'elementsB' after the last row, false 
! otherwise
! transpose = true if the transpose of the sparse matrix is to be copied, false
! if the untransposed matrix is to be copied
        logical, intent (in) :: append
        logical, intent (in) :: transpose

! Output
        integer, intent (out), dimension (ncolsBmax) :: numB
        integer, intent (out), dimension (nBmax, ncolsBmax) :: listB

        real, intent (out), dimension (nBmax, ncolsBmax) :: elementsB

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu, inu
        integer index
        integer jmu
        integer mmu
        integer ncdim
        integer nrdim
        integer nnu

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize matrix B to zero.
        if (.not. append) then
         numB = 0
         listB = 0
         elementsB = 0.0d0
        end if

! FIXME should we OpenMP-ize the loops?
! Do a straight copy if the dimensions are the same, else copy (or append)   
! segment desired into new array. 
        if (nAmax .eq. nBmax .and. ncolsA .eq. ncolsB .and. (.not. append)   &
     &                       .and. (.not. transpose)) then
         numB(1:ncolsB) = numA(1:ncolsA)
         do imu = 1, ncolsB
          do inu = 1, numB(imu)
           listB(inu,imu) = listA(inu,imu)
           elementsB(inu,imu) = mult*elementsA(inu,imu)
          end do
         end do
        else
         do imu = 1, ncolsA
          do jmu = 1, numA(imu)
           index = listA(jmu,imu)
           if (index .gt. istartrowB + ncolsB - 1) then
            exit
           else if (index .ge. istartrowB) then
            index = index + 1 - istartrowB
            if (transpose) then
             mmu = imu
             nnu = index
            else
             mmu = index
             nnu = imu
            end if
            numB(nnu) = numB(nnu) + 1
            if (append) then
             listB(numB(nnu),nnu) = mmu + isparsesize + istartrowA - 1
            else
             listB(numB(nnu),nnu) = mmu + istartrowA - 1
            end if
            elementsB(numB(nnu),nnu) = mult*elementsA(jmu,imu)
           end if
          end do
         end do
        end if

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
