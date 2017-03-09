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


! set_maxdimensions.f90
! Program Description
! ===========================================================================
!       This routine determines the maximum dimension for a given parameter.
! This maximum is the ultimate maximum having considered the information 
! from all processors.  Information for a set of indices must be communicated 
! to the processor myrank.  This information is then used locally to 
! evaluate a local maximum.  The ultimate maximum is the maximum of all
! local maximums.
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
        subroutine set_maxdimension (nactualprocs, myrank, nprows, nAmax,    &
     &                               nBmax, ncolsAmax, ncolsB, ncolsBmax,    &
     &                               numA, listA, istartrowA, numB, listB,   &
     &                               istartrowB, ioption, nCmax)
        use interactions
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ioption
        integer, intent (in) :: istartrowA
        integer, intent (in) :: istartrowB
        integer, intent (in) :: myrank
        integer, intent (in) :: nactualprocs
        integer, intent (in) :: nprows
        integer, intent (in) :: nAmax
        integer, intent (in) :: nBmax
        integer, intent (in) :: ncolsAmax
        integer, intent (in) :: ncolsB
        integer, intent (in) :: ncolsBmax

        integer, intent (in), dimension (ncolsAmax) :: numA
        integer, intent (in), dimension (ncolsBmax) :: numB
        integer, intent (in), dimension (nAmax, ncolsAmax) :: listA
        integer, intent (in), dimension (nBmax, ncolsBmax) :: listB

! Output
        integer, intent (out) :: nCmax

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror, request
        integer imu, inu
        integer index, blockno
        integer iproc
        integer irecrows
        integer irecstart
        integer packpos_sparse
        integer packsize_sparse
        integer ndiv, nmod
        integer isendproc, irecproc
        logical more_sendrecv
        integer sendsize

        integer*1, dimension (:,:), allocatable :: packarray_sparse
        integer, dimension (MPI_STATUS_SIZE) :: status

        integer, dimension (:),allocatable :: numC
        integer, dimension (:,:),allocatable :: listC

! Temporary array communicated from neighboring processor.
        integer, dimension (:),allocatable :: numT
        integer, dimension (:,:),allocatable :: listT

! Bit matrix indicating which blocks of a product matrix are nonzero.
        logical, dimension (:, :), allocatable :: bit_matrix
        logical, dimension (:), allocatable :: bit_column

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Processor 0:
        allocate(numC(nprows))
        allocate(listC(norbitals,nprows))
        allocate(numT(ncolsBmax))
        allocate(listT(nBmax,ncolsBmax))
        numC = 0
        call sparse_getdimension (nprows, nAmax, nBmax, norbitals, ncolsAmax,&
     &                            ncolsB, ncolsBmax, numA, listA, istartrowA,&
     &                            numB, listB, istartrowB, .false., numC,    &
     &                            listC, nCmax)

! Local section sizes
        if (ioption .eq. 0) then
         ndiv = norbitals / nactualprocs
         nmod = mod(norbitals,nactualprocs)
        else
         ndiv = nbands / nactualprocs
         nmod = mod(nbands,nactualprocs)
        end if

! Send local coefficients to the next processor.
        if (nactualprocs .gt. 1) then

! Allocate the bit matrix.
! FIXME: Can we reuse this bit matrix in the code following the 
! set_maxdimension call?
         allocate (bit_matrix (nactualprocs, nactualprocs))
         allocate (bit_column (nactualprocs))
         bit_column = .false.
         do imu = 1, nprows
          do inu = 1, numA(imu)
           index = listA(inu,imu)
           if (index .gt. (ndiv+1) * nmod) then
            blockno = (index - (ndiv+1) * nmod - 1) / ndiv + nmod + 1
           else
            blockno = (index - 1) / (ndiv + 1) + 1
           end if
           bit_column (blockno) = .true.
          end do
         end do

! Communicate the bit vector over all processors.

         call MPI_ALLGATHER (bit_column, nactualprocs, MPI_LOGICAL,          &
     &                       bit_matrix, nactualprocs, MPI_LOGICAL,          &
     &                       MPI_BTN_WORLD, ierror)

         deallocate (bit_column)

! Allocates space for packed sparse matrices.
         call sparse_getpacksize (nBmax, ncolsBmax, packsize_sparse)
         allocate (packarray_sparse (packsize_sparse,2))

! FIXME We may want to check if these elements need to be sent at all.
         packpos_sparse = 0
         call sparse_pack_indices (packsize_sparse, ncolsB, nBmax, ncolsBmax,&
     &                             numB, listB, packarray_sparse,            &
     &                             packpos_sparse)

! Now determine nCmax based on each of the products from the other processors. 
! Thus the value of nCmax may increase as more products are considered.

! Repeatedly call getsendrecv to get a processor to send to and a processor
! to receive from, until everything has been sent and received.
         sendsize = packpos_sparse
         do
          call getsendrecv (myrank, .false., bit_matrix,      &
     &                      isendproc, irecproc, more_sendrecv)
          if (.not. more_sendrecv) exit

! Send to the processor returned by getsendrecv.
          if (isendproc .ne. myrank) then
           call MPI_ISEND (packarray_sparse, sendsize, MPI_PACKED, &
     &                     isendproc, 0, MPI_BTN_WORLD, request, ierror)
          end if
    
! Receive from the processor returned by getsendrecv.
          if (irecproc .ne. myrank) then
           call MPI_RECV (packarray_sparse (:,2), packsize_sparse,       &
     &                    MPI_PACKED, irecproc, 0, MPI_BTN_WORLD, status,    &
     &                    ierror)
          end if

          if (isendproc .ne. myrank) call MPI_WAIT (request, status, ierror)

          if (irecproc .ne. myrank) then
             irecrows = ndiv
             if (irecproc .lt. nmod) then
                irecrows = irecrows + 1
                irecstart = irecrows*irecproc + 1
             else
                irecstart = (irecrows + 1)*nmod + irecrows*(irecproc - nmod) + 1
             end if

! Now unpack the packet containing the indices that were received.
             packpos_sparse = 0
             call sparse_unpack_indices (packarray_sparse (:,2),           &
                  &                                 packsize_sparse, nBmax, ncolsBmax,    &
               &                                 numT, listT, packpos_sparse)
             call sparse_getdimension (nprows, nAmax, nBmax, norbitals,        &
                  &                               ncolsAmax, irecrows, ncolsBmax, numA,   & 
                  &                               listA, istartrowA, numT, listT,         &
                  &                               irecstart, .true., numC, listC, nCmax)
          end if
         end do
         deallocate (packarray_sparse)
         deallocate (bit_matrix)
        end if
        
 
! Deallocate Arrays
! ===========================================================================
        deallocate(numC,listC,numT,listT)
! Format Statements
! ===========================================================================
 
        return
        end
