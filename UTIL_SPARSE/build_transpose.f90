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

! build_transpose.f90
! Program Description
! ===========================================================================
!       This routine builds the transpose of a matrix. To do this requires
! communicating pieces of a matrix to the different processors in order to
! obtain all the necessary information for building the transpose.
!
!       This routine was origanally written for determining the transpose of 
! the c matrix, but can be used for determining the tranpose of any matrix. 
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
        subroutine build_transpose (nactualprocs, myrank, bm, return_bm,     &
     &                              ncmax, nctmax,     &
     &                              nprows, nprowsmax, ncrows, ncrowsmax,    &
     &                              numc_local, listc_local, c_compact_local,&
     &                              ipstart, numct_local, listct_local,      &
     &                              ct_compact_local, icstart, rows, cols)
        use interactions
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icstart
        integer, intent (in) :: ipstart
        integer, intent (in) :: myrank
        integer, intent (in) :: nactualprocs
        integer, intent (in) :: ncmax
        integer, intent (in) :: nctmax
        integer, intent (in) :: ncrows
        integer, intent (in) :: ncrowsmax
        integer, intent (in) :: nprows
        integer, intent (in) :: nprowsmax
        integer, intent (in) :: rows, cols
        logical, intent (in) :: return_bm

        integer, intent (in), dimension (nprowsmax) :: numc_local
        integer, intent (in), dimension (ncmax, nprowsmax) :: listc_local
        logical, intent (inout), dimension (nactualprocs, nactualprocs) :: bm

        real, intent (in), dimension (ncmax, nprowsmax) :: c_compact_local

! Output
        integer, intent (out), dimension (ncrowsmax) :: numct_local
        integer, intent (out), dimension (nctmax, ncrowsmax) :: listct_local

        real, intent (out), dimension (nctmax, ncrowsmax) :: ct_compact_local
!$ volatile numct_local,listct_local,ct_compact_local
  
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer request, ierror
        integer imu, inu
        integer index
        integer iproc
        integer irecproc, isendproc
        integer irecrows
        integer irecstart
        integer jndex
        integer kndex
        integer list_temp
        integer mmu
        integer nnu
        integer packpos_sparse
        integer packsize_sparse
        integer nrdiv, nrmod, ncdiv, ncmod
        integer sndcount, blockno
        logical more_sendrecv
        integer sendsize

        integer*1, dimension (:, :), allocatable :: packarray_sparse

! Bit matrix and column for communication reduction
        logical, dimension (nactualprocs, nactualprocs) :: bit_matrix
        logical, dimension (nactualprocs) :: bit_column

! Temporary storage vectors. 
        integer, dimension (:),allocatable :: numT1_local
        integer, dimension (:,:),allocatable :: listT1_local

        real, dimension (:,:),allocatable :: T1_compact_local

! MPI receive status
        integer, dimension (MPI_STATUS_SIZE) :: status

! number of processors that need to be received from
        integer recv_num
! number of processors received from so far
        integer recv_index
! list of processors received from
        integer, dimension (:), allocatable :: recv_ids

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ============================================================================
 
! Procedure
! ============================================================================
        numct_local = 0

! Now pass the local piece of the c matrix to the other processors in order
! to build the remaining piece of the c transpose matrix. 
! Send local coefficients to the next processor.
        if (nactualprocs .gt. 1) then

         allocate(numT1_local(nprowsmax))
         allocate(listT1_local(ncmax,nprowsmax))
         allocate(T1_compact_local(ncmax,nprowsmax))

! Local section sizes
         nrdiv = rows / nactualprocs
         nrmod = mod(rows,nactualprocs)
         ncdiv = cols / nactualprocs
         ncmod = mod(cols,nactualprocs)

         if (return_bm) then
! Compute the bit vector for this section of C.
            bit_column = .false.
            do imu = 1, nprows
               do inu = 1, numc_local(imu)
                  index = listc_local(inu,imu)
                  if (index .gt. (ncdiv+1) * ncmod) then
                     blockno = (index - (ncdiv+1) * ncmod - 1) / ncdiv + ncmod + 1
                  else
                     blockno = (index - 1) / (ncdiv + 1) + 1
                  end if
                  bit_column (blockno) = .true.
               end do
            end do

! Communicate the bit vector over all processors.

            call MPI_ALLGATHER (bit_column, nactualprocs, MPI_LOGICAL,          &
                 &                       bit_matrix, nactualprocs, MPI_LOGICAL,          &
                 &                       MPI_BTN_WORLD, ierror)
            bm = bit_matrix
         else
            bit_matrix = bm
         end if


! Determines the number of packets that need to be received
! (nonzero entries in row myrank).
         recv_num = 0
         do imu = 1, nactualprocs
          if (imu.ne.myrank+1 .and. bit_matrix(myrank+1,imu)) recv_num = recv_num + 1
         end do
         allocate (recv_ids(recv_num))

! Allocate space for packed sparse matrices.
         call sparse_getpacksize (ncmax, nprowsmax, packsize_sparse)
         allocate (packarray_sparse (packsize_sparse,0:recv_num))
         packpos_sparse = 0
! FIXME We may want to check if these coefficients need to be sent at all.
         call sparse_pack (packsize_sparse, nprows, ncmax, nprowsmax,        &
     &                     numc_local, listc_local, c_compact_local,         &
     &                     packarray_sparse, packpos_sparse)

! Repeatedly call getsendrecv to get a processor to send to and a processor
! to receive from, until everything has been sent and received.
         sendsize = packpos_sparse
         recv_index = 0
         do
          call getsendrecv (myrank, .true., bit_matrix,       &
     &                      isendproc, irecproc, more_sendrecv)
          if (.not. more_sendrecv) exit

! Send to the processor returned by getsendrecv.
          if (isendproc .ne. myrank) then
           call MPI_ISEND (packarray_sparse, sendsize, MPI_PACKED, &
     &                     isendproc, 0, MPI_BTN_WORLD, request, ierror)
          end if
    
! Receive from the processor returned by getsendrecv.
          if (irecproc .ne. myrank) then
             recv_index = recv_index + 1
! save the processor index that this was received from
             recv_ids(recv_index) = irecproc
             call MPI_RECV (packarray_sparse (:,recv_index), packsize_sparse,       &
                  &                    MPI_PACKED, irecproc, 0, MPI_BTN_WORLD, status,    &
                  &                    ierror)
          end if
          if (isendproc .ne. myrank) call MPI_WAIT (request, status, ierror)
         end do
        end if

! Calculate c transpose for the local c matrix contained on this processor.
        call sparse_copy (ncmax, nctmax, nprows, nprowsmax, ncrows,          &
     &                    ncrowsmax, numc_local, listc_local,                &
     &                    c_compact_local, ipstart, icstart, 0, .false.,     &
     &                    .true., 1.0, numct_local, listct_local,            &
     &                    ct_compact_local)


        if (nactualprocs .gt. 1) then
! perform the transposes with the received packets
         do inu = 1, recv_num
          irecproc = recv_ids(inu)
          irecrows = nrdiv
          if (irecproc .lt. nrmod) then
             irecrows = irecrows + 1
             irecstart = irecrows*irecproc + 1
          else 
             irecstart = (irecrows + 1)*nrmod + irecrows*(irecproc - nrmod) + 1
          end if

! Now unpack the packet containing the coeffiecients that was received.
          packpos_sparse = 0
          call sparse_unpack (packarray_sparse (:,inu), packsize_sparse,  &
     &                         ncmax, nprowsmax, numT1_local, listT1_local,  &
     &                         T1_compact_local, packpos_sparse)
          call sparse_copy (ncmax, nctmax, nprows, nprowsmax, ncrows,       &
     &                       ncrowsmax, numT1_local, listT1_local,           &
     &                       T1_compact_local, irecstart, icstart, 0, .true.,&
     &                       .true., 1.0, numct_local, listct_local,         &
     &                       ct_compact_local)
         end do

         deallocate (packarray_sparse,recv_ids,numT1_local,listT1_local,T1_compact_local)
        end if

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
 
        return
        end
