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


! getsendrecv.f90
! Program Description
! ===========================================================================
!      
!
! ===========================================================================
! Written by:
! Spencer Shellman
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! ===========================================================================
!
! Program Declaration
! ===========================================================================

subroutine getsendrecv (myrank, reverse, bitmat, send, recv, more)
  use ordern
  implicit none

! getsendrecv is used to reduce communication overhead by using a block matrix to decide which
! processors need to send to/receive from which others.
! With each call, starting with the initial call (init=true), the input block matrix (bitmat)
! is updated to indicate the new global communication state after the required operations are performed.
! Also, each processor gets back send, recv to indicate which processors it should send to
! and/or receive from.

! myrank = rank of this processor (used to select a row and column in bitmat)
! reverse = true if columns indicate processors that must be sent to and rows indicate processors
!   that must be received from, false if vice-versa
! bitmat = nactualprocs x nactualprocs bit matrix;
!   row myrank (or column if reverse=true) has a .true. for each processor that processor myrank
!   should send its packet to.
!   column myrank (or row if reverse=true) has a .true. for each processor that processor myrank should
!   receive a packet from.
! send = on return, the processor that processor myrank should send its
!   packet to after this call, or myrank if it needs to send nothing
! recv = on return, the processor that processor myrank should receive a
!   packet from after this call, or myrank if it needs to receive nothing
! more = on return, .true. if getsendrecv needs to be called again to get
!   info about more processors to send to/receive from,
!   or .false. if no more work needs to be done by processor myrank.
!
! How to use:  The processors involved in a matrix multiplication should set
! up a bit matrix such that row k (or column if reverse=true) contains a .true. for each processor
! that processor k needs to send its packet to, and column k (or row if reverse=true) contains
! a .true. for each processor that processor k needs to receive a packet
! from, to execute the multiplication correctly.  (The contents
! of element (k,k) of the bit matrix are ignored and assumed .false.
! since processor k never needs to send to or receive from itself.)
! The processors involved should call getsendrecv repeatedly, passing in
! the bit matrix each time (this causes the bit matrix to be modified),
! until more = false.  As soon as more = false on a processor, that
! processor can stop calling getsendrecv, since it has no more work to do.
! Each time processor k calls getsendrecv, it gets back information about
! which processor to send to, which processor to receive from, and the
! original source of the packet it will receive.
! On every call (that is, every time processors call getsendrecv together)
! the processors will not be instructed to send more than one section to
! the same processor in one call.  This is how getsendrecv is able to
! instruct a processor to receive from only one processor in one call.
! In addition, getsendrecv zeroes out the entries of the bit matrix
! that correspond to its last output values, so that the next time
! getsendrecv is called it will not instruct a processor to send to
! or receive from a processor that it has already sent to or
! received from.

  integer, intent (in) :: myrank
  logical, intent (in) :: reverse
  logical, dimension (0:nactualprocs-1, 0:nactualprocs-1) , intent (inout) :: bitmat
  integer, intent (out) :: send, recv
  logical, intent (out) :: more

  integer i, j, x, y
  logical anysend, mustsend (0:nactualprocs-1), willrecv (0:nactualprocs-1)

  send = myrank
  recv = myrank
  more = .false.

  if (nactualprocs .gt. 1) then
     anysend = .false.
! Fill in mustsend with flags indicating which processors need to send.
     do i = 0, nactualprocs-1
        mustsend (i) = .false.
        do j = 0, nactualprocs-1
           if (i .ne. j) then
              if (reverse) then
                 x = j
                 y = i
              else
                 x = i
                 y = j
              end if
              if (bitmat (x, y)) then
                 mustsend (i) = .true.
                 anysend = .true.
                 exit
              end if
           end if
        end do
     end do
! If no processor needs to send then the cycle is complete.
     if (anysend) then
        ! FIXME We may want to consider setting more to false as soon as this processor has
        ! no hope of being used again in this cycle.  This would save on busy waiting.
        more = .true.
        willrecv = .false.
        do i = 0, nactualprocs-1
           if (mustsend (i)) then
              j = mod(i+1, nactualprocs)
! We move forward from the location of processor i's packet until we find a processor
! that it must be sent to.
              do
                 if (reverse) then
                    x = j
                    y = i
                 else
                    x = i
                    y = j
                 end if
                 if ((.not. willrecv(j)) .and. bitmat (x, y)) then
! Processor i's packet must be sent to processor j.
! Update the block matrix to indicate that the send is taking place.
                    bitmat (x, y) = .false.
                    willrecv(j) = .true.
! The calling processor is processor i.  It is instructed to send to j.
                    if (myrank .eq. i) send = j
! The calling processor is j.  It is instructed to receive from processor i.
                    if (myrank .eq. j) recv = i
                    exit
                 end if
                 j = mod(j+1, nactualprocs)
                 if (j .eq. i) exit
              end do
           end if
        end do
     end if
  end if
  return
end subroutine getsendrecv





subroutine sendrecv (myrank, iteration, ndiv, nmod, bm, return_bm, &
     & numA, listA, A_compact, nrowsA, ncolsA, nrowsAmax, iAstart, &
     & numB, listB, B_compact, nrowsB, ncolsB, nrowsBmax, iBstart, &
     & numC, listC, C_compact, chooseA, chooseB, countC, ncolsC, nrowsCmax)
  use ordern
  implicit none

! Performs a specified number of distributed matrix multiplications.

! myrank = rank of calling proc
! iteration = 0 if this is being called during the first iteration, when the index sets
!   need to be constructed, > 0 otherwise.
! ndiv = # of rows in A / # of procs
! nmod = # of rows in A mod # of procs
! bm = a nactualprocs X nactualprocs block matrix for A.
!   If return_bm = true then this is generated and returned
!   by sendrecv, otherwise sendrecv uses the input value without modifying it.
! numA,listA = index set of matrices in A (B,C similar)
!   (we assume the same index set applies to all A matrices, or B, or C)
! A_compact = elements of matrices A (B,C similar)
! ncolsA = column dimension of A (B,C similar)
! nrowsA,nrowsAmax = number of rows to operate on in A, and actual row dimension of A
!   (B,C similar; nrowsC is determined implicitly from nrowsA and nrowsB)
! chooseA = array of countC indices, indicating which matrix A should be multiplied to
!   produce the corresponding product matrix C (chooseB similar)
! countC = # of product matrices C that will be generated.
!   The numbers of A and B matrices are fewer than countC.
!   The i'th product matrix C = A(chooseA(i)) * B(chooseB(i)).
! iAstart = currently ignored, should be eliminated
! iBstart = the column in A at which to start each multiplication (that's why it's called iBstart)

  include 'mpif.h'

  integer, intent (in) :: myrank, iteration, ndiv, nmod, iAstart, iBstart, countC
  logical, intent (in) :: return_bm
  logical, dimension (nactualprocs, nactualprocs), intent (inout) :: bm

  integer, intent (in) :: nrowsA, ncolsA, nrowsAmax
  integer, dimension (nrowsAmax), intent (in) :: numA
  integer, dimension (ncolsA, nrowsAmax), intent (in) :: listA
  real, dimension (ncolsA, nrowsAmax, countC), intent (in) :: A_compact

  integer, intent (in) :: nrowsB, ncolsB, nrowsBmax
  integer, dimension (nrowsBmax), intent (in) :: numB
  integer, dimension (ncolsB, nrowsBmax), intent (in) :: listB
  real, dimension (ncolsB, nrowsBmax, countC), intent (in) :: B_compact

  integer, intent (in) :: ncolsC, nrowsCmax
  integer, dimension (nrowsCmax), intent (inout) :: numC
  integer, dimension (ncolsC, nrowsCmax), intent (inout) :: listC
  real, dimension (ncolsC, nrowsCmax, countC), intent (out) :: C_compact
  integer, dimension (countC), intent (in) :: chooseA, chooseB

!
! Local variables
!

  integer countB
  integer packpos_sparse, packsize_sparse
  integer isendproc, irecproc, irecrows, irecstart
  integer request, ierror
  integer imu, inu, index, blockno
  integer sendsize
  logical more_sendrecv

! number of processors that need to be received from
  integer recv_num
! number of processors received from so far
  integer recv_index
! list of processors received from
  integer, dimension (:), allocatable :: recv_ids

! Bit matrix indicating which blocks of a product matrix are nonzero.
  logical, dimension (:, :), allocatable :: bit_matrix
  logical, dimension (:), allocatable :: bit_column

! Temporary storage vectors. 
  integer, dimension (:), allocatable :: numT
  integer, dimension (:, :), allocatable :: listT
  real, dimension (:, :, :), allocatable :: T_compact

  integer*1, dimension (:,:), allocatable :: packarray_sparse

! MPI receive status
  integer, dimension (MPI_STATUS_SIZE) :: status

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! calculate the maximum index into the B matrix array
  countB = 0
  do imu = 1, countC
     countB = max (countB, chooseB (imu))
  end do

! initialize the C index set to 0 if this is the first iteration
  if (iteration .eq. 0) then
     numC = 0
  end if

! The multiplication is performed by multiplying a block of the A matrix by the
! corresponding distributed section of the global B matrix, then adding the
! results together to get the C matrix.
! We receive all the packets first, then multiply, for better scaling.

  if (countB .gt. 0 .and. countC .gt. 0) then

! Send local elements to the next processor.
     if (nactualprocs .gt. 1) then

! Allocate the bit matrix.
        allocate (bit_matrix (nactualprocs, nactualprocs))

! The bit matrix indicates which blocks of the global A matrix are zero.
! If a block of A is zero then we do not need to receive the corresponding B section.
! FIXME should also possibly take into account which local sections of B are zero.
        if (return_bm) then
! FIXME this is an inefficient copy
! Compute the bit vector for A.
           allocate (bit_column (nactualprocs))
           bit_column = .false.
           do imu = 1, nrowsA
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

           call MPI_ALLGATHER (bit_column, nactualprocs, MPI_LOGICAL, bit_matrix, &
                & nactualprocs, MPI_LOGICAL, MPI_BTN_WORLD, ierror)
           deallocate (bit_column)
           bm = bit_matrix
        else
           bit_matrix = bm
        end if

! Determines the number of packets that need to be received
! (nonzero entries in column myrank).
        recv_num = 0
        do imu = 1, nactualprocs
           if (imu.ne.myrank+1 .and. bit_matrix(imu,myrank+1)) recv_num = recv_num + 1
        end do
        allocate (recv_ids(recv_num))

! Allocate some temporary arrays.
        allocate (numT (nrowsBmax))
        allocate (listT (ncolsB, nrowsBmax))
        allocate (T_compact (ncolsB, nrowsBmax, countB))

! Allocates space for packed sparse matrices.
        call sparse_getpacksize (ncolsB, nrowsBmax, packsize_sparse)
        packsize_sparse = countB * packsize_sparse
        allocate (packarray_sparse (packsize_sparse,0:recv_num))

! FIXME We may want to check if these elements need to be sent at all.
        packpos_sparse = 0
        if (countB .gt. 1) then
           call sparse_pack_indices (packsize_sparse, nrowsB, ncolsB, nrowsBmax, &
                & numB, listB, packarray_sparse, packpos_sparse)
           do imu = 1, countB
              call sparse_pack_elements (packsize_sparse, nrowsB, ncolsB, nrowsBmax, &
                   & numB, B_compact(1,1,imu), packarray_sparse, packpos_sparse)
           end do
        else
           call sparse_pack (packsize_sparse, nrowsB, ncolsB, nrowsBmax, &
                & numB, listB, B_compact(1,1,1), packarray_sparse, packpos_sparse)
        end if

! Repeatedly call getsendrecv to get a processor to send to and a processor
! to receive from, until everything has been sent and received.
        sendsize = packpos_sparse
        recv_index = 0
        do
           call getsendrecv (myrank, .false., bit_matrix, &
                & isendproc, irecproc, more_sendrecv)
           if (.not. more_sendrecv) exit

! Send to the processor returned by getsendrecv.
           if (isendproc .ne. myrank) then
              call MPI_ISEND (packarray_sparse, sendsize, MPI_PACKED, &
                   & isendproc, 0, MPI_BTN_WORLD, request, ierror)
           end if
    
! Receive from the processor returned by getsendrecv.
           if (irecproc .ne. myrank) then
              recv_index = recv_index + 1
! save the processor index that this was received from
              recv_ids(recv_index) = irecproc
              call MPI_RECV (packarray_sparse (:,recv_index), packsize_sparse, MPI_PACKED, &
                   & irecproc, 0, MPI_BTN_WORLD, status, ierror)
           end if
! finish the latest send
           if (isendproc .ne. myrank) call MPI_WAIT (request, status, ierror)
        end do
     end if

! For each product the local block of the A matrix is multiplied by the B matrix
! to generate part of the C matrix.
     do imu = 1, countC
        call sparse_mult (nrowsA, ncolsA, ncolsB, ncolsC, nrowsAmax, nrowsB,   &
             & nrowsBmax, nrowsCmax, numA, listA, A_compact(1,1,chooseA(imu)), iAstart,      &
             & numB, listB, B_compact(1,1,chooseB(imu)), iBstart, .false., &
             & imu .eq. 1 .and. iteration .eq. 0, numC, listC, C_compact(1,1,imu))
     end do

     if (nactualprocs .gt. 1) then
! perform the multiplications with the received packets
        do inu = 1, recv_num
           irecproc = recv_ids(inu)
           irecrows = ndiv
           if (irecproc .lt. nmod) then
              irecrows = irecrows + 1
              irecstart = irecrows*irecproc + 1
           else 
              irecstart = (irecrows + 1)*nmod + irecrows*(irecproc - nmod) + 1
           end if

! Now unpack the packet containing the elements that were received.
           packpos_sparse = 0
           if (countB .gt. 1) then
              call sparse_unpack_indices (packarray_sparse (:,inu), packsize_sparse, &
                   & ncolsB, nrowsBmax, numT, listT, packpos_sparse)
              do imu = 1, countB
                 call sparse_unpack_elements (packarray_sparse (:,inu), packsize_sparse, &
                      & ncolsB, nrowsBmax, numT, T_compact(:,:,imu), packpos_sparse)
              end do
           else
              call sparse_unpack (packarray_sparse (:,inu), packsize_sparse, &
                   & ncolsB, nrowsBmax, numT, listT, T_compact, packpos_sparse)
           end if

           do imu = 1, countC
              call sparse_mult (nrowsA, ncolsA, ncolsB, ncolsC, nrowsAmax, irecrows,   &
                   & nrowsBmax, nrowsCmax, numA, listA, A_compact(1,1,chooseA(imu)), iAstart,      &
                   & numT, listT, T_compact(1,1,chooseB(imu)), irecstart, .true., &
                   & imu .eq. 1 .and. iteration .eq. 0, numC, listC, C_compact(1,1,imu))
           end do
        end do

        deallocate (numT, listT, T_compact)
        deallocate (packarray_sparse)
        deallocate (bit_matrix)
        deallocate (recv_ids)
     end if
  end if
  return
end subroutine sendrecv
