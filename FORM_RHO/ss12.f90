! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
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

 
! ss12.f90
! Program Description
! ===========================================================================
! 	This routine calculates S^1/2, where S is the overlap matrix of the
! strong covalent interactions. To keep within the spirit of Order N, this
! calculation is done by taking and expanding S^1/2 in a Chebyshev series.
! 
!                   f(x) ~= sum [c_k*T_(k-1)(x)] - 0.5*c_1 
!
! Chebyshev polynomials are bound between [-1,1]. So we must scale the 
! eigenvalues of the overlap to fit into this range, thus [a,b] maps onto
! [-1,1].
!	Also to make things converge nicely, the range of the sparse form
! must be increased. This is because S^1/2 has more interactions between 
! neighbors than does S.
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
        subroutine ss12 (nprows, ipstart, ipower)
        use interactions
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
	integer, intent(in) :: ipower	     ! the power of s calculating +-1
	integer, intent(in) :: ipstart       
	integer, intent(in) :: nprows        
 
! Local Parameters and Data Declaration
! ===========================================================================
	real, parameter :: a = 0.10d0        ! lower bound of range
	real, parameter :: b = 4.20d0        ! upper bound of range
	integer, parameter :: nchebs = 50    ! number Cheby coefficients
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror
	integer imu                     ! counter over basis
	integer icheb                   ! degree of Chebyshev polynomial
        integer ichooseA
        integer ichooseB
        integer inu, jnu                ! counters over indices
        integer myrank
	integer nS12max_local 
        integer nodiv, nomod
	integer nterms                  ! number of Chebyshev terms to take

	integer, dimension (:), allocatable :: index1 

	real, external :: func12        ! function for s^1/2
	real, external :: funcm12       ! function for s^-1/2

        integer, dimension (:), allocatable :: numt2_local
        integer, dimension (:, :), allocatable :: listt2_local

	real, dimension (0:nchebs - 1) :: c            ! Chebyshev coefficients
	real, dimension (:,:), allocatable :: delta
	real, dimension (:,:), allocatable :: t0 
	real, dimension (:,:), allocatable :: t1
	real, dimension (:,:), allocatable :: t2

        character (len = 40) append_string
        character (len = 10) extension
        character (len = 40) filename
        character (len = 40) root

        logical, dimension (:, :), allocatable :: dummy_bit_matrix

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ===========================================================================
        allocate (delta (nhmax, nprowsmax))
        allocate (index1 (nprowsmax))

! Procedure
! ===========================================================================
! Find out which processor this is.
        call MPI_COMM_RANK (MPI_BTN_WORLD, myrank, ierror)

! Sizes of local sections.
        nodiv = norbitals/nactualprocs
        nomod = mod(norbitals,nactualprocs)

! Dummy bit matrix for holding sendrecv return value
        allocate (dummy_bit_matrix (nactualprocs,nactualprocs))

! Determine all the coefficients up to the nchebs'th term.
	if (ipower .eq. 1) then
         call chebft (a, b, c, nchebs, func12)
         nterms = 25
        else if (ipower .eq. -1) then
         call chebft (a, b, c, nchebs, funcm12)
         nterms = 50
        end if

! Scale the range for [a,b] by mapping it into [-1,1]
! So scale S into delta, which fits into new range [-1,1]
	delta = 0.0d0
        do imu = 1, nprows
	 do inu = 1, numh(imu)
          delta(inu,imu) = (2.0d0/(b - a))*s_compact(inu,imu)
	  if (listh(inu,imu) .eq. imu + ipstart - 1)                         & 
     &     delta(inu,imu) = delta(inu,imu) - (b + a)/(b - a)
	 end do
 	end do

        call set_maxdimension (nactualprocs, myrank, nprows, nhmax, nhmax,   &
     &                         nprowsmax, nprows, nprowsmax, numh, listh,    &
     &                         ipstart, numh, listh, ipstart, 0, nS12max_local) 
        call MPI_ALLREDUCE (nS12max_local, nS12max, 1, MPI_INTEGER, MPI_MAX, &
     &                      MPI_BTN_WORLD, ierror)


! Now that we have the parameter nS12max, allocate arrays nums12 and lists12
! and initialize to zero. Also the arrays t0, t1, and t2.
        if (allocated (nums12_local))                                        &
         deallocate (nums12_local, lists12_local, s12_compact_local)
        allocate (nums12_local(nprowsmax))
        allocate (lists12_local(nS12max, nprowsmax))
        allocate (s12_compact_local(nS12max, nprowsmax))
        allocate (t0(nS12max, nprowsmax))
        allocate (t1(nS12max, nprowsmax))
        allocate (t2(nS12max, nprowsmax))
        nums12_local = 0
        lists12_local = 0
        s12_compact_local = 0.0d0
 	t0 = 0.0d0
	t1 = 0.0d0
	t2 = 0.0d0

! Now, calculate the S^2 term. This essentially will yield t2 = 2.0*S*t1 - t0,
! the second term in the Chebyshev polynomial expansion.  
        ichooseA = 1
        ichooseB = 1
! FIXME we should reuse the bit matrix generated when computing C^T from C.
        call sendrecv (myrank, 0, nodiv, nomod, dummy_bit_matrix, .true.,    &
     &                 numh, listh, delta, nprows, nhmax, nprowsmax,         &
     &                 ipstart, numh, listh, delta, nprows, nhmax,           &
     &                 nprowsmax, ipstart, nums12_local, lists12_local,      &
     &                 t2, ichooseA, ichooseB, 1, nS12max, nprowsmax)

! Determine the maximum number of elements in the sparse representation for s12.
! Subtract off the identity matrix. 
        t2 = 2.0d0*t2
 	do imu = 1, nprows
	 do inu = 1, nums12_local(imu)
	  if (lists12_local(inu,imu) .eq. imu + ipstart - 1) then 
           t2(inu,imu) = t2(inu,imu) - 1.0d0 
           s12_compact_local(inu,imu) = 0.5*c(0)
          end if
          do jnu = 1, numh(imu)
           if (lists12_local(inu,imu) .eq. listh(jnu,imu)) then 
            s12_compact_local(inu,imu) =                                     &
     &       s12_compact_local(inu,imu) + c(1)*delta(jnu,imu) 
            t0(inu,imu) = delta(jnu,imu)
           end if
          end do
	 end do
        end do
        s12_compact_local = s12_compact_local + c(2)*t2
        t1 = t2 
        t2 = 0.0d0

! Next terms, i.e. 3 - nterms, in Cheby's expansion
      	do icheb = 3, nterms - 1

! Now, calculate the S^2 term. This essentially will yield t2 = 2.0*S*t1 - t0,
! the second term in the Chebyshev polynomial expansion.  
         ichooseA = 1
         ichooseB = 1
! FIXME we should reuse the bit matrix generated when computing C^T from C.
         call sendrecv (myrank, icheb, nodiv, nomod, dummy_bit_matrix,       &
     &                  .true., numh, listh, delta, nprows, nhmax, nprowsmax,&
     &                  ipstart, nums12_local, lists12_local, t1, nprows,    &
     &                  nS12max, nprowsmax, ipstart, nums12_local,           &
     &                  lists12_local, t2, ichooseA, ichooseB, 1, nS12max,   &
     &                  nprowsmax)
         t2 = 2.0d0*t2 - t0
         s12_compact_local = s12_compact_local + c(icheb)*t2
         t0 = t1
         t1 = t2
         t2 = 0.0d0
        end do


! Deallocate Arrays
! ===========================================================================
        deallocate (dummy_bit_matrix)
        deallocate (delta)
        deallocate (index1)
        deallocate (t0, t1, t2)

! Format Statements
! ===========================================================================
 
        return
        end
