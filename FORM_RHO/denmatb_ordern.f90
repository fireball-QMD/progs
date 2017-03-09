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


! denmatb_ordern.f90
! Program Description
! ===========================================================================
!       This routine calculates the charges from the given density matrices.
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
        subroutine denmatb_ordern (natoms, nprows, ipstart, ncrows, icstart, &
     &                             ifixcharge, iqout)
        use charges
        use constants_fireball
        use density
        use interactions
        use neighbor_map
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icstart
        integer, intent (in) :: ifixcharge
        integer, intent (in) :: ipstart
        integer, intent (in) :: iqout
        integer, intent (in) :: natoms
        integer, intent (in) :: ncrows
        integer, intent (in) :: nprows

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ichooseA
        integer ichooseB
        integer ierror
        integer imu
        integer in1
        integer ineigh
        integer inextproc
        integer inu
        integer iprevproc
        integer iproc
        integer irecproc
        integer irecrows
        integer irecstart
        integer irecvpos
        integer irequest
        integer isendpos
        integer issh
        integer itemppos
        integer jnu
        integer mqn
        integer myrank
        integer nbdiv
        integer nbmod
        integer nodiv
        integer nomod
        integer nnu
        integer packpos_sparse
        integer packsize_sparse

        integer*1, dimension (:), allocatable :: packarray_sparse
        integer, dimension (MPI_STATUS_SIZE) :: istatus

        real aux3

! Temporary storage for Lowdin transformation of the coefficients.
        integer nLmax
        integer nLmax_local
        integer, dimension (:), allocatable :: numL_local 
        integer, dimension (:, :), allocatable :: listL_local 

        real, dimension (:,:), allocatable :: L_compact_local
!$volatile nLmax, nLmax_local, numL_local, listL_local, L_compact_local      

! TESTING
        character (len = 40) append_string
        character (len = 10) extension
        character (len = 40) filename
        character (len = 40) root

        logical, dimension (:, :), allocatable :: dummy_bit_matrix

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Determine which processor this one is.
        call MPI_COMM_RANK (MPI_BTN_WORLD, myrank, ierror)

! Sizes of local sections.
        nbdiv = nbands/nactualprocs
        nbmod = mod(nbands,nactualprocs)
        nodiv = norbitals/nactualprocs
        nomod = mod(norbitals,nactualprocs)

! Dummy bit matrix for holding sendrecv return value
        allocate (dummy_bit_matrix (nactualprocs,nactualprocs))

! ****************************************************************************
!
!  C O M P U T E    L O W D I N    C H A R G E S
! ****************************************************************************
! Initialize to the input charges.
        if (iqout .eq. 1) then
         if (ifixcharge .eq. 1) then
          do iatom = 1, natoms
           QLowdin_TOT(iatom) = 0.0d0
           in1 = imass(iatom)
           do issh = 1, nssh(in1)
            Qout(issh,iatom) = Qin(issh,iatom)
            QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + Qin(issh,iatom) 
           end do
          end do
         else
 
! Calculate the matrix product of ct_compact*s12_compact_local.
          call set_maxdimension (nactualprocs, myrank, nprows, nS12max,      &
     &                           ncmax, nprowsmax, nprows, nprowsmax,        &
     &                           numS12_local, listS12_local, ipstart,       &
     &                           numc_local, listc_local, ipstart, 0,        &
     &                           nLmax_local)
          call MPI_ALLREDUCE (nLmax_local, nLmax, 1, MPI_INTEGER, MPI_MAX,   &
     &                        MPI_BTN_WORLD, ierror)

          ichooseA = 1
          ichooseB = 1
! FIXME we should reuse the bit matrix generated when computing C^T from C.
          allocate (numL_local (nprowsmax))
          allocate (listL_local (nLmax, nprowsmax))
          allocate (L_compact_local (nLmax, nprowsmax))
          call sendrecv (myrank, 0, nodiv, nomod, dummy_bit_matrix, .true.,   &
     &                   numS12_local, listS12_local, S12_compact_local,     &
     &                   nprows, nS12max, nprowsmax, ipstart, numc_local,    &
     &                   listc_local, c_compact_local, nprows, ncmax,        &
     &                   nprowsmax, ipstart, numL_local, listL_local,        &
     &                   L_compact_local, ichooseA, ichooseB, 1, nLmax,      &
     &                   nprowsmax)


! ****************************************************************************
! Lowdin sum 
! ****************************************************************************
! First, determine the Lowdin transformation quantity S^-1/2 - this was done 
! in ss12.f90 called earlier.
! Next, transform the coefficients according to S^-1/2 (Lowdin) transformation.
! Only do this section for the master processor. 

! Collect all of the Lowdin matrix pieces from each of the processors and
! determine the charges for each orbital, Qout.  
! First collect the piece from this processor, then communicate the other
! pieces from the other processors.  Combine to get one Qout, containing all 
! pieces, on each processor.  For now just do this on the master processor, 
! but eventually we may chose to do this on every processor,
! so that the forces evaluations are distributed over several processors.

! Initialize everything to zero.
          do iatom = 1, natoms
           in1 = imass(iatom)
           QLowdin_TOT(iatom) = 0.0d0 
           do issh = 1, nssh(in1)
            Qout(issh,iatom) = 0.0d0
           end do
          end do

          do iatom = 1, natoms
           in1 = imass(iatom)
           do imu = 1, nprows
            aux3 = 0.0d0
            do inu = 1, numL_local(imu)
             aux3 = aux3 + spin*L_compact_local(inu,imu)**2
            end do
            jnu = 0
            do issh = 1, nssh(in1)  
             do mqn = 1, 2*lssh(issh,in1) + 1
              jnu = jnu + 1
              nnu = jnu + degelec(iatom)
              if (imu + ipstart - 1 .eq. nnu) then
               Qout(issh,iatom) = Qout(issh,iatom) + aux3 
               QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux3 
              end if
             end do
            end do
           end do
          end do

! FIXME:  Do we want to optimize the communications in this part?
! Is the functionality going to change entirely?
! Send local coefficients to the next processor.
          if (nactualprocs .gt. 1) then

! Allocates space for packed sparse matrices.
! One packet is being sent here - L_compact
           call sparse_getpacksize (nLmax, nprowsmax, packsize_sparse)
           allocate (packarray_sparse (2*packsize_sparse))
           packpos_sparse = 0
           call sparse_pack_indices (packsize_sparse, nprows, nLmax,         &
     &                               nprowsmax, numL_local, listL_local,     &
     &                               packarray_sparse, packpos_sparse)
           call sparse_pack_elements (packsize_sparse, nprows, nLmax,        &
     &                                nprowsmax, numL_local,                 &
     &                                L_compact_local, packarray_sparse,     &
     &                                packpos_sparse)
           iprevproc = mod(myrank + nactualprocs - 1, nactualprocs)
           inextproc = mod(myrank + 1, nactualprocs)
           isendpos = 1
           irecvpos = packsize_sparse + 1

! Get the Lowdin arrays from the other processors.  Upon receiving one set of
! vectors from the previous processor, send it to the next processor.
           do iproc = 1, nactualprocs

! Send packet containing coefficients.
            call MPI_ISEND (packarray_sparse(isendpos), packpos_sparse,      &
     &                      MPI_PACKED, inextproc, 0, MPI_BTN_WORLD,         &
     &                      irequest, ierror)
            call MPI_RECV (packarray_sparse(irecvpos), packsize_sparse,      &
     &                     MPI_PACKED, iprevproc, 0, MPI_BTN_WORLD, istatus, &
     &                     ierror)
            call MPI_WAIT (irequest, istatus, ierror)

! From which processor did we receive the vector?
            irecproc = mod(myrank + nactualprocs - iproc, nactualprocs)
            if (irecproc .ne. myrank) then
             irecrows = norbitals/nactualprocs
             if (irecproc .lt. mod(norbitals,nactualprocs)) then
              irecrows = irecrows + 1
              irecstart = irecrows*irecproc + 1
             else
              irecstart = (irecrows + 1)*mod(norbitals,nactualprocs)         &
     &                   + irecrows*(irecproc - mod(norbitals,nactualprocs)) + 1
             end if

! Now unpack the packet containing the coeffiecients that was received.
             packpos_sparse = 0
             call sparse_unpack_indices (packarray_sparse(irecvpos),         &
     &                                   packsize_sparse, nLmax, nprowsmax,  & 
     &                                   numL_local, listL_local,            &
     &                                   packpos_sparse)
             call sparse_unpack_elements (packarray_sparse(irecvpos),        &
     &                                    packsize_sparse, nLmax, nprowsmax, &
     &                                    numL_local, L_compact_local,       &
     &                                    packpos_sparse)

! Calculate Lowdin charges.
             do iatom = 1, natoms
              in1 = imass(iatom)
              do imu = 1, nprows
               aux3 = 0.0d0
               do inu = 1, numL_local(imu)
                aux3 = aux3 + spin*L_compact_local(inu,imu)**2
               end do
               jnu = 0
               do issh = 1, nssh(in1)  
                do mqn = 1, 2*lssh(issh,in1) + 1
                 jnu = jnu + 1
                 nnu = jnu + degelec(iatom)
                 if (imu + irecstart - 1 .eq. nnu) then
                  Qout(issh,iatom) = Qout(issh,iatom) + aux3 
                  QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux3 
                 end if
                end do
               end do
              end do
             end do
            end if

! End loop over processors.
            itemppos = isendpos
            isendpos = irecvpos
            irecvpos = itemppos
           end do
           deallocate (packarray_sparse)
          end if
         end if
         deallocate (numL_local, listL_local, L_compact_local)
        end if
 
! Deallocate Arrays
! ===========================================================================
        deallocate (dummy_bit_matrix)
        deallocate (hs_bit_matrix)
        deallocate (numh, listh, h_s_compact)
        deallocate (numct_local, listct_local, ct_compact_local)
 
! Format Statements
! ===========================================================================
400     format (9f9.4)
500     format (70('='))
501     format (2x, ' Atom # ', 2x, ' Type ', 2x, ' Shells ', 1x,' Charges ')
502     format (3x, i5, 7x, a2, 5x, i2, 4x, 8(1x, f5.2))
503     format (3x, i5, 7x, a2, 4x, f10.4)
 
        return
        end
