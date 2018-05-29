! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University of Utah - James P. Lewis, Chair
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


! denmata_ordern.f90
! Program Description
! ===========================================================================
!       This routine calculates the density matrices and the band-structure
! energy, as well as similar density matrices.
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
        subroutine denmata_ordern (natoms, nprows, ipstart, ncrows, icstart) 
        use constants_fireball
        use density
        use interactions
        use mpi_declarations
        use neighbor_map
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icstart
        integer, intent (in) :: ipstart
        integer, intent (in) :: natoms
        integer, intent (in) :: ncrows
        integer, intent (in) :: nprows

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irequest, ierror
        integer imu, inu
        integer index
        integer iproc
        integer irecproc
        integer irecrows
        integer irecstart
        integer myrank
        integer nbdiv
        integer nbmod
        integer nFpmax, nFpmax_local
        integer nHijmax, nHijmax_local
        integer nodiv
        integer nomod
        integer npmax, npmax_local
        integer packpos_sparse
        integer packsize_sparse
        integer iprevproc, inextproc
        integer isendpos, irecvpos, itemppos

        integer, dimension (2) :: ichooseA
        integer, dimension (2) :: ichooseB
        integer*1, dimension (:), allocatable :: packarray_sparse
        integer, dimension (MPI_STATUS_SIZE) :: status


! ****************************************************************************
! Temporary storage vectors.
! ****************************************************************************
        integer, dimension (:), allocatable :: numT_local
        integer, dimension (:), allocatable :: numTf_local
        integer, dimension (:, :), allocatable :: listT_local
        integer, dimension (:, :), allocatable :: listTf_local

        real, dimension (:, :), allocatable :: T_compact_local
        real, dimension (:, :), allocatable :: Tf_compact_local
        real, dimension (:, :), allocatable :: Tfs_compact_local

! ****************************************************************************
! F and Fs matrix declaration - local pieces. Defined by F = h*C and Fs = s*C.
! ****************************************************************************
        integer, dimension (:), allocatable, save :: numF_local
        integer, dimension (:, :), allocatable, save :: listF_local

        real, dimension (:, :, :), allocatable, target :: F_Fs_compact_local
        real, dimension (:, :), pointer :: F_compact_local, Fs_compact_local

! ****************************************************************************
! Fp and Fps matrix declaration - local pieces. Defined by F = Hij*C and 
! Fps = Sij*C.
! ****************************************************************************
        integer, dimension (:), allocatable :: numFp_local
        integer, dimension (:, :), allocatable :: listFp_local

        real, dimension (:, :, :), allocatable, target :: Fp_Fps_compact_local
        real, dimension (:, :), pointer :: Fp_compact_local, Fps_compact_local

! ****************************************************************************
! Hij and Sij matrix declaration - local pieces. Defined by Hij = Ct*F and
! Sij = Ct*Fs.
! ****************************************************************************
        integer, dimension (:), allocatable :: numHij_local
        integer, dimension (:, :), allocatable :: listHij_local

        real, dimension (:, :, :), allocatable, target :: Hij_Sij_compact_local
        real, dimension (:, :), pointer :: Hij_compact_local, Sij_compact_local

        logical, dimension (:, :), allocatable :: dummy_bit_matrix

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ============================================================================

! Procedure
! ===========================================================================
! Initialize some things
        call MPI_COMM_RANK (MPI_BTN_WORLD, myrank, ierror)

! sizes of local sections
        nodiv = norbitals/nactualprocs
        nomod = mod(norbitals,nactualprocs)
        nbdiv = nbands/nactualprocs
        nbmod = mod(nbands,nactualprocs)

!        if (myrank .eq. 0) then
!         write (*,*) '  '
!         write (*,*) '  '
!         write (*,*) ' ****************************************************** '
!         write (*,*) '  '
!         write (*,*) '                   Welcome to denmat --              '
!         write (*,*) '        This is the parallel linear-scaling version. '
!         write (*,*) '  '
!         write (*,*) ' ****************************************************** '
!        end if

! Dummy bit matrix for holding sendrecv return value 
        allocate(dummy_bit_matrix(nactualprocs,nactualprocs))

! ****************************************************************************
!
!                      C O M P U T E    D E N S I T I E S
! ****************************************************************************
! Compute Hij and Sij. This section is necessary in order to convert the result
! of the Order N minimization into the desired rho and cape form. This section
! is basically carrying out a transformation.
! Calculate c transpose for the local c matrix contained on this processor.
        call build_transpose (nactualprocs, myrank, dummy_bit_matrix, .true., &
     &                        ncmax, nctmax, nprows,   &
     &                        nprowsmax, ncrows, ncrowsmax, numc_local,      &
     &                        listc_local, c_compact_local, ipstart,         &
     &                        numct_local, listct_local, ct_compact_local,   &
     &                        icstart, norbitals, nbands)


! F/Fs Matrices
! ****************************************************************************
! Allocate arrays. 
        allocate (numF_local (nprowsmax))
        allocate (listF_local (nFmax, nprowsmax))
        allocate (F_Fs_compact_local (nFmax, nprowsmax, 2))
        F_compact_local => F_Fs_compact_local (:, :, 1)
        Fs_compact_local => F_Fs_compact_local (:, :, 2)

! Computes local part of the F and Fs matrices from the coefficient vectors.
!       write (*,*) ' denmat: Building F and Fs matrices, myrank = ', myrank
        ichooseA(1) = 1
        ichooseA(2) = 2
        ichooseB(1) = 1
        ichooseB(2) = 1
        call sendrecv (myrank, 0, nodiv, nomod, hs_bit_matrix, .false., numh, &
     &                 listh, h_s_compact, nprows, nhmax, nprowsmax, ipstart,&
     &                 numc_local, listc_local, c_compact_local, nprows,     &
     &                 ncmax, nprowsmax, ipstart, numF_local, listF_local,   &
     &                 F_Fs_compact_local, ichooseA, ichooseB, 2, nFmax,     &
     &                 nprowsmax)

! H/S Matrices                
! ****************************************************************************
! First Determine the maximum value for the parameter nHijmax.
! Allocate appropriate arrays.
!       write (*,*) ' denmat: Building Hij and Sij matrices, myrank = ', myrank
        call set_maxdimension (nactualprocs, myrank, ncrows, nctmax, nFmax,  &
     &                         ncrowsmax, nprows, nprowsmax, numct_local,    &
     &                         listct_local, icstart, numF_local,            &
     &                         listF_local, ipstart, 0, nHijmax_local)
        call MPI_ALLREDUCE (nHijmax_local, nHijmax, 1, MPI_INTEGER, MPI_MAX, &
     &                      MPI_BTN_WORLD, ierror)

        allocate (numHij_local (ncrowsmax))
        allocate (listHij_local (nHijmax, ncrowsmax))
        allocate (Hij_Sij_compact_local (nHijmax, ncrowsmax, 2))
        Hij_compact_local => Hij_Sij_compact_local (:, :, 1)
        Sij_compact_local => Hij_Sij_compact_local (:, :, 2)

        ichooseA(1) = 1
        ichooseA(2) = 1
        ichooseB(1) = 1
        ichooseB(2) = 2
! FIXME we should reuse the bit matrix generated when computing C^T from C.
        call sendrecv (myrank, 0, nodiv, nomod, dummy_bit_matrix, .true.,      &
     &                 numct_local, listct_local, ct_compact_local, ncrows,  &
     &                 nctmax, ncrowsmax, icstart, numF_local, listF_local,  &
     &                 F_Fs_compact_local, nprows, nFmax, nprowsmax, ipstart,&
     &                 numHij_local, listHij_local, Hij_Sij_compact_local,   &
     &                 ichooseA, ichooseB, 2, nHijmax, ncrowsmax)
        deallocate (numF_local, listF_local, F_Fs_compact_local)

! Redefine Sij = spin*(2.0d0*I - Sij) and Hij = spin*Hij for density matrix 
! purposes.
        do imu = 1, ncrows 
         do inu = 1, numHij_local(imu)
          index = listHij_local(inu,imu)
          if (index .eq. imu + icstart - 1) then 
           Sij_compact_local(inu,imu) =                                      &
     &      spin*(2.0d0 - Sij_compact_local(inu,imu))
          else 
           Sij_compact_local(inu,imu) = spin*Sij_compact_local(inu,imu)
          end if
         end do 
        end do
        Hij_compact_local = spin*Hij_compact_local

! Fp/Fps Matrices
! ****************************************************************************
! Now compute the predecesor to rho and cape - cape_compact_local and 
! rho_compact_local, which is just rho and cape in the compact form for this 
! local processor.
! This time Fp = Hij*c and Fps = Sij*c, and rho_compact_local and 
! cape_compact_local are the final answers.  Note that rho_compact_local and 
! cape_compact_local should be as sparse as the original h_compact. 

! First Determine the maximum value for the parameter nFpmax.
! Allocate appropriate arrays.
!       write (*,*) ' denmat: Building Fp and Fps matrices, myrank = ', myrank
        call set_maxdimension (nactualprocs, myrank, nprows, ncmax, nHijmax, &
     &                         nprowsmax, ncrows, ncrowsmax, numc_local,     &
     &                         listc_local, ipstart, numHij_local,           &
     &                         listHij_local, icstart, 1, nFpmax_local)
        call MPI_ALLREDUCE (nFpmax_local, nFpmax, 1, MPI_INTEGER, MPI_MAX,   &
     &                      MPI_BTN_WORLD, ierror)

        allocate (numFp_local (nprowsmax))
        allocate (listFp_local (nFpmax, nprowsmax))
        allocate (Fp_Fps_compact_local (nFpmax, nprowsmax, 2))
        Fp_compact_local => Fp_Fps_compact_local (:, :, 1)
        Fps_compact_local => Fp_Fps_compact_local (:, :, 2)

        ichooseA(1) = 1
        ichooseA(2) = 1
        ichooseB(1) = 1
        ichooseB(2) = 2
        call sendrecv (myrank, 0, nbdiv, nbmod, dummy_bit_matrix, .true.,      &
     &                 numc_local, listc_local, c_compact_local, nprows,     &
     &                 ncmax, nprowsmax, ipstart, numHij_local,              &
     &                 listHij_local, Hij_Sij_compact_local, ncrows, nHijmax,&
     &                 ncrowsmax, icstart, numFp_local, listFp_local,        &
     &                 Fp_Fps_compact_local, ichooseA, ichooseB, 2, nFpmax,  &
     &                 nprowsmax)
        deallocate (numHij_local, listHij_local, Hij_Sij_compact_local)

! Final Result - rho_compact/cape_compact Matrices         
! ****************************************************************************
! First Determine the maximum value for the parameter npmax.
! Allocate appropriate arrays.
!       write (*,*) ' Building cape and rho matrices, myrank =  ', myrank
        call set_maxdimension (nactualprocs, myrank, nprows, nFpmax, nctmax, &
     &                         nprowsmax, ncrows, ncrowsmax, numFp_local,    &
     &                         listFp_local, ipstart, numct_local,           & 
     &                         listct_local, icstart, 1, npmax_local)
        call MPI_ALLREDUCE (npmax_local, npmax, 1, MPI_INTEGER, MPI_MAX,     &
     &                      MPI_BTN_WORLD, ierror)

        allocate (numcape_local (nprowsmax))
        allocate (listcape_local (npmax, nprowsmax))
        allocate (numrho_local (nprowsmax))
        allocate (listrho_local (npmax, nprowsmax))
        allocate (cape_rho_compact_local (npmax, nprowsmax, 2))
        cape_compact_local => cape_rho_compact_local (:, :, 1)
        rho_compact_local => cape_rho_compact_local (:, :, 2)

        ichooseA(1) = 1
        ichooseA(2) = 2
        ichooseB(1) = 1
        ichooseB(2) = 1
        call sendrecv (myrank, 0, nbdiv, nbmod, dummy_bit_matrix, .true.,      &
     &                 numFp_local, listFp_local, Fp_Fps_compact_local,      &
     &                 nprows, nFpmax, nprowsmax, ipstart, numct_local,      &
     &                 listct_local, ct_compact_local, ncrows, nctmax,       &
     &                 ncrowsmax, icstart, numcape_local, listcape_local,    &
     &                 cape_rho_compact_local, ichooseA, ichooseB, 2, npmax,&
     &                 nprowsmax)

! FIXME is it really necessary to have two separate index lists for 
! cape and rho?
        numrho_local = numcape_local
        listrho_local = listcape_local

        deallocate (numFp_local, listFp_local)
        deallocate (Fp_Fps_compact_local)

! ****************************************************************************
! DENSITY MATRIX -  
! ****************************************************************************
! Collect all of the density matrix pieces from each of the processors and
! store in rho and cape.
! First collect the piece from this processor, then communicate the other
! pieces from the other processors.  Combine to get one rho and cape,
! containing all pieces, on each processor.  For now just do this on the
! master processor, but eventually we may chose to do this on every processor, 
! so that the forces evaluations are distributed over several processors.
!       write (*,*) ' Build sparse density matrix, myrank = ', myrank
        if (myrank .eq. 0) then
         call formrho_sparse (natoms, nprows + ipstart - 1, ipstart) 
        end if

! FIXME:  Do we want to optimize the communications in this part?
! Is the functionality going to change entirely?

! Send local coefficients to the next processor.
        if (nactualprocs .gt. 1) then

! Allocates space for packed sparse matrices.
! Two packets are being sent here - rho_compact and cape_compact, so double
! the packsize_sparse parameter.
         call sparse_getpacksize (npmax, nprowsmax, packsize_sparse)
         packsize_sparse = 2*packsize_sparse
         allocate (packarray_sparse (2*packsize_sparse))
         packpos_sparse = 0
         call sparse_pack (packsize_sparse, nprows, npmax, nprowsmax,        &
     &                     numrho_local, listrho_local, rho_compact_local,   &
     &                     packarray_sparse, packpos_sparse)
         call sparse_pack (packsize_sparse, nprows, npmax, nprowsmax,        &
     &                     numcape_local, listcape_local,                    &
     &                     cape_compact_local, packarray_sparse,             &
     &                     packpos_sparse)

         iprevproc = mod(myrank + nactualprocs - 1, nactualprocs)
         inextproc = mod(myrank + 1, nactualprocs)
         isendpos = 1
         irecvpos = packsize_sparse+1

! Get the density arrays from the other processors.  Upon receiving one set of
! vectors from the previous processor, send it to the next processor.
         do iproc = 1, nactualprocs

! Send packet containing coefficients.
          call MPI_ISEND (packarray_sparse(isendpos), packpos_sparse,        &
     &                    MPI_PACKED, inextproc, 0, MPI_BTN_WORLD, irequest, &
     &                    ierror)
          call MPI_RECV (packarray_sparse(irecvpos), packsize_sparse,        &
     &                   MPI_PACKED, iprevproc, 0, MPI_BTN_WORLD, status,    &
     &                   ierror)
          call MPI_WAIT (irequest, status, ierror)

! From which processor did we receive the vector?
          irecproc = mod(myrank + nactualprocs - iproc, nactualprocs)
          if (irecproc .ne. myrank) then
           irecrows = norbitals/nactualprocs
           if (irecproc .lt. mod(norbitals,nactualprocs)) then
            irecrows = irecrows + 1
            irecstart = irecrows*irecproc + 1
           else
            irecstart = (irecrows + 1)*mod(norbitals,nactualprocs)           &
     &                 + irecrows*(irecproc - mod(norbitals,nactualprocs)) + 1
           end if
 
! Now unpack the packet containing the coeffiecients that was received.
           packpos_sparse = 0
           call sparse_unpack (packarray_sparse(irecvpos), packsize_sparse,  &
     &                         npmax, nprowsmax, numrho_local, listrho_local,&
     &                         rho_compact_local, packpos_sparse)
           call sparse_unpack (packarray_sparse(irecvpos), packsize_sparse,  &
     &                         npmax, nprowsmax, numcape_local,              &
     &                         listcape_local, cape_compact_local,           &
     &                         packpos_sparse)

           if (myrank .eq. 0) then
            call formrho_sparse (natoms, irecrows + irecstart - 1, irecstart) 
           end if
          end if
          itemppos = isendpos
          isendpos = irecvpos
          irecvpos = itemppos
         end do
         deallocate (packarray_sparse)
        end if

! Distribute cape and rho matrices over the processors.
        call MPI_BCAST(rho,numorb_max*numorb_max*natoms*neigh_max, &
             & mpi_whatever_real,0,MPI_BTN_WORLD,ierror)
        call MPI_BCAST(cape,numorb_max*numorb_max*natoms*neigh_max, &
             & mpi_whatever_real,0,MPI_BTN_WORLD,ierror)

! Deallocate Statements
! ===========================================================================
        deallocate (numcape_local, listcape_local)
        deallocate (numrho_local, listrho_local)
        deallocate (cape_rho_compact_local)
        deallocate (dummy_bit_matrix)

! Format Statements
! ===========================================================================
 
        return
        end
