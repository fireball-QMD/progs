! copyright info:
!
!                             @Copyright 2002
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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

 
! kspace_ordern_slave.f90
! Program Description
! ===========================================================================
!       This routine evaluates the band-structure energy using a linear-
! scaling approach rather than through exact diagonalization which is O(N^3).
! An energy functional is minimized with repect to the coefficients of the
! wavefunctions.  The minimum of this functional yields the band-structure
! energy and a set of orthonormal basis functions.
!
!       This is the slave process routine.
 
! ===========================================================================
! Code written by:
! Spencer Shellman 
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine kspace_ordern_slave (nprocs, my_proc)
        use charges
        use configuration
        use interactions
        use neighbor_map
        use forces
        use mpi_declarations
        use ordern
        use integrals
        use constants_fireball
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatomstart
        integer iblockno
        integer ibsize
        integer pkpos,pksize
        integer icluster
        integer icolumn
        integer icstart
        integer ierror
        integer ifixcharge
        integer ifixneigh
        integer iforce
        integer igauss
        integer imu
        integer index
        integer inu
        integer, save :: ioptionlwf
        integer iordern
        integer ipstart
        integer iqout
        integer irbsize
        integer istepno
        integer iteration
        integer itestrange
        integer itheory
        integer itheory_xc
        integer itry
        integer, save :: itime_step
        integer ivdw
        integer iwrtcharges
        integer iwrtdensity
        integer iwrtneigh
        integer iwrtneigh_com
        integer Kscf
        integer max_scf_iterations
        integer, save :: natoms
        integer natomsp
        integer ncdiv
        integer ncmod
        integer ncrows
        integer nprocs_max
        integer nprows
        integer, save :: nspecies
        integer, save :: nstepi

        integer, dimension (:), allocatable :: ibcastbuffer
        integer*1, dimension (:), allocatable :: ibcb
        integer, dimension (1:24) :: ioff2c           ! for diagnostic purposes
!        integer, dimension (0:14) :: ioff2c         ! for diagnostic purposes
        integer*1, dimension (:), allocatable :: packarray_sparse
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer, dimension (:), allocatable :: displs, recvcounts

        real ebs
        real ebs_local
        real ebs_new
        real eta
        real gnorm
        real gnorm_new
!       real, save :: gnorm_limit
        real testrange
        real sigma
        real sigmaold

        real, dimension (3) :: a1vec, a2vec, a3vec
        integer*1, dimension (:), allocatable :: rbcastbuffer
        real, dimension (:, :), allocatable :: r2_tmp
        real, dimension (:, :, :), allocatable :: r3_tmp
        real, dimension (:, :, :, :), allocatable :: r4_tmp
        real, dimension (:, :, :, :, :), allocatable :: r5_tmp
        real, dimension (nspec_max, nsh_max) :: rcutoff

        character (len = 40) append_string
        character (len = 40) basisfile
        character (len = 10) extension
        character (len = 40) filename
        character (len = 40) root

        logical coefficients
        logical isnan
        logical ordern_stop
        logical ordern_success
        logical scf_achieved
        logical retried

        logical, dimension (:), allocatable :: hs_bit_column

        integer, dimension (:,:), allocatable :: listc_temp
        real, dimension (:,:), allocatable :: c_compact_temp

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize some constants_fireball
        call initconstants (sigma, sigmaold, scf_achieved)

        iordern = 1
        ioptionlwf = 1
        iwrtneigh = 0
        iwrtneigh_com = 0

        ibsize = 19
        allocate (ibcastbuffer(ibsize))

! Broadcast the number of orbitals, etc. to all processors.
        call MPI_BCAST (ibcastbuffer, ibsize, MPI_INTEGER, 0, MPI_COMM_WORLD,&
     &                  ierror)
        numorb_max = ibcastbuffer(1)
        natoms = ibcastbuffer(2)
        nspecies = ibcastbuffer(3)
        norbitals = ibcastbuffer(4)
        nprowsmax = ibcastbuffer(5)
        nstepi = ibcastbuffer(6)
        itheory = ibcastbuffer(7)
        itheory_xc = ibcastbuffer(8)
        icluster = ibcastbuffer(9)
        ifixcharge = ibcastbuffer(10)
        ifixneigh = ibcastbuffer(11)
        iqout = ibcastbuffer(12)
        isorpmax = ibcastbuffer(13)
        ideriv_max = ibcastbuffer(14)
        itestrange = ibcastbuffer(15)
        max_scf_iterations = ibcastbuffer(16)
        igauss = ibcastbuffer(17)
        ivdw = ibcastbuffer(18)
        interactions2c_max = ibcastbuffer(19)
        deallocate (ibcastbuffer)

! The actual number of processors that will work on the matrix.
! For best results the number of actual processors used should be something
! that evenly divides the Hamiltonian matrix into equal sized blocks. 
        do nprocs_max = natoms, 1, -1
         if (mod(natoms,nprocs_max) .eq. 0) exit
        end do
        nactualprocs = nprocs_max
        do nprocs_max = nprocs, 1, -1
         if (mod(norbitals,nprocs_max) .eq. 0) exit
        end do
        if (nactualprocs .gt. nprocs_max) nactualprocs = nprocs_max 

! Create the communication domain consisting of relevant processors.
        call MPI_COMM_SPLIT (MPI_COMM_WORLD, min(my_proc/nactualprocs,1),    &
     &                       my_proc, MPI_BTN_WORLD, ierror)

        if (my_proc .ge. nactualprocs) then
! Do nothing if this processor is irrelevant
         call MPI_FINALIZE ()
999      go to 999
        end if

! Allocate the imass and nelectrons arrays.  Then receive the information 
! contained in these arrays from the master process.
        allocate (degelec (natoms))
        allocate (imass (natoms))
        allocate (lssh (nsh_max, nspecies))
        allocate (nelectron (natoms))
        allocate (nssh (nspecies))
        allocate (num_orb (nspecies))
        allocate (nzx (nspecies))
        allocate (ratom (3, natoms))

        irbsize = 11 + 3*natoms + 3*125 + nspec_max*nsh_max
        call MPI_PACK_SIZE (irbsize, mpi_whatever_real, MPI_BTN_WORLD,       &
     &                      pksize, ierror)
        allocate (rbcastbuffer(pksize))

        call MPI_BCAST (pksize, 1, MPI_INTEGER, 0, MPI_BTN_WORLD, ierror)
        call MPI_BCAST (rbcastbuffer, pksize, MPI_PACKED, 0, MPI_BTN_WORLD,  &
     &                  ierror)

        pkpos = 0
        call MPI_UNPACK (rbcastbuffer, pksize, pkpos, a1vec, 3,              &
     &                   mpi_whatever_real, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (rbcastbuffer, pksize, pkpos, a2vec, 3,              &
     &                   mpi_whatever_real, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (rbcastbuffer, pksize, pkpos, a3vec, 3,              &
     &                   mpi_whatever_real, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (rbcastbuffer, pksize, pkpos, testrange, 1,          &
     &                   mpi_whatever_real, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (rbcastbuffer, pksize, pkpos, range_vdw, 1,          &
     &                   mpi_whatever_real, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (rbcastbuffer, pksize, pkpos, ratom, 3*natoms,       &
     &                   mpi_whatever_real, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (rbcastbuffer, pksize, pkpos, xl, 3*125,             &
     &                   mpi_whatever_real, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (rbcastbuffer, pksize, pkpos, rcutoff,               &
     &                   nspec_max*nsh_max, mpi_whatever_real, MPI_BTN_WORLD,&
     &                   ierror)
        deallocate (rbcastbuffer)

        ibsize = natoms + nsh_max*nspecies + natoms + nspecies + nspecies    &
     &          + natoms + nspecies
        call MPI_PACK_SIZE (ibsize, MPI_INTEGER, MPI_BTN_WORLD, pksize, ierror)
        allocate (ibcb(pksize))

        call MPI_BCAST(pksize,1,MPI_INTEGER,0,MPI_BTN_WORLD,ierror)
        call MPI_BCAST (ibcb, pksize, MPI_PACKED, 0, MPI_BTN_WORLD, ierror)

        pkpos = 0
        call MPI_UNPACK (ibcb, pksize, pkpos, imass, natoms, MPI_INTEGER,    &
     &                   MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (ibcb, pksize, pkpos, nssh, nspecies, MPI_INTEGER,   &
     &                   MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (ibcb, pksize, pkpos, degelec, natoms, MPI_INTEGER,  &
     &                   MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (ibcb, pksize, pkpos, lssh, nsh_max*nspecies,        &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (ibcb, pksize, pkpos, nelectron, natoms, MPI_INTEGER,&
     &                   MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (ibcb, pksize, pkpos, num_orb, nspecies, MPI_INTEGER,&
     &                   MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (ibcb, pksize, pkpos, nzx, nspecies, MPI_INTEGER,    &
     &                   MPI_BTN_WORLD, ierror)
        deallocate (ibcb)

! Calculate neigh_max in parallel
        call allocate_neigh (nactualprocs, my_proc, iordern,         &
     &                             icluster, ivdw, ifixneigh, rcutoff)

! Allocate arrays (since reallocate_f may reallocate them)
        call allocate_f (natoms, neigh_max, numorb_max, itheory, igauss, ivdw)
        call allocate_h (natoms, neigh_max, itheory, igauss)
        call allocate_rho (natoms, neigh_max, numorb_max)

! ME2c_max and ME3c_max must be bcasted before the allocations.
        call ME_max_bcast
        call allocate_ordern (natoms, nspecies)
        call readdata_ordern_init (nspecies, ioff2c)

        istepno = 0
        do
         istepno = istepno + 1
         Kscf = mod(istepno - 1,max_scf_iterations) + 1
         if (itheory .eq. 0) Kscf = 1
! The new value of ratom is needed for neighbors
         call ratom_bcast (natoms, ratom)
! Update neigh_max
         if (Kscf .eq. 1) then
          if (ifixneigh .eq. 0) then
           call reallocate_neigh (nactualprocs, my_proc, iordern,    &
     &                            itheory, igauss, icluster, ivdw, rcutoff)
 
           call neighbors (natoms, nactualprocs, my_proc, iordern, icluster, &
     &                     iwrtneigh, ivdw, basisfile, rcutoff)

          else
           call initneighbors (natoms, ivdw)
          end if
          call common_neighbors (natoms, nactualprocs, my_proc, iordern,     &
     &                           iwrtneigh_com, rcutoff)
         end if

         ibsize = 9
         allocate (ibcastbuffer(ibsize))
         call MPI_BCAST (ibcastbuffer, ibsize, MPI_INTEGER, 0, MPI_BTN_WORLD,&
     &                   ierror)
         itime_step = ibcastbuffer(1)
         nbands = ibcastbuffer(2)
         ncrowsmax = ibcastbuffer(3)
         nhmax = ibcastbuffer(4)
         ncmax = ibcastbuffer(5)
         nctmax = ibcastbuffer(6)
         nFmax = ibcastbuffer(7)
         nFtmax = ibcastbuffer(8)
         iforce = ibcastbuffer(9)
         deallocate (ibcastbuffer)

         if (itheory .eq. 1 .or. itheory .eq. 2) then
          call get_ewald (natoms, nspecies, nactualprocs, my_proc, 0,        &
     &                    icluster, itheory, 1, a1vec, a2vec, a3vec)
         end if

! Assemble two- and three-center interactions in parallel
! Distribute variables for 2 center interactions. 
         iordern = 1
         if (Kscf .eq. 1) then
          call assemble_2c_ordern_init (natoms, nspecies, ivdw)
          call assemble_2c (natoms, nactualprocs, iforce, iordern, ioff2c)
          call assemble_2c_ordern_final (natoms)
         end if

         call Qin_bcast (natoms)
         call Qneutral_bcast (nspecies)

         if (itheory_xc .eq. 0) then
          call assemble_hxc_2c (natoms, nprocs, Kscf, iordern, itheory,      &
     &                          igauss, rcutoff)
         else if (itheory_xc .eq. 1) then
!         call assemble_olsxc_2c (natoms, nprocs, my_proc, Kscf, iordern)
         end if

         if (itheory .eq. 1) then
          call assemble_ca_2c (natoms, nspecies, nactualprocs, iforce,       &
     &                         itheory, iordern, rcutoff)
          call assemble_ca_2c_ordern_final (natoms)
         else if (itheory .eq. 2) then
          call assemble_eh_2c (natoms, nspecies, nactualprocs, iordern)
         end if

! Distribute variables for 3-center interactions.
         if (Kscf .eq. 1)                                                    &
     &    call assemble_3c (natoms, nspecies, nactualprocs, iordern, igauss, &
     &                      rcutoff)

         if (itheory .eq. 1) then
          call assemble_ca_3c (natoms, nspecies, nactualprocs, iordern,      &
     &                         rcutoff)
          call assemble_lr (natoms, nspecies, nactualprocs, iordern)
         end if

         if (igauss .eq. 0 .and. itheory_xc .eq. 0) then
          call assemble_hxc_3c (natoms, nspecies, nactualprocs, Kscf,        &
     &                          iordern, itheory, igauss, rcutoff)
         else if (igauss .eq. 1 .or. itheory_xc .eq. 1) then
!         call assemble_olsxc_3c (natoms, nactualprocs, Kscf, iordern)
         end if

         if (Kscf .eq. 1) call assemble_3c_ordern_final (natoms)
         if (itheory .eq. 1)  call assemble_ca_3c_ordern_final (natoms)

         call buildh (natoms, nactualprocs, itheory, iordern, itestrange,    &
     &                testrange)

! Obtain the Hamiltonian and overlap matrix pieces assigned to this processor.
! Determine the indices of the matrix rows assigned to this processor.
         nprows = norbitals/nactualprocs
         if (my_proc .lt. mod(norbitals,nactualprocs)) then
          nprows = nprows + 1
          ipstart = nprows*my_proc + 1
         else
          ipstart = (nprows + 1)*mod(norbitals,nactualprocs)                 &
                   + nprows*(my_proc - mod(norbitals,nactualprocs)) + 1
         end if

         ncrows = nbands/nactualprocs
         if (my_proc .lt. mod(nbands,nactualprocs)) then
          ncrows = ncrows + 1
          icstart = ncrows*my_proc + 1
         else
          icstart = (ncrows + 1)*mod(nbands,nactualprocs)                    &
                   + ncrows*(my_proc - mod(nbands,nactualprocs)) + 1
         end if

! Allocate the necessary compact arrays for the Hamiltonian and overlap.
! FIXME do we have a memory leak here?
         allocate (numh (nprowsmax))
         allocate (listh (nhmax, nprowsmax))
         allocate (h_s_compact (nhmax, nprowsmax, 2))
         h_compact => h_s_compact (:, :, 1)
         s_compact => h_s_compact (:, :, 2)

! Allocate the bit matrix that will indicate positions of nonzero blocks
! in the h/s matrices, and the temporary matrix for reduction.
         allocate (hs_bit_matrix (nactualprocs, nactualprocs))
         allocate (hs_bit_column (nactualprocs))

         call formsh_compact (natoms, nprows + ipstart - 1, ipstart)

         ncdiv = norbitals/nactualprocs
         ncmod = mod(norbitals,nactualprocs)

! Compute the bit vector for this part of the h/s matrix.
         hs_bit_column = .false.
         do imu = 1, nprows
          do inu = 1, numh(imu)
          index = listh(inu,imu)
           if (index .gt. (ncdiv+1) * ncmod) then
            iblockno = (index - (ncdiv+1) * ncmod - 1) / ncdiv + ncmod + 1
           else
            iblockno = (index - 1) / (ncdiv + 1) + 1
           end if
           hs_bit_column(iblockno) = .true.
          end do
         end do

! Communicate the bit vector over all processors.
! Each entry has a 1 where the corresponding block in the h/s transpose
! is nonzero, or 0 if the block is zero.
! The resulting bit matrix is (theoretically) symmetric.
         call MPI_ALLGATHER (hs_bit_column, nactualprocs, MPI_LOGICAL,       &
     &                       hs_bit_matrix, nactualprocs, MPI_LOGICAL,       &
     &                       MPI_BTN_WORLD, ierror)
         deallocate (hs_bit_column)

! Establish the name of the file to contain the coefficients.
         root = 'coefficients/c_compact'
         write (extension,'(''.'',i3.3)') my_proc
         filename = append_string (root, extension)

! Dimension of the c-vector is number of orbitals by the number of occupied
! bands. Allocate the c-vector and initially randomize the coeffients.
         if (itime_step .eq. nstepi .and. Kscf .eq. 1) then
          allocate (numc_local (nprowsmax))
          allocate (listc_local (ncmax, nprowsmax))
          allocate (c_compact_local (ncmax, nprowsmax))
          inquire (file = filename, exist = coefficients)
!         if (nstepi .eq. 1 .or. .not. coefficients) then 
          if (.not. coefficients) then 
           call formc_compact (natoms, nspecies, ioptionlwf, ratom,          &
     &                         nprows + ipstart - 1, ipstart)
          else
           open (unit = 31, file = filename, status = 'unknown')
           do imu = 1, nprows
            read (31,*) icolumn, numc_local(imu)
            read (31,*) (listc_local(inu,imu), inu = 1, numc_local(imu))
            read (31,*) (c_compact_local(inu,imu), inu = 1, numc_local(imu))
           end do
           close (unit = 31)
          end if
         else
! Reallocate the c matrix if ncmax has changed.
          if (ubound(listc_local,dim=1) .ne. ncmax) then
           allocate(listc_temp(ncmax,nprowsmax),c_compact_temp(ncmax,nprowsmax))
           listc_temp = 0
           c_compact_temp = 0.0d0
           do imu = 1, nprowsmax
            listc_temp(1:numc_local(imu),imu) =                              &
     &       listc_local(1:numc_local(imu),imu)
            c_compact_temp(1:numc_local(imu),imu) =                          &
     &       c_compact_local(1:numc_local(imu),imu)
           end do
           deallocate(listc_local,c_compact_local)
           allocate(listc_local(ncmax,nprowsmax))
           allocate(c_compact_local(ncmax,nprowsmax))
           listc_local = listc_temp
           c_compact_local = c_compact_temp
           deallocate (listc_temp,c_compact_temp)
          end if
         end if

         retried = .false.
         itry = 0

! Come back here if the optimization blows up
! Start iteration of energy and gradients. 
4567     ordern_success = .false.
         iteration = 0
         eta = 0.0d0
         call eandg (nprows, ipstart, ncrows, icstart, iteration, ebs, gnorm,&
     &               0, eta, ordern_stop)
         if (ordern_stop) then
          ordern_success = .true.
         else
          do iteration = 1, max_ordern_iterations
           call eandg (nprows, ipstart, ncrows, icstart, iteration, ebs_new, &
     &                 gnorm_new, 0, eta, ordern_stop)
           if (ordern_stop) then
            ordern_success = .true.
            ebs = ebs_new
            gnorm = gnorm_new
            exit
           end if

! Test the relative error.
           if (abs(ebs - ebs_new) .le. ordern_tolerance*abs(ebs_new)         &
     &         .and. (gnorm_new .le. ordern_grad_tolerance                   &
     &                .or. abs(gnorm - gnorm_new) .le. ordern_tolerance)) then
            ordern_success = .true.
            ebs = ebs_new
            gnorm = gnorm_new
            exit
           end if
           ebs = ebs_new
           gnorm = gnorm_new
          end do
         end if
         call eandg_cleanup

! Did the optimization blow up?
         if (abs(ebs) .gt. ebs_blowup .or. ebs .gt. 0.0d0 .or. isnan(ebs)) then
! Do not retry more than once.
          if (retried) stop
          call formc_compact (natoms, nspecies, ioptionlwf, ratom,           &
     &                        nprows + ipstart - 1, ipstart)
          itry = itry + 1
          if (itry .eq. 10) retried = .true.
          go to 4567
         end if

! save gradient norm to use as threshold at next step
!        gnorm_limit = max(gnorm, ordern_grad_tolerance)

! Write out the coefficients into a storage file - for restart capabilities.
         open (unit = 31, file = filename, status = 'unknown')
         do imu = 1, nprows
          write (31,*) imu + ipstart - 1, numc_local(imu)
          write (31,*) (listc_local(inu,imu), inu = 1, numc_local(imu))
          write (31,*) (c_compact_local(inu,imu), inu = 1, numc_local(imu))
         end do
         close (unit = 31)

         call denmata_ordern (natoms, nprows, ipstart, ncrows, icstart) 
         if (ifixcharge .ne. 1 .and. iqout .ne. 2 .and. Kscf .eq. 1)         &
     &    call ss12 (nprows, ipstart, 1)
         call denmatb_ordern (natoms, nprows, ipstart, ncrows, icstart,      &
     &                        ifixcharge, iqout)
! learn from processor 0 if scf is achieved
         if (itheory .eq. 0) then
          scf_achieved = .true.
         else
          call MPI_BCAST (scf_achieved, 1, MPI_LOGICAL, 0, MPI_BTN_WORLD,    &
     &                    ierror)
         end if
         if (scf_achieved) then
! exit the scf loop
          istepno = 0
          Kscf = max_scf_iterations
         end if
         if (itheory .eq. 0 .or. Kscf .eq. max_scf_iterations) then
! end of the scf loop
          if (itheory .eq. 1 .or. itheory .eq. 2) then 
           call get_ewald (natoms, nspecies, nactualprocs, my_proc, iforce,  &
     &                     icluster, itheory, 1, a1vec, a2vec, a3vec)

          end if
          if (iforce .eq. 1) then
           call Dassemble_2c (natoms, nactualprocs, iordern, igauss)

           call Qin_bcast (natoms)
           call Qneutral_bcast (nspecies)

           if (itheory_xc .eq. 0) then
            call Dassemble_hxc_2c (natoms, nprocs, iordern, itheory, igauss, &
     &                             rcutoff)
           else if (itheory_xc .eq. 1) then
!           call Dassemble_olsxc_2c (natoms, nprocs, my_proc, iordern, jspin)
           end if
           call Dassemble_2c_ordern_final (natoms)

           if (itheory .eq. 1) then
            call Dassemble_ca_2c (natoms, nspecies, nactualprocs, iordern,       &
     &                            rcutoff)
            call Dassemble_ca_2c_ordern_final (natoms)
           else if (itheory .eq. 2) then
            call Dassemble_eh_2c (natoms, nspecies, nactualprocs, my_proc,   &
     &                            iordern)
           end if

           call Dassemble_3c (natoms, nspecies, nactualprocs, iordern,       &
     &                        igauss, rcutoff)

           if (itheory .eq. 1) then
            call Dassemble_ca_3c (natoms, nspecies, nactualprocs, iordern,   &
     &                            rcutoff)
            call Dassemble_lr (natoms, nspecies, nactualprocs, iordern)
           end if

           if (igauss .eq. 0 .and. itheory_xc .eq. 0) then
            call Dassemble_hxc_3c (natoms, nspecies, nactualprocs, iordern,  &
     &                             itheory, igauss, rcutoff)
           else if (igauss .eq. 1 .or. itheory_xc .eq. 1) then
!           call Dassemble_olsxc_3c (natoms, nactualprocs, Kscf, iordern)
           end if

           call Dassemble_3c_ordern_final (natoms)
           if (itheory .eq. 1) call Dassemble_ca_3c_ordern_final (natoms)
          end if
         end if
        end do

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
        return
        end subroutine kspace_ordern_slave
