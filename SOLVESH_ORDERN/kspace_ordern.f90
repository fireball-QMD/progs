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
 
! kspace_ordern.f90
! Program Description
! ===========================================================================
!       This routine evaluates the band-structure energy using a linear-
! scaling approach rather than through exact diagonalization which is O(N^3).
! An energy functional is minimized with repect to the coefficients of the
! wavefunctions.  The minimum of this functional yields the band-structure
! energy and a set of orthonormal basis functions.
 
! ===========================================================================
! Original Order-N compilation written by Spencer Shellman.

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
        subroutine kspace_ordern (natoms, nspecies, nprocs, my_proc, Kscf,   &
     &                            itime_step, nstepi, ratom, nprows,   &
     &                            ipstart, ncrows, icstart, ebs) 
        use charges
        use interactions
        use neighbor_map
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step
        integer, intent (in) :: Kscf
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies
        integer, intent (in) :: nstepi
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc
 

        real, intent (in), dimension (3, natoms) :: ratom

! Output
        integer, intent (out) :: icstart
        integer, intent (out) :: ipstart
        integer, intent (out) :: ncrows
        integer, intent (out) :: nprows
        real, intent (out) :: ebs

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iblockno
        integer icolumn
        integer ierror
        integer imu
        integer index
        integer inu
        integer ioptionlwf
        integer iproc
        integer isendrows
        integer isendstart
        integer iteration
        integer itry
        integer ncdiv 
        integer ncmod
        integer nprocs_max
 
        integer*1, dimension (:), allocatable :: packarray_sparse

        real ebs_local
        real ebs_new
        real gnorm
        real gnorm_new
!       real, save :: gnorm_limit
        real eta

        character (len = 40) append_string
        character (len = 10) extension
        character (len = 40) filename
        character (len = 40) root

        logical coefficients
        logical isnan
        logical ordern_success
        logical ordern_stop
        logical retried

        logical, dimension (:), allocatable :: hs_bit_column

        integer, dimension (:,:), allocatable :: listc_temp
        real, dimension (:,:), allocatable :: c_compact_temp

! date/time variables
        integer, dimension (8) :: iv
        character (len = 13) cd, ct, cz

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
        write (*,*) '  '
        write (*,*) ' ****************************************************** '
        write (*,*) '  '
        write (*,*) '                   Welcome to kspace --              '
        write (*,*) '        This is the parallel linear-scaling version. '
        write (*,*) '  '
        write (*,*) ' ****************************************************** '
        ioptionlwf = 1

! BUILD COMPACT S AND H OVER NUMBER OF PROCESSORS.
! ****************************************************************************
! We build H and S going through the list of atoms and their neighbors.
! How do we know where the particular (i,j) spot belongs in the grand
! matrix?  This should help you understand the degelec shifting business:
!
!                  atom 1         atom2            atom  3
!
!              nu1 nu2 nu3   nu1 nu2 nu3 nu4     nu1 nu2 nu3
!
!                           _________________
!         mu1               |               |
!  atom 1 mu2               |    H(1,2)     |
!         mu3               |               |
!                           -----------------
!         mu1
!  atom 2 mu2
!         mu3
!         mu4
!
!         mu1
!  atom3  mu2
!         mu3
!
! to the (1,2) portion at the right place we use degelec(iatom), which is
! passed, it remembers how many orbitals one must skip to get to the spot
! reserved for iatom, e.g. in this example degelec(1)=0, degelc(2)=3.

! Determine the number of the matrix rows assigned to this processor.
        nprows = norbitals/nactualprocs
        if (mod(norbitals,nactualprocs) .gt. 0) nprows = nprows + 1

        ncrows = nbands/nactualprocs
        if (mod(nbands,nactualprocs) .gt. 0) ncrows = ncrows + 1

! Allocate the necessary compact arrays for the Hamiltonian and overlap.
        allocate (numh (nprowsmax))
        allocate (listh (nhmax, nprowsmax))
        allocate (h_s_compact (nhmax, nprowsmax, 2))
        h_compact => h_s_compact (:, :, 1)
        s_compact => h_s_compact (:, :, 2)
 
! Allocate the bit matrix that will indicate positions of nonzero blocks
! in the h/s matrices, and the temporary matrix for reduction.
        allocate (hs_bit_matrix (nactualprocs, nactualprocs))
        allocate (hs_bit_column (nactualprocs))

! For each section of the Hamiltonian and overlap matrices, computes the sparse
! representation of the section and sends it to a processor.  Note that the 
! segments are transposed; this requires less storage and facilitates sparse 
! matrix multiplication (see eandg).
        write (*,*) '  '
        write (*,*) ' Building H and S; pieces built on different processors. ' 

! Determine the indices of the matrix rows assigned to this processor.
        ipstart = 1
        icstart = 1

! Formulate h_compact and s_compact for the 0th processor
        call formsh_compact (natoms, nprows + ipstart - 1, ipstart)

        ncdiv = norbitals/nactualprocs
        ncmod = mod(norbitals,nactualprocs)

! Compute the bit vector for this section of h/s.
        hs_bit_column = .false.
        do imu = 1, nprows
         do inu = 1, numh(imu)
          index = listh(inu,imu)
          if (index .gt. (ncdiv+1) * ncmod) then
           iblockno = (index - (ncdiv + 1)*ncmod - 1)/ncdiv + ncmod + 1
          else
           iblockno = (index - 1)/(ncdiv + 1) + 1
          end if
          hs_bit_column (iblockno) = .true.
         end do
        end do

! Communicate the bit vector over all processors.
! Each entry has a 1 where the corresponding block in h/s is nonzero, or 0 if 
! the block is zero. The resulting bit matrix is (theoretically) symmetric.
        call MPI_ALLGATHER (hs_bit_column, nactualprocs, MPI_LOGICAL,        &
     &                      hs_bit_matrix, nactualprocs, MPI_LOGICAL,        &
     &                      MPI_BTN_WORLD, ierror)
        deallocate (hs_bit_column)

!
! "DIAGONALIZE" THE HAMILTONIAN.
! ****************************************************************************
! The NEW way of doing things.
! The "diagonalization" is performed via the methods of conjugate gradients
! in the following manner:  An energy functional is formulated with
! corresponding gradients of the functional.  This functional is then
! minimized with respect to the calculated gradients which are taken with
! respect to the cooefficients.

! Establish the name of the file to contain the coefficients.
        root = 'coefficients/c_compact'
        write (extension,'(''.'',i3.3)') my_proc
        filename = append_string (root, extension)

! First allocate the c-vector and initially randomize the coefficients.
! Dimension of the c-vector is number of orbitals by the number of occupied
! bands. 
        if (itime_step .eq. nstepi .and. Kscf .eq. 1) then
         allocate (numc_local (nprowsmax))
         allocate (listc_local (ncmax, nprowsmax))
         allocate (c_compact_local (ncmax, nprowsmax))
         inquire (file = filename, exist = coefficients)
!        if (nstepi .eq. 1 .or. .not. coefficients) then
         if (.not. coefficients) then
          write (*,*) ' Building random wavefunction coefficients (c''s) '
          call formc_compact (natoms, nspecies, ioptionlwf, ratom,           &
     &                        nprows + ipstart - 1, ipstart)
         else
          write (*,*) ' Reading wavefunction coefficients (c''s) from files. '
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
                 listc_temp(1:numc_local(imu),imu) = listc_local(1:numc_local(imu),imu)
                 c_compact_temp(1:numc_local(imu),imu) = c_compact_local(1:numc_local(imu),imu)
              end do
              deallocate(listc_local,c_compact_local)
              allocate(listc_local(ncmax,nprowsmax),c_compact_local(ncmax,nprowsmax))
              listc_local = listc_temp
              c_compact_local = c_compact_temp
              deallocate(listc_temp,c_compact_temp)
           end if
        end if

        retried = .false.
        itry = 0
! come back here if the optimization blows up

! Start iteration of energy and gradients.
4567    write (*,100) 
        write (*,*) ' Call order-N minimization routine for Hamiltonian. '
        write (*,100) 
 
! Reset gradient norm threshold, otherwise use the final gradient norm
! from the previous step
!       gnorm_limit = ordern_grad_tolerance

! Show time at start of run.
!       call date_and_time (cd, ct, cz, iv)
!       write (*,*) 'Processor 0: time at start of order-N = ', ct
        ordern_success = .false.
        iteration = 0
        eta = 0.0d0
        call eandg (nprows, ipstart, ncrows, icstart, iteration, ebs, gnorm, &
     &              0, eta, ordern_stop)
        write (*,200) iteration, ebs, gnorm
        if (ordern_stop) then
         ordern_success = .true.
        else
         do iteration = 1, max_ordern_iterations
          call eandg (nprows, ipstart, ncrows, icstart, iteration, ebs_new,  &
     &                gnorm_new, 0, eta, ordern_stop)
          write (*,200) iteration, ebs_new, gnorm_new

! Show time at end of iteration.
!         call date_and_time (cd, ct, cz, iv)
!         write (*,*) 'Processor 0: time at end of iteration = ', ct
          if (ordern_stop) then
           ordern_success = .true.
           ebs = ebs_new
           gnorm = gnorm_new
           exit
          end if

! Write out the coefficients into a storage file - for restart capabilities.
!         if (.not. (abs(ebs) .gt. ebs_blowup .or. ebs .gt. 0.0d0)) then 
!          open (unit = 31, file = filename, status = 'unknown')
!          do imu = 1, nprows
!           write (31,*) imu + ipstart - 1, numc_local(imu)
!           write (31,*) (listc_local(inu,imu), inu = 1, numc_local(imu))
!           write (31,*) (c_compact_local(inu,imu), inu = 1, numc_local(imu))
!          end do
!          close(unit = 31)
!         end if

! JPL 2000 FIXME
! Find the RMS of the gradient, this should also fall under some tolerance
! before ending the iteration cycle.
! Tests the relative error.
          if (abs (ebs - ebs_new) .le. ordern_tolerance*abs(ebs_new)         &
     &        .and. (gnorm_new .le. ordern_grad_tolerance                    &
     &               .or. abs(gnorm - gnorm_new) .le. ordern_tolerance)) then
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
        if (.not. ordern_success) then
         write (*,*) ' ERROR: order-N routine made ', max_ordern_iterations, &
     &               ' iterations, not successful'
        end if
        write (*,100) 

! Did the optimization blow up?
        if (abs(ebs) .gt. ebs_blowup .or. ebs .gt. 0.0d0 .or. isnan(ebs)) then
         if (retried) then
! Do not retry more than once. 
          write (*,*) ' FATAL ERROR: Retried optimization step and it blew '
          write (*,*) ' up again.  Exiting...'
          stop
         end if
         write (*,*) ' ERROR: Optimization blew up! '
         write (*,*) ' Randomizing coefficients and trying again... '
         call formc_compact (natoms, nspecies, ioptionlwf, ratom,            &
     &                       nprows + ipstart - 1, ipstart)
         itry = itry + 1
         if (itry .eq. 10) retried = .true.
         go to 4567
        end if

! Save gradient norm to use as the threshold at the next step,
! to avoid excessive iterations
!       gnorm_limit = max(gnorm, ordern_grad_tolerance)

! Write out the coefficients into a storage file - for restart capabilities.
        open (unit = 31, file = filename, status = 'unknown')
        do imu = 1, nprows
         write (31,*) imu + ipstart - 1, numc_local(imu)
         write (31,*) (listc_local(inu,imu), inu = 1, numc_local(imu))
         write (31,*) (c_compact_local(inu,imu), inu = 1, numc_local(imu))
        end do
        close (unit = 31)
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (2x, ' Processor0: Iteration ', i3, 1x, ' energy = ', f15.4, &
                1x, ' gnorm = ', f11.6) 
301     format (16f8.3)
302     format (24f8.3)
 
        return
        end



        subroutine readdata_ordern_init (nspecies, ioff2c)
        use dimensions
        use forces
        use interactions
        use integrals
        use charges
        use constants_fireball
        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nspecies

! 17 different two-center interactions
        integer, intent (inout), dimension (1:24) :: ioff2c
!        integer, intent (inout), dimension (0:14) :: ioff2c
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror, krbsize

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

!         call bcast(ioff2c,15,MPI_INTEGER)
        call bcast(ioff2c,17,MPI_INTEGER)
        call bcast(icon3c,nspecies**3,MPI_INTEGER)
!        call bcast(ind2c,15*9,MPI_INTEGER)
!        call bcast(ind2c,17*9,MPI_INTEGER)
        call bcast(ind2c,20*9,MPI_INTEGER)
        call bcast(numx3c_bcna,(isorpmax+1)*(nspecies**3),MPI_INTEGER)
        call bcast(numy3c_bcna,(isorpmax+1)*(nspecies**3),MPI_INTEGER)
        call bcast(numz2c,interactions2c_max*(nspecies**2),MPI_INTEGER)
        call bcast(numx3c_xc3c,(ideriv_max+1)*(nspecies**3),MPI_INTEGER)
        call bcast(numy3c_xc3c,(ideriv_max+1)*(nspecies**3),MPI_INTEGER)
        call bcast(index_max2c,nspecies**2,MPI_INTEGER)
        call bcast(index_max3c,nspecies**2,MPI_INTEGER)
        call bcast(mu,ME3c_max*(nspecies**2),MPI_INTEGER)
        call bcast(mvalue,ME3c_max*(nspecies**2),MPI_INTEGER)
        call bcast(nu,ME3c_max*(nspecies**2),MPI_INTEGER)
        call bcast(index_maxPP,nspecies**2,MPI_INTEGER)
        call bcast(muPP,ME3c_max*(nspecies**2),MPI_INTEGER)
        call bcast(nuPP,ME3c_max*(nspecies**2),MPI_INTEGER)
        call bcast(lsshPP,nsh_max*nspecies,MPI_INTEGER)
        call bcast(num_orbPP,nspecies,MPI_INTEGER)
        call bcast(nsshPP,nspecies,MPI_INTEGER)
        call bcast(mu2shell,numorb_max*nspecies,mpi_integer)
        call bcast(dq,nspecies,mpi_whatever_real)
        call bcast(hx_bcna,(isorpmax+1)*(nspecies**3),mpi_whatever_real)
        call bcast(hy_bcna,(isorpmax+1)*(nspecies**3),mpi_whatever_real)
        call bcast(x3cmax_bcna,(isorpmax+1)*(nspecies**3),mpi_whatever_real)
        call bcast(y3cmax_bcna,(isorpmax+1)*(nspecies**3),mpi_whatever_real)
        if (.not. superspline) &
             & call bcast(xintegral_2c,ME2c_max*nfofx*interactions2c_max*(nspecies**2),mpi_whatever_real)
        call bcast(z2cmax,interactions2c_max*(nspecies**2),mpi_whatever_real)
        if (superspline) &
             & call bcast(splineint_2c,4*ME2c_max*nfofx*interactions2c_max*(nspecies**2),mpi_whatever_real)
        call bcast(hx_xc3c,(ideriv_max+1)*(nspecies**3),mpi_whatever_real)
        call bcast(hy_xc3c,(ideriv_max+1)*(nspecies**3),mpi_whatever_real)
        call bcast(x3cmax_xc3c,(ideriv_max+1)*(nspecies**3),mpi_whatever_real)
        call bcast(y3cmax_xc3c,(ideriv_max+1)*(nspecies**3),mpi_whatever_real)
        call bcast(cl_PP,nsh_max*nspecies,mpi_whatever_real)

        krbsize = ME3c_max*numXmax*numYmax*(isorpmax+1)*(nspecies**3)

        call bcast (bcna_01, krbsize, mpi_whatever_real)
        call bcast (bcna_02, krbsize, mpi_whatever_real)
        call bcast (bcna_03, krbsize, mpi_whatever_real)
        call bcast (bcna_04, krbsize, mpi_whatever_real)
        call bcast (bcna_05, krbsize, mpi_whatever_real)

        krbsize = ME3c_max*numXmax*numYmax*(ideriv_max+1)*(nspecies**3)

        call bcast (xc3c_01, krbsize, mpi_whatever_real)
        call bcast (xc3c_02, krbsize, mpi_whatever_real)
        call bcast (xc3c_03, krbsize, mpi_whatever_real)
        call bcast (xc3c_04, krbsize, mpi_whatever_real)
        call bcast (xc3c_05, krbsize, mpi_whatever_real)

        end subroutine readdata_ordern_init


        subroutine assemble_ordern_sub_ewald (natoms, sub_ewald)
        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        real, dimension (natoms), intent (inout) :: sub_ewald
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror

        real, dimension (:), allocatable :: r1_tmp

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ========================================================================

        allocate (r1_tmp (natoms))

        call MPI_ALLREDUCE (sub_ewald, r1_tmp, natoms, &
             & mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        sub_ewald = r1_tmp

        deallocate (r1_tmp)
 
        return
        end




        subroutine assemble_ordern_sub_dewald (natoms, sub_ewald, sub_dewald)
        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        real, dimension (natoms), intent (inout) :: sub_ewald
        real, dimension (3, natoms), intent (inout) :: sub_dewald
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror

        real, dimension (:), allocatable :: r1_tmp
        real, dimension (:,:), allocatable :: r2_tmp

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ========================================================================

        allocate (r1_tmp (natoms), r2_tmp(3,natoms))

        call MPI_ALLREDUCE (sub_ewald, r1_tmp, natoms, &
             & mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        sub_ewald = r1_tmp
        call MPI_ALLREDUCE (sub_dewald, r2_tmp, 3*natoms, &
             & mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        sub_dewald = r2_tmp

        deallocate (r1_tmp,r2_tmp)
 
        return
        end



        subroutine Dassemble_eh_2c_ordern_final (natoms)
        use mpi_declarations
        use ordern
        use charges
        use forces
        use neighbor_map
        use interactions
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror

        real, dimension (:,:), allocatable :: r2_tmp

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ========================================================================

        allocate (r2_tmp (3, natoms))

        call MPI_ALLREDUCE (fcoulomb, r2_tmp, 3*natoms, &
             & mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        fcoulomb = r2_tmp
        call MPI_ALLREDUCE (fxcnu, r2_tmp, 3*natoms, &
             & mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        fxcnu = r2_tmp
        call MPI_ALLREDUCE (flrew, r2_tmp, 3*natoms, &
             & mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        flrew = r2_tmp

        deallocate (r2_tmp)
 
        return
        end



        subroutine ewald_energy_ordern_final (natoms, iforce, icluster, itheory)
        use charges
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: iforce
        integer, intent (in) :: itheory
        integer, intent (in) :: natoms
 
! Local Parameters and Data Declaration
! ===========================================================================
        real, dimension (:, :), allocatable :: r2_tmp
        real, dimension (:, :, :), allocatable :: r3_tmp
        integer ierror

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================

! sum the ewald matrices over all processors
        allocate (r2_tmp (natoms, natoms))

        call MPI_ALLREDUCE (ewald, r2_tmp, natoms**2, mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)

        ewald = r2_tmp
        deallocate (r2_tmp)
        if (iforce .eq. 1) then
           allocate (r2_tmp (3, natoms))

           call MPI_ALLREDUCE (fewald, r2_tmp, 3*natoms, mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
           fewald = r2_tmp
           deallocate (r2_tmp)
           allocate (r3_tmp (3, natoms, natoms))
           call MPI_ALLREDUCE (dewald, r3_tmp, 3*(natoms**2), mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)

           dewald = r3_tmp
           deallocate (r3_tmp)
        end if

        return
        end


        subroutine get_allgatherv_arrays (nelements,nprocs,my_proc,c,nelementsp,displs,recvcounts)
          implicit none

          integer, intent (in) :: nelements,nprocs,my_proc,c
          integer, intent (out) :: nelementsp
          integer, dimension (0:nprocs-1), intent (out) :: displs,recvcounts

          integer imu

          do imu = 0, nprocs-1
             nelementsp = nelements/nprocs
             if (imu .lt. mod(nelements,nprocs)) then
                nelementsp = nelementsp + 1
                displs(imu) = nelementsp*imu
             else
                displs(imu) = (nelementsp + 1)*mod(nelements,nprocs)                   &
                     + nelementsp*(imu - mod(nelements,nprocs))
             end if
             recvcounts(imu) = nelementsp
          end do
          nelementsp = nelements/nprocs
          if (my_proc .lt. mod(nelements,nprocs)) nelementsp = nelementsp + 1

          displs = displs * c
          recvcounts = recvcounts * c
          nelementsp = nelementsp * c

          return
        end subroutine get_allgatherv_arrays


        subroutine allgatherv_i (elements,my_proc,total,subtotal,displs,recvcounts)
          use ordern
          implicit none

          include 'mpif.h'

          integer, intent (in) :: my_proc,total,subtotal
          integer, dimension (total), intent (inout) :: elements
          integer, dimension (*), intent (in) :: displs,recvcounts
!$ volatile elements
          integer, dimension (:), allocatable :: temp
          integer ierror

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

          allocate (temp(total))
          call MPI_ALLGATHERV(elements(displs(my_proc+1)+1),subtotal,MPI_INTEGER, &
               & temp,recvcounts,displs,MPI_INTEGER,MPI_BTN_WORLD,ierror)

          elements = temp
          deallocate(temp)
          return
        end subroutine allgatherv_i


        subroutine allgatherv_r (elements,my_proc,total,subtotal,displs,recvcounts)
          use mpi_declarations
          use ordern
          implicit none

          include 'mpif.h'

          integer, intent (in) :: my_proc,total,subtotal
          real, dimension (total), intent (inout) :: elements
          integer, dimension (*), intent (in) :: displs,recvcounts
!$ volatile elements
          real, dimension (:), allocatable :: temp
          integer ierror

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

          allocate(temp(total))
          call MPI_ALLGATHERV(elements(displs(my_proc+1)+1),subtotal,mpi_whatever_real, &
               & temp,recvcounts,displs,mpi_whatever_real,MPI_BTN_WORLD,ierror)

          elements = temp
          deallocate(temp)
          return
        end subroutine allgatherv_r


        subroutine buildh_ordern_final (natoms,nprocs,my_proc,itheory)
        use charges
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use ordern
        use neighbor_map
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms,nprocs,my_proc,itheory
  
! Local Parameters and Data Declaration
! ===========================================================================
        integer, dimension (0:nprocs-1) :: displs, recvcounts
        integer c
        integer natomsp
        integer ierror

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================

! gather the buildh matrices over all processors
           c = numorb_max*numorb_max*neigh_max
           call get_allgatherv_arrays (natoms,nprocs,my_proc,c,natomsp,displs,recvcounts)
           call allgatherv_r (h_mat,my_proc,natoms*c,natomsp,displs,recvcounts)
           call allgatherv_r (s_mat,my_proc,natoms*c,natomsp,displs,recvcounts)
           call allgatherv_r (t_mat,my_proc,natoms*c,natomsp,displs,recvcounts)
           call allgatherv_r (vna,my_proc,natoms*c,natomsp,displs,recvcounts)
           call allgatherv_r (vnl,my_proc,natoms*c,natomsp,displs,recvcounts)
           call allgatherv_r (vxc,my_proc,natoms*c,natomsp,displs,recvcounts)
           if (itheory.eq.1 .or. itheory.eq.2) then
              call allgatherv_r (ewaldlr,my_proc,natoms*c,natomsp,displs,recvcounts)
              call allgatherv_r (ewaldsr,my_proc,natoms*c,natomsp,displs,recvcounts)
           end if

        return
        end


        subroutine find_neigh_max_ordern_final ()
        use charges
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use ordern
        use neighbor_map
        implicit none

        include 'mpif.h'

        integer neigh_max_temp, ierror
 
! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Maximize neigh_max over all processors
        call MPI_ALLREDUCE (neigh_max, neigh_max_temp, 1, MPI_INTEGER,       &
     &                      MPI_MAX, MPI_BTN_WORLD, ierror)
        neigh_max = neigh_max_temp

        call MPI_ALLREDUCE (neigh_max_vdw, neigh_max_temp, 1, MPI_INTEGER,   &
     &                      MPI_MAX, MPI_BTN_WORLD, ierror)
        neigh_max_vdw = neigh_max_temp

        return
        end



        subroutine neighbors_ordern_final (natoms, nprocs, my_proc, ivdw)
        use charges
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use ordern
        use neighbor_map
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: my_proc
        integer, intent (in) :: natoms
        integer, intent (in) :: nprocs, ivdw
  
! Local Parameters and Data Declaration
! ===========================================================================
        integer, dimension (0:nprocs - 1) :: displs, recvcounts
        integer c
        integer natomsp

! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
! Gather the neighj/neighb/neighn arrays over all processors.
! Gather the neighj_vdw/neighb_vdw/neighn_vdw arrays over all processors.
        call get_allgatherv_arrays (natoms, nprocs, my_proc, 1, natomsp,     &
     &                              displs, recvcounts)
        call allgatherv_i (neighn, my_proc, natoms, natomsp, displs, recvcounts)
        if (ivdw.eq.1) call allgatherv_i (neighn_vdw, my_proc, natoms, natomsp, displs,     &
     &                     recvcounts)
        c = neigh_max
        call allgatherv_i (neigh_b, my_proc, natoms*c, natomsp*c, displs*c,      &
     &                     recvcounts*c)
        call allgatherv_i (neigh_j, my_proc, natoms*c, natomsp*c, displs*c,      &
     &                     recvcounts*c)
        if (ivdw.eq.1) then
           c = neigh_max_vdw
           call allgatherv_i (neigh_b_vdw, my_proc, natoms*c, natomsp*c, displs*c,  &
                &                     recvcounts)
           call allgatherv_i (neigh_j_vdw, my_proc, natoms*c, natomsp*c, displs*c,  &
                &                     recvcounts*c)
        end if
        return
        end


        subroutine common_neighbors_ordern_final (natoms,nprocs,my_proc)
        use charges
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use ordern
        use neighbor_map
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms,nprocs,my_proc
  
! Local Parameters and Data Declaration
! ===========================================================================
        integer, dimension (0:nprocs-1) :: displs, recvcounts
        integer c
        integer natomsp

! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ==========================================================================

! gather the neigh_com arrays over all processors
        call get_allgatherv_arrays (natoms,nprocs,my_proc,1,natomsp,displs,recvcounts)
        call allgatherv_i (neigh_comn,my_proc,natoms,natomsp,displs,recvcounts)
        c = neigh_max**2
        natomsp = natomsp * c
        displs = displs * c
        recvcounts = recvcounts * c
        call allgatherv_i (neigh_comm,my_proc,natoms*c,natomsp,displs,recvcounts)
        c = c * 2
        natomsp = natomsp * 2
        displs = displs * 2
        recvcounts = recvcounts * 2
        call allgatherv_i (neigh_comb,my_proc,natoms*c,natomsp,displs,recvcounts)
        call allgatherv_i (neigh_comj,my_proc,natoms*c,natomsp,displs,recvcounts)

        return
        end
