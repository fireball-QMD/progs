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

! eandg.f90
! Program Description
! ===========================================================================
!       This routine calculates the energy and gradient of the functional
! which describes the band-structure energy.  The equation of the
! functional is:
!
!        E = 2*sum{matrix}(H_matrix - eta*S_matrix)
!                         *(2*delta_matrix - S_matrix) + eta*N
!
! For now use OLD functional:
!
!        E = 2*(sum{ii}H_ii - sum{matrix}H_ji*(S_matrix - delta_matrix)
!
! ===========================================================================
! Original Order-N compilation writtten by Spencer Shellman.
!
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
        subroutine eandg (nprows0, ipstart0, ncrows0, icstart0, iteration,       &
     &                    ebs_local, gnorm, N, eta, stop_now)
        use interactions
        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'

! eandg computes the energy and gradient when the coeff matrix = c,
! and computes an optimal step size along the CG (or L-BFGS) direction.
! The energy value returned applies to the value of c at the end of the call;
! c is not advanced until the next call.
! Each processor holds the transpose of a subset of the rows of c, F, Fs, H, S, etc.
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icstart0
        integer, intent (in) :: ipstart0
        integer, intent (in) :: iteration
        integer, intent (in) :: ncrows0
        integer, intent (in) :: nprows0
        integer, intent (in) :: N         ! total number of electrons in system
        real, intent (in) :: eta

! Output
        real, intent (out) :: ebs_local, gnorm
        logical, intent (out) :: stop_now
  
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ichoosesize
        integer ierror
        integer imu, inu, jmu, index, indexH, indexS
        integer nodiv, nomod, nbdiv, nbmod
        integer myrank
        integer, save :: nGmax
        integer nGtmax
        integer nGmax_local
        integer, save :: nHijmax
        integer nHijmax_local
        integer, save :: reset_iter
        integer, save :: nprows, ncrows, ipstart, icstart
        integer, dimension (4) :: rankdep

        integer, dimension (4) :: ichooseA
        integer, dimension (4) :: ichooseB

        double precision ebs1, ebs2
        real eta_f1
        real eta_f2
        real eta_factor
        real, save :: eta_old
        double precision ebs_temp
        double precision ebs_temp2 

        double precision, dimension (2) :: prod(2)
        double precision, dimension (2) :: prodt(2)

        real multA 
        real multB
        real, save :: stepsize, ebs_old
        real gnorm_local

! Holds the C matrix for the next step.
        real, dimension (:, :), allocatable, target, save :: altC_compact_local

        logical need_list
        logical isnan

! Bit matrix indicating which blocks of a product matrix are nonzero.
        logical, dimension (:, :), allocatable, save :: Hij_Sij_bit_matrix, FXSX_bit_matrix, &
             & c_bit_matrix, CX_bit_matrix, HX_SX_bit_matrix

! bit matrices or'ed together to define an MPI graph
        logical, dimension (:,:), allocatable :: accum_bit_matrix

! ****************************************************************************
! Temporary storage vectors. 
! ****************************************************************************
        integer, dimension (:), allocatable :: numTg1_local
        integer, dimension (:), allocatable :: numTg2_local
        integer, dimension (:), allocatable :: numTg3_local
        integer, dimension (:), allocatable :: numTg4_local
        integer, dimension (:, :), allocatable :: listTg1_local
        integer, dimension (:, :), allocatable :: listTg2_local
        integer, dimension (:, :), allocatable :: listTg3_local
        integer, dimension (:, :), allocatable :: listTg4_local

        real, dimension (:, :, :), allocatable, target ::                    &
     &   Tg1_Tg2_Tg3_Tg4_compact_local
        real, dimension (:, :), pointer :: Tg1_compact_local,                &
     &   Tg2_compact_local, Tg3_compact_local, Tg4_compact_local

! ****************************************************************************
! CX matrix declaration - local pieces. Defined by CX = C * X. 
! ****************************************************************************
        integer, dimension (:), allocatable, save :: numCX_local
        integer, dimension (:, :), allocatable, save :: listCX_local

        real, dimension (:, :), allocatable :: CX_compact_local

! ****************************************************************************
! F and Fs matrix declaration - local pieces. Defined by F = h*C and Fs = s*C. 
! ****************************************************************************
        integer, dimension (:), allocatable, save :: numF_local
        integer, dimension (:, :), allocatable, save :: listF_local

        real, dimension (:, :, :), allocatable, target :: F_Fs_compact_local
        real, dimension (:, :), pointer :: F_compact_local, Fs_compact_local

! ****************************************************************************
! FX and FsX matrix declaration - local pieces. Defined by FX = h*C*X and FsX = s*C*X.
! ****************************************************************************
        integer, dimension (:), allocatable, save :: numFX_local
        integer, dimension (:, :), allocatable, save :: listFX_local

        real, dimension (:, :, :), allocatable, target :: FX_FsX_compact_local
        real, dimension (:, :), pointer :: FX_compact_local, FsX_compact_local

! ****************************************************************************
! Hij and Sij matrix declaration - local pieces. Defined by Hij = Ct*F and
! Sij = Ct*Fs.
! ****************************************************************************
        integer, dimension (:), allocatable, target, save :: numHij_local
        integer, dimension (:, :), allocatable, target, save :: listHij_local

        real, dimension (:, :, :), allocatable, target :: Hij_Sij_compact_local
        real, dimension (:, :), pointer :: Hij_compact_local, Sij_compact_local

! ****************************************************************************
! HX and SX matrix declaration - local pieces. Defined by HX = Hij*X and
! SX = Sij*X.
! ****************************************************************************
        integer, dimension (:), allocatable, target, save :: numHX_local
        integer, dimension (:, :), allocatable, target, save :: listHX_local

        real, dimension (:, :, :), allocatable, target :: HX_SX_compact_local
        real, dimension (:, :), pointer :: HX_compact_local, SX_compact_local

! pointers whose targets depend on whether the X matrix is used
        integer, dimension (:), pointer :: numHPtr_local
        integer, dimension (:, :), pointer :: listHPtr_local
        real, dimension (:, :, :), pointer :: HPtr_SPtr_compact_local
        real, dimension (:, :), pointer :: HPtr_compact_local, SPtr_compact_local

        integer, dimension (:), allocatable :: degrees
        integer, dimension (:), allocatable :: edges
        integer edgecount

        integer FXmax, FXmax_local, HXmax, HXmax_local, CXmax, CXmax_local

! X optimization variables
        integer xiteration
        logical x_success, x_stop
        real xgnorm, xtrace, xtrace_new

! Boolean which controls whether the X matrix is used
        logical, parameter :: xuse = .false.

!$ volatile numHij_local, listHij_local, Hij_Sij_compact_local
!$ volatile Hij_compact_local, Sij_compact_local

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

        integer maxiter

! Allocate Arrays
! ============================================================================
 
! Procedure
! ============================================================================
        stop_now = .false.
! save these variables; they may be reordered when the MPI graph is created
        if (iteration.eq.0) then
           nprows = nprows0
           ncrows = ncrows0
           ipstart = ipstart0
           icstart = icstart0
        end if

! Initialize the variable need_list
! On the first iteration this is set to true, so that the index sets of the 
! various matrices may be constructed.  Once constructed they do not change 
! (theoretically).
        need_list = .false.
        if (iteration .eq. 0) then
         need_list = .true.
         if (allocated(altC_compact_local)) deallocate(altC_compact_local)
         allocate (altC_compact_local (ncmax, nprowsmax))
        else
! Here we set c to its new value computed during the last call.
! FIXME inefficient copy
         c_compact_local = altC_compact_local
        end if
! If this is the first iteration or the value of eta has changed, 
! reset the CG accumulator.
        if (iteration .eq. 0 .or. abs(eta_old - eta) .gt. 1.0d-3) reset_iter = 0

! Find out which processor this is.
        call MPI_COMM_RANK (MPI_BTN_WORLD, myrank, ierror)

! sizes of local sections
        nodiv = norbitals/nactualprocs
        nomod = mod(norbitals,nactualprocs)
        nbdiv = nbands/nactualprocs
        nbmod = mod(nbands,nactualprocs)

! Allocate some arrays for transpose of coefficients.
        if (iteration .eq. 0) then 
         if (allocated(numct_local))                                         & 
     &    deallocate (numct_local, listct_local, ct_compact_local)
         allocate (numct_local (ncrowsmax))
         allocate (listct_local (nctmax, ncrowsmax))
         allocate (ct_compact_local (nctmax, ncrowsmax))
         if (allocated(c_bit_matrix)) deallocate(c_bit_matrix)
         allocate (c_bit_matrix (nactualprocs, nactualprocs))
         allocate (accum_bit_matrix (nactualprocs, nactualprocs))
         accum_bit_matrix = hs_bit_matrix
        end if

! Calculate c transpose for the local c matrix contained on this processor.
        call build_transpose (nactualprocs, myrank, c_bit_matrix, iteration.eq.0, &
     &                        ncmax, nctmax, nprows,   &
     &                        nprowsmax, ncrows, ncrowsmax, numc_local,      &
     &                        listc_local, c_compact_local, ipstart,      &
     &                        numct_local, listct_local, ct_compact_local,   &
     &                        icstart, norbitals, nbands)

        if (iteration.eq.0) accum_bit_matrix = accum_bit_matrix .or. c_bit_matrix

! ****************************************************************************
!                            F/Fs   M A T R I C E S
! ****************************************************************************
! Computes local part of the F and Fs matrices from the coefficient vectors.
!       write (*,*) ' Building F and Fs matrices, myrank = ', myrank
        if (iteration .eq. 0) then
         if (allocated (numF_local)) deallocate (numF_local, listF_local)
         allocate (numF_local (nprowsmax))
         allocate (listF_local (nFmax, nprowsmax))
        end if
        allocate (F_Fs_compact_local (nFmax, nprowsmax, 2))
        F_compact_local => F_Fs_compact_local (:, :, 1)
        Fs_compact_local => F_Fs_compact_local (:, :, 2)

! Here we multiply h and s by c to get F and Fs.
        ichoosesize = 2
        ichooseA(1) = 1
        ichooseA(2) = 2
        ichooseB(1) = 1
        ichooseB(2) = 1
        call sendrecv (myrank, iteration, nodiv, nomod, hs_bit_matrix,       &
     &                 .false., numh, listh, h_s_compact, nprows, nhmax,      &
     &                 nprowsmax, ipstart, numc_local, listc_local,          &
     &                 c_compact_local, nprows, ncmax, nprowsmax, ipstart,&
     &                 numF_local, listF_local, F_Fs_compact_local, ichooseA,&
     &                 ichooseB, ichoosesize, nFmax, nprowsmax)


! ****************************************************************************
!                            H/S   M A T R I C E S
! ****************************************************************************
! First Determine the maximum value for the parameter nHijmax.
! Allocate appropriate arrays.
!       write (*,*) ' Building Hij and Sij matrices, myrank = ', myrank
        if (iteration .eq. 0) then
         call set_maxdimension (nactualprocs, myrank, ncrows, nctmax, nFmax, &
     &                          ncrowsmax, nprows, nprowsmax, numct_local,   &
     &                          listct_local, icstart, numF_local,           &
     &                          listF_local, ipstart, 0, nHijmax_local)

         call MPI_ALLREDUCE (nHijmax_local, nHijmax, 1, MPI_INTEGER, MPI_MAX,&
     &                       MPI_BTN_WORLD, ierror)

         if (allocated (numHij_local)) deallocate (numHij_local, listHij_local)
         allocate (numHij_local (ncrowsmax))
         allocate (listHij_local (nHijmax, ncrowsmax))
         if (allocated(Hij_Sij_bit_matrix)) deallocate(Hij_Sij_bit_matrix)
         allocate (Hij_Sij_bit_matrix (nactualprocs, nactualprocs))
        end if
        allocate (Hij_Sij_compact_local (nHijmax, ncrowsmax, 2))
        Hij_compact_local => Hij_Sij_compact_local (:, :, 1)
        Sij_compact_local => Hij_Sij_compact_local (:, :, 2)

! Here we multiply c^T by F and Fs to get H and S.         
        ichoosesize = 2
        ichooseA(1) = 1
        ichooseA(2) = 1
        ichooseB(1) = 1
        ichooseB(2) = 2
! FIXME we should reuse the bit matrix generated when computing C^T from C.
        call sendrecv (myrank, iteration, nodiv, nomod, Hij_Sij_bit_matrix,  &
     &                 iteration.eq.0, numct_local, listct_local, ct_compact_local,  &
     &                 ncrows, nctmax, ncrowsmax, icstart, numF_local,       &
     &                 listF_local, F_Fs_compact_local, nprows, nFmax,       &
     &                 nprowsmax, ipstart, numHij_local, listHij_local,      &
     &                 Hij_Sij_compact_local, ichooseA, ichooseB,            &
     &                 ichoosesize, nHijmax, ncrowsmax)

        if (iteration.eq.0) accum_bit_matrix = accum_bit_matrix .or. Hij_Sij_bit_matrix

        if (xuse) then
! ****************************************************************************
!                            X MATRIX INITIALIZATION
! ****************************************************************************

        if (iteration .eq. 0) then
! FIXME we use nbands for the row dimension; can we do better?
         Xmax = nbands
         if (allocated(numX_local)) deallocate(numX_local,listX_local,X_compact_local)
         allocate (numX_local (ncrowsmax))
         allocate (listX_local (Xmax, ncrowsmax))
         allocate (X_compact_local (Xmax, ncrowsmax))

! initialize to identity matrix
         numX_local = Xmax
         X_compact_local = 0.0d0
         do imu = 1, ncrows
            X_compact_local(imu+icstart-1,imu) = 1.0d0
            do inu = 1, Xmax
               listX_local(inu,imu) = inu
            end do
         end do
        end if
! perform minimization of Tr[(-S)(2X-XSX)]

        if (iteration.ge.0) then
         if (myrank.eq.0) write (*,*) 'Performing X matrix optimization.'
         x_success = .false.
         xiteration = 0
         call xeandg (ncrows, icstart, xiteration, nHijmax, numHij_local, &
              & listHij_local, Sij_compact_local, xtrace, xgnorm, x_stop)
         if (x_stop) then
            x_success = .true.
         else
            if (iteration.eq.0) then
               maxiter = 500
            else
               maxiter = 50
            end if
            do xiteration = 1, maxiter
               call xeandg (ncrows, icstart, xiteration, nHijmax, numHij_local, &
                    & listHij_local, Sij_compact_local, xtrace_new, xgnorm, x_stop)
               if (myrank.eq.0) write (*,201) xiteration, xtrace_new
               if (x_stop) then
                  x_success = .true.
                  xtrace = xtrace_new
                  exit
               end if
               if (abs (xtrace - xtrace_new) .le. ordern_tolerance*abs(xtrace_new)         &
                    &        .and. xgnorm .le. ordern_grad_tolerance) then
                  x_success = .true.
                  xtrace = xtrace_new
                  exit
               end if
               xtrace = xtrace_new
            end do
           end if
          end if
! ****************************************************************************
!                            CX MATRIX
! ****************************************************************************

!       write (*,*) ' Building CX  matrix, myrank = ', myrank
        if (iteration .eq. 0) then
! FIXME is there a faster way to compute CXmax?
         call set_maxdimension (nactualprocs, myrank, nprows, ncmax, Xmax, &
     &                          nprowsmax, ncrows, ncrowsmax, numc_local,   &
     &                          listc_local, ipstart, numX_local,           &
     &                          listX_local, icstart, 1, CXmax_local)

         call MPI_ALLREDUCE (CXmax_local, CXmax, 1, MPI_INTEGER, MPI_MAX,&
     &                       MPI_BTN_WORLD, ierror)

        if (allocated (numCX_local)) deallocate (numCX_local, listCX_local)
         allocate (numCX_local (nprowsmax))
         allocate (listCX_local (CXmax, nprowsmax))
         if (allocated(CX_bit_matrix)) deallocate(CX_bit_matrix)
         allocate (CX_bit_matrix (nactualprocs, nactualprocs))
        end if
        allocate (CX_compact_local (CXmax, nprowsmax))

        ichoosesize = 1
        ichooseA(1) = 1
        ichooseB(1) = 1
        call sendrecv (myrank, iteration, nbdiv, nbmod, CX_bit_matrix,       &
     &                 iteration.eq.0,                                       &
     &                 numc_local, listc_local, c_compact_local, nprows,  &
     &                 ncmax, nprowsmax, ipstart, numX_local,              &
     &                 listX_local, X_compact_local, ncrows, Xmax,&
     &                 ncrowsmax, icstart, numCX_local, listCX_local,      &
     &                 CX_compact_local, ichooseA, ichooseB, ichoosesize,     &
     &                 CXmax, nprowsmax)
        if (iteration.eq.0) accum_bit_matrix = accum_bit_matrix .or. CX_bit_matrix

! ****************************************************************************
!                            FX/FsX   M A T R I C E S
! ****************************************************************************
! Computes local part of the FX and FsX matrices from the coefficient vectors.
!       write (*,*) ' Building FX and FsX matrices, myrank = ', myrank
        if (iteration .eq. 0) then
! FIXME is there a faster way to compute FXmax?
         call set_maxdimension (nactualprocs, myrank, nprows, nhmax, CXmax, &
     &                          nprowsmax, nprows, nprowsmax, numh,   &
     &                          listh, ipstart, numCX_local,           &
     &                          listCX_local, ipstart, 0, FXmax_local)
         call MPI_ALLREDUCE (FXmax_local, FXmax, 1, MPI_INTEGER, MPI_MAX,&
     &                       MPI_BTN_WORLD, ierror)
         if (allocated (numFX_local)) deallocate (numFX_local, listFX_local)
         allocate (numFX_local (nprowsmax))
         allocate (listFX_local (FXmax, nprowsmax))
        end if
        allocate (FX_FsX_compact_local (FXmax, nprowsmax, 2))
        FX_compact_local => FX_FsX_compact_local (:, :, 1)
        FsX_compact_local => FX_FsX_compact_local (:, :, 2)
! Here we multiply h and s by CX to get FX and FsX.
        ichoosesize = 2
        ichooseA(1) = 1
        ichooseA(2) = 2
        ichooseB(1) = 1
        ichooseB(2) = 1
        call sendrecv (myrank, iteration, nodiv, nomod, hs_bit_matrix,       &
     &                 .false., numh, listh, h_s_compact, nprows, nhmax,      &
     &                 nprowsmax, ipstart, numCX_local, listCX_local,          &
     &                 CX_compact_local, nprows, CXmax, nprowsmax, ipstart,&
     &                 numFX_local, listFX_local, FX_FsX_compact_local, ichooseA,&
     &                 ichooseB, ichoosesize, FXmax, nprowsmax)

! ****************************************************************************
!                            HX/SX   M A T R I C E S
! ****************************************************************************
! Allocate appropriate arrays.
!       write (*,*) ' Building HX and SX matrices, myrank = ', myrank
        if (iteration .eq. 0) then
         call set_maxdimension (nactualprocs, myrank, ncrows, nctmax, FXmax, &
     &                          ncrowsmax, nprows, nprowsmax, numct_local,   &
     &                          listct_local, icstart, numFX_local,           &
     &                          listFX_local, ipstart, 0, HXmax_local)

         call MPI_ALLREDUCE (HXmax_local, HXmax, 1, MPI_INTEGER, MPI_MAX,&
     &                       MPI_BTN_WORLD, ierror)

         if (allocated (numHX_local)) deallocate (numHX_local, listHX_local)
         allocate (numHX_local (ncrowsmax))
         allocate (listHX_local (HXmax, ncrowsmax))
         if (allocated(HX_SX_bit_matrix)) deallocate(HX_SX_bit_matrix)
         allocate (HX_SX_bit_matrix (nactualprocs, nactualprocs))
        end if
        allocate (HX_SX_compact_local (HXmax, ncrowsmax, 2))
        HX_compact_local => HX_SX_compact_local (:, :, 1)
        SX_compact_local => HX_SX_compact_local (:, :, 2)

! Here we multiply c^T by FX and FsX to get HX and SX.
        ichoosesize = 2
        ichooseA(1) = 1
        ichooseA(2) = 1
        ichooseB(1) = 1
        ichooseB(2) = 2
! FIXME we should reuse the bit matrix generated when computing C^T from C.
        call sendrecv (myrank, iteration, nodiv, nomod, HX_SX_bit_matrix,  &
     &                 iteration.eq.0, numct_local, listct_local, ct_compact_local,  &
     &                 ncrows, nctmax, ncrowsmax, icstart, numFX_local,       &
     &                 listFX_local, FX_FsX_compact_local, nprows, FXmax,       &
     &                 nprowsmax, ipstart, numHX_local, listHX_local,      &
     &                 HX_SX_compact_local, ichooseA, ichooseB,            &
     &                 ichoosesize, HXmax, ncrowsmax)

        if (iteration.eq.0) accum_bit_matrix = accum_bit_matrix .or. HX_SX_bit_matrix

        end if

        if (xuse) then
           numHptr_local => numHX_local
           listHptr_local => listHX_local
           Hptr_Sptr_compact_local => HX_SX_compact_local
           Hptr_compact_local => HX_compact_local
           Sptr_compact_local => SX_compact_local
        else
           numHptr_local => numHij_local
           listHptr_local => listHij_local
           Hptr_Sptr_compact_local => Hij_Sij_compact_local
           Hptr_compact_local => Hij_compact_local
           Sptr_compact_local => Sij_compact_local
        end if

! Compute the energy value(s) at the current point(s), using the functional 
! formula. We add positive and negative values separately for greater stability.
        ebs1 = 0.0d0
        ebs2 = 0.0d0
        eta_f1 = N
        eta_f2 = 0.0d0
!$omp parallel do private(index) reduction(+:ebs1,ebs2,eta_f1,eta_f2)
        do imu = 1, ncrows
         do inu = 1, numHptr_local(imu)
          index = listHptr_local(inu,imu)
          if (index .eq. imu + icstart - 1) then
           if (Hptr_Sptr_compact_local(inu,imu,1) .ge. 0) then
            ebs1 = ebs1 + 4.0d0*Hptr_Sptr_compact_local(inu,imu,1)
           else
            ebs2 = ebs2 + 4.0d0*Hptr_Sptr_compact_local(inu,imu,1)
           end if
           if (Hptr_Sptr_compact_local(inu,imu,2) .le. 0) then
            eta_f1 = eta_f1 + (-4.0d0)*Hptr_Sptr_compact_local(inu,imu,2)
           else
            eta_f2 = eta_f2 + (-4.0d0)*Hptr_Sptr_compact_local(inu,imu,2)
           end if
          end if
          if (Hptr_Sptr_compact_local(inu,imu,1)                         &
               &         *Hptr_Sptr_compact_local(inu,imu,2) .le. 0) then
           ebs1 = ebs1                                                      &
                &            + (-2.0d0)*Hptr_Sptr_compact_local(inu,imu,1)              &
                &                   *Hptr_Sptr_compact_local(inu,imu,2)
          else
           ebs2 = ebs2                                                      &
                &            + (-2.0d0)*Hptr_Sptr_compact_local(inu,imu,1)              &
                &                   *Hptr_Sptr_compact_local(inu,imu,2)
          end if
          eta_f1 = eta_f1                                                   &
               &             + 2.0d0*Hptr_Sptr_compact_local(inu,imu,2)             &
               &                    *Hptr_Sptr_compact_local(inu,imu,2)
         end do
        end do
        eta_factor = eta_f1 + eta_f2
        ebs_temp = ebs1 + ebs2 + eta * eta_factor

! Compute energy value at this point.
! FIXME I wish I could combine this with the reduction for the energy values 
! at the nearby points, but I have to test whether the energy value has gone 
! up or down first.
        call MPI_ALLREDUCE (ebs_temp, ebs_temp2, 1,              &
     &                      mpi_whatever_double, MPI_SUM, MPI_BTN_WORLD,    &
     &                      ierror)

        ebs_local = ebs_temp2

! Reset CG (take a SD step) if the energy value has increased.
        if (reset_iter .gt. 0 .and. ebs_local .ge. ebs_old) then
         if (myrank.eq.0) write (*,*)                                        &
     &    ' Warning: energy has not decreased, taking a steepest descent step! '
         reset_iter = 0
        end if


! ****************************************************************************
!                     C A L C U L A T E   G R A D I E N T S
! ****************************************************************************
! First Determine the maximum value for the parameter nGmax.
! Allocate appropriate arrays.
!       write (*,*) ' Building gradient matrices, myrank = ', myrank
        if (iteration .eq. 0) then
         if (xuse) then
         call set_maxdimension (nactualprocs, myrank, nprows, FXmax, HXmax,&
     &                          nprowsmax, ncrows, ncrowsmax, numFX_local,    &
     &                          listFX_local, ipstart, numHX_local,          &
     &                          listHX_local, icstart, 1, nGmax_local)
         else
            call set_maxdimension (nactualprocs, myrank, nprows, nFmax, nHijmax, &
                 & nprowsmax, ncrows, ncrowsmax, numF_local, &
                 & listF_local, ipstart, numHij_local, &
                 & listHij_local, icstart, 1, nGmax_local)
         end if

         call MPI_ALLREDUCE (nGmax_local, nGmax, 1, MPI_INTEGER, MPI_MAX,    &
     &                       MPI_BTN_WORLD, ierror)
! take into account all terms of the gradient in determining nGmax
         if (xuse) then
            nGmax = nGmax + FXmax
         else
            nGmax = nGmax + nFmax
         end if

! Allocate the gradient matrix
         if (allocated (numG_local))                                         &
     &    deallocate (numG_local, listG_local, G_compact_local)
         allocate (numG_local (nprowsmax))
         allocate (listG_local (nGmax, nprowsmax))
         allocate (G_compact_local (nGmax, nprowsmax))

         if (allocated(FXSX_bit_matrix)) deallocate(FXSX_bit_matrix)
         allocate (FXSX_bit_matrix (nactualprocs, nactualprocs))

! Allocate the CG step direction, which accumulates from one iteration to 
! the next
         if (allocated (numRG_local))                                        &
     &    deallocate (numRG_local, listRG_local, RG_compact_local)
         allocate (numRG_local (nprowsmax))
         allocate (listRG_local (nGmax, nprowsmax))
         allocate (RG_compact_local (nGmax, nprowsmax))

        end if
        allocate (numTg1_local (nprowsmax))
        allocate (listTg1_local (nGmax, nprowsmax))
        allocate (Tg1_Tg2_Tg3_Tg4_compact_local (nGmax, nprowsmax, 4))
        Tg1_compact_local => Tg1_Tg2_Tg3_Tg4_compact_local (:, :, 1)
        Tg2_compact_local => Tg1_Tg2_Tg3_Tg4_compact_local (:, :, 2)
        Tg3_compact_local => Tg1_Tg2_Tg3_Tg4_compact_local (:, :, 3)
        Tg4_compact_local => Tg1_Tg2_Tg3_Tg4_compact_local (:, :, 4)

! Compute FsXHX, FXSX, FsXSX
        ichoosesize = 3
        ichooseA(1) = 2
        ichooseA(2) = 1
        ichooseA(3) = 2
        ichooseB(1) = 1
        ichooseB(2) = 2
        ichooseB(3) = 2
        if (.not. xuse) then
           call sendrecv (myrank, 0, nbdiv, nbmod, FXSX_bit_matrix, &
                & iteration.eq.0, &
                & numF_local, listF_local, F_Fs_compact_local, nprows, &
                & nFmax, nprowsmax, ipstart, numHij_local, &
                & listHij_local, Hij_Sij_compact_local, ncrows, nHijmax, &
                & ncrowsmax, icstart, numTg1_local, listTg1_local, &
                & Tg1_Tg2_Tg3_Tg4_compact_local, ichooseA, ichooseB, ichoosesize, &
                & nGmax, nprowsmax)
        else
        call sendrecv (myrank, 0, nbdiv, nbmod, FXSX_bit_matrix,       &
     &                 iteration.eq.0,                                       &
     &                 numFX_local, listFX_local, FX_FsX_compact_local, nprows,  &
     &                 FXmax, nprowsmax, ipstart, numHX_local,              &
     &                 listHX_local, HX_SX_compact_local, ncrows, HXmax,&
     &                 ncrowsmax, icstart, numTg1_local, listTg1_local,      &
     &                 Tg1_Tg2_Tg3_Tg4_compact_local, ichooseA, ichooseB, ichoosesize,     &
     &                 nGmax, nprowsmax)
        end if
        if (iteration.eq.0) accum_bit_matrix = accum_bit_matrix .or. FXSX_bit_matrix

! Allocate some temporary arrays.
        allocate (numTg2_local (nprowsmax))
        allocate (listTg2_local (nGmax, nprowsmax))
        allocate (numTg3_local (nprowsmax))
        allocate (listTg3_local (nGmax, nprowsmax))
        allocate (numTg4_local (nprowsmax))
        allocate (listTg4_local (nGmax, nprowsmax))

! Add terms of the gradient expression
        multA = 1.0d0
        multB = 1.0d0
        call sparse_add (nprows, nGmax, nGmax, nGmax, nprowsmax, nprowsmax,  &
     &                   nprowsmax, numTg1_local, listTg1_local,             &
     &                   Tg1_compact_local, multA, ipstart, numTg1_local,    &
     &                   listTg1_local, Tg2_compact_local, multB, ipstart,   &
     &                   numTg4_local, listTg4_local, Tg4_compact_local)

        multA = 1.0d0
        multB = -1.0d0
        if (xuse) then
        call sparse_add (nprows, nGmax, FXmax, nGmax, nprowsmax, nprowsmax,  &
     &                   nprowsmax, numTg1_local, listTg1_local,             &
     &                   Tg3_compact_local, multA, ipstart, numFX_local,      &
     &                   listFX_local, FsX_compact_local, multB, ipstart,      &
     &                   numTg2_local, listTg2_local, Tg2_compact_local)
        else
           call sparse_add (nprows, nGmax, nFmax, nGmax, nprowsmax, nprowsmax, &
                & nprowsmax, numTg1_local, listTg1_local, &
                & Tg3_compact_local, multA, ipstart, numF_local, &
                & listF_local, Fs_compact_local, multB, ipstart, &
                & numTg2_local, listTg2_local, Tg2_compact_local)
        end if

        multA = 1.0d0
        multB = eta
        if (xuse) then
        call sparse_add (nprows, FXmax, nGmax, nGmax, nprowsmax, nprowsmax,  &
     &                   nprowsmax, numFX_local, listFX_local, FX_compact_local,&
     &                   multA, ipstart, numTg2_local, listTg2_local,        &
     &                   Tg2_compact_local, multB, ipstart, numTg3_local,    &
     &                   listTg3_local, Tg3_compact_local)
        else
           call sparse_add (nprows, nFmax, nGmax, nGmax, nprowsmax, nprowsmax, &
                & nprowsmax, numF_local, listF_local, F_compact_local, &
                & multA, ipstart, numTg2_local, listTg2_local, &
                & Tg2_compact_local, multB, ipstart, numTg3_local, &
                & listTg3_local, Tg3_compact_local)
        end if

        multA = 8.0d0
        multB = -4.0d0
        call sparse_add (nprows, nGmax, nGmax, nGmax, nprowsmax, nprowsmax,  &
     &                   nprowsmax, numTg3_local, listTg3_local,             &
     &                   Tg3_compact_local, multA, ipstart, numTg4_local,    &
     &                   listTg4_local, Tg4_compact_local, multB, ipstart,   &
     &                   numTg1_local, listTg1_local, Tg1_compact_local)

        deallocate (numTg3_local, listTg3_local)
        deallocate (numTg4_local, listTg4_local)

! constrain the gradient by the c matrix index set
        call sparse_mask (nprows, nGmax, ncmax, nGmax, nprowsmax, nprowsmax,  &
     &                   nprowsmax, numTg1_local, listTg1_local,             &
     &                   Tg1_compact_local, ipstart, numc_local,    &
     &                   listc_local, ipstart,   &
     &                   numTg2_local, listTg2_local, Tg2_compact_local)

        if (reset_iter .gt. 0) then

! Saves the difference between the old and new gradients.
         multA = 1.0d0
         multB = -1.0d0
         call sparse_add (nprows, nGmax, nGmax, nGmax, nprowsmax, nprowsmax, &
     &                    nprowsmax, numTg2_local, listTg2_local,            &
     &                    Tg2_compact_local, multA, ipstart, numG_local,     &
     &                    listG_local, G_compact_local, multB, ipstart,      &
     &                    numTg1_local, listTg1_local, Tg1_compact_local)

         call sparse_norm2 (nprows, nGmax, nprowsmax, numG_local,            &
              &                      listG_local, G_compact_local, 1, .false., prodt(2))
        end if

! Saves the new gradients for use by the next iteration.
        numG_local = numTg2_local
        listG_local = listTg2_local
        G_compact_local = Tg2_compact_local

        if (reset_iter .gt. 0) then
! Conjugate gradient (Polak-Ribiere; see Nocedal/Wright, Numerical
! Optimization, pp. 120-122)
         call sparse_vecprod (nprows, nGmax, nGmax, nprowsmax, nprowsmax,    &
     &                        numTg1_local, listTg1_local, Tg1_compact_local,&
     &                        1, numTg2_local, listTg2_local,                &
     &                        Tg2_compact_local, 1, .false., prodt(1))

         call MPI_ALLREDUCE (prodt, prod, 2, mpi_whatever_double, MPI_SUM,  &
     &                       MPI_BTN_WORLD, ierror)



! FIXME make sure elsewhere that G and GX are nonzero so the division will be safe
         multA = 1.0d0
         multB = prod(1) / prod(2)
         call sparse_add (nprows, nGmax, nGmax, nGmax, nprowsmax, nprowsmax, &
     &                    nprowsmax, numTg2_local, listTg2_local,            &
     &                    Tg2_compact_local, multA, ipstart, numRG_local,    &
     &                    listRG_local, RG_compact_local, multB, ipstart,    &
     &                    numTg1_local, listTg1_local, Tg1_compact_local)

! FIXME inefficient copy
         multA = 1.0d0
         call sparse_copy (nGmax, nGmax, nprows, nprowsmax, nprows,          &
     &                     nprowsmax, numTg1_local, listTg1_local,           &
     &                     Tg1_compact_local, ipstart, ipstart, 0, .false.,  &
     &                     .false., multA, numRG_local, listRG_local,        &
     &                     RG_compact_local)

        else
! Initialize the CG search direction with the current gradient.
! FIXME inefficient copy
         numRG_local = numTg2_local
         listRG_local = listTg2_local
         RG_compact_local = Tg2_compact_local
        end if

! Calculate the RMS of the gradient
! FIXME the 2-norm is computed earlier; reuse it
        gnorm_local = 0.0d0
!$omp parallel do reduction(+:gnorm_local)
        do imu = 1, nprows
         do inu = 1, numG_local(imu)
          gnorm_local = gnorm_local + G_compact_local(inu,imu)**2
         end do
        end do
        call MPI_ALLREDUCE (gnorm_local, gnorm, 1, mpi_whatever_real,        &
     &                      MPI_SUM, MPI_BTN_WORLD, ierror)
        gnorm = sqrt(gnorm / (norbitals*nbands))
        deallocate (numTg1_local, listTg1_local)
        deallocate (numTg2_local, listTg2_local)
        deallocate (Tg1_Tg2_Tg3_Tg4_compact_local)


! ****************************************************************************
! ****************************************************************************
!                   D E T E R M I N E   S T E P   S I Z E
!                                  A N D 
!        C A L C U L A T E    L O C A L   B A N D - S T R U C T U R E
! ****************************************************************************

        if (xuse) then
        call getstepsize (myrank, iteration, reset_iter, nprows, ipstart, ncrows, icstart, N,  &
             &                         eta, nGmax, CXmax, HXmax, FXmax,  &
             &                         numCX_local, listCX_local, numFX_local, listFX_local,  &
             &                         numHX_local, listHX_local, numct_local,     &
             &                         listct_local, ct_compact_local, ebs_local,    &
             &                         stepsize, stop_now, xuse)
        else
           call getstepsize (myrank, iteration, reset_iter, nprows, ipstart, ncrows, icstart, N, &
                & eta, nGmax, ncmax, nHijmax, nFmax, &
                & numc_local, listc_local, numF_local, listF_local, &
                & numHij_local, listHij_local, numct_local, &
                & listct_local, ct_compact_local, ebs_local, &
                & stepsize, stop_now, xuse)
        end if
        deallocate (Hij_Sij_compact_local)
        deallocate (F_Fs_compact_local)
        if (xuse) then
        deallocate (HX_SX_compact_local)
        deallocate (FX_FsX_compact_local)
        deallocate (CX_compact_local)
        end if

        if (.not. stop_now) then
! FIXME: Maybe we should change sparse_add so we can avoid changing the index 
! sets.
         allocate (numTg1_local (nprowsmax))
         allocate (listTg1_local (ncmax, nprowsmax))

! Saves the changed coefficient vectors and X matrix.  The coeff vectors and X matrix
! are not changed until the next iteration.
! We assume here that the gradient has the same or less sparsity as C/X.
         multA = 1.0d0
         multB = stepsize
         call sparse_add (nprows, ncmax, nGmax, ncmax, nprowsmax, nprowsmax, &
     &                    nprowsmax, numc_local, listc_local,                &
     &                    c_compact_local, multA, ipstart,                   &
     &                    numRG_local, listRG_local,                           &
     &                    RG_compact_local, multB, ipstart,   &
     &                    numTg1_local, listTg1_local, altC_compact_local)

         deallocate (numTg1_local, listTg1_local)
        end if

! Save this value of eta to compare it later.
        eta_old = eta

! Save the energy value to compare it later.
        ebs_old = ebs_local

! Update the CG iteration index.
! FIXME: Do we want to periodically set reset_iter to zero?
        reset_iter = reset_iter + 1

        if (iteration.eq.0) then
! create MPI graph to optimize communications
           allocate(degrees(nactualprocs),edges(nactualprocs*nactualprocs))
           edges = 0
           edgecount=0
           degrees(1) = 0
           do imu = 1,nactualprocs
              if (imu.gt.1) degrees(imu) = degrees(imu-1)
              do inu = 1,nactualprocs
                 if (imu.ne.inu .and. (accum_bit_matrix(imu,inu) .or. accum_bit_matrix(inu,imu))) then
                    edgecount = edgecount + 1
                    edges(edgecount) = inu-1
                    degrees(imu) = degrees(imu)+1
                 end if
              end do
           end do
           call mpi_graph_create(MPI_BTN_WORLD, nactualprocs, degrees, edges, &
                & .true., MPI_OPT_WORLD, ierror)
           deallocate(degrees,edges)

! reorder matrix segments according to the optimized graph
! this needs to be done to saved index sets & element sets
           call reorder_matrix (nhmax, nprowsmax, &
                & numh, listh, h_s_compact,&
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           call reorder_matrix_noindex (nhmax, nprowsmax, &
                & h_s_compact(:,:,2),&
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           call reorder_matrix (ncmax, nprowsmax, &
                & numc_local, listc_local, c_compact_local,&
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           if (xuse) call reorder_matrix (Xmax, ncrowsmax, &
                & numX_local, listX_local, X_compact_local,&
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
! FIXME is this call necessary?  Will C^T be used later?
           call reorder_matrix (nctmax, ncrowsmax, &
                & numct_local, listct_local, ct_compact_local,&
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           call reorder_matrix_noindex (ncmax, nprowsmax, &
                & altC_compact_local,&
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           call reorder_matrix (nGmax, nprowsmax, &
                & numG_local, listG_local, G_compact_local,&
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           call reorder_matrix_indexonly (nFmax, nprowsmax, &
                & numF_local, listF_local, &
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           call reorder_matrix_indexonly (nHijmax, ncrowsmax, &
                & numHij_local, listHij_local, &
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           if (xuse) call reorder_matrix_indexonly (CXmax, ncrowsmax, &
                & numCX_local, listCX_local, &
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           if (xuse) call reorder_matrix_indexonly (FXmax, ncrowsmax, &
                & numFX_local, listFX_local, &
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
           if (xuse) call reorder_matrix_indexonly (HXmax, ncrowsmax, &
                & numHX_local, listHX_local, &
                & MPI_BTN_WORLD, MPI_OPT_WORLD)
! Set the CG accumulator to the reordered gradient.
           numRG_local = numG_local
           listRG_local = listG_local
           RG_compact_local = G_compact_local

! reorder rank-dependent variables
           rankdep(1) = ipstart
           rankdep(2) = icstart
           rankdep(3) = nprows
           rankdep(4) = ncrows
           call reorder_int_array (rankdep, 4, MPI_BTN_WORLD, MPI_OPT_WORLD)
           ipstart = rankdep(1)
           icstart = rankdep(2)
           nprows = rankdep(3)
           ncrows = rankdep(4)

! save the old communicator and replace it
           MPI_BTN_WORLD_SAVE = MPI_BTN_WORLD
           MPI_BTN_WORLD = MPI_OPT_WORLD
        end if

! Deallocate Arrays
! ===========================================================================

        if (iteration.eq.0) deallocate(accum_bit_matrix)

! Format Statements
! ===========================================================================
	 
201     format (2x, ' Processor    0: Iteration ', i5, ' rank = ', f18.6) 
 
        return
        end subroutine eandg




! cleans up stuff done by eandg
! MUST be called after final eandg iteration!!!

      subroutine eandg_cleanup ()

        use ordern
        use interactions
        implicit none

        include 'mpif.h'

        integer ierror

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! reorder matrix segments to the original order
! FIXME will numh/listh/h_s be needed at all?
        call reorder_matrix (nhmax, nprowsmax, &
             & numh, listh, h_s_compact,&
             & MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE)
        call reorder_matrix_noindex (nhmax, nprowsmax, &
             & h_s_compact(:,:,2),&
             & MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE)
        call reorder_matrix (ncmax, nprowsmax, &
             & numc_local, listc_local, c_compact_local,&
             & MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE)

! restore old communicator
        call mpi_comm_free(MPI_OPT_WORLD,ierror)
        MPI_BTN_WORLD = MPI_BTN_WORLD_SAVE

        return
      end subroutine eandg_cleanup



! Reorders the distributed segments of a matrix when the MPI communicator is changed
! from comm_from to comm_to, without the index sets.
! Processor k gets the segment from processor l under the
! old ranks, where k is the old rank and l is the new rank.
        subroutine reorder_matrix_noindex (ncolsmax, nrowsmax, &
     &                              M_compact_local,&
     &                              comm_from, comm_to)

        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: comm_from, comm_to
        integer, intent (in) :: ncolsmax
        integer, intent (in) :: nrowsmax

        real, intent (inout), dimension (ncolsmax, nrowsmax) :: M_compact_local

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror, request
        integer rank_from,rank_to,recv_from,send_to
! MPI receive status
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer group_from,group_to
        real, dimension (:,:), allocatable :: M_compact_temp

! Allocate Arrays
! ============================================================================
 
! Procedure
! ============================================================================

        allocate(M_compact_temp(ncolsmax,nrowsmax))
        M_compact_temp = M_compact_local

        call MPI_COMM_RANK (comm_from, rank_from, ierror)
        call MPI_COMM_RANK (comm_to, rank_to, ierror)
        call mpi_comm_group(comm_from,group_from,ierror)
        call mpi_comm_group(comm_to,group_to,ierror)
        call mpi_group_translate_ranks(group_to,1,rank_from,&
             &group_from,send_to,ierror)
        recv_from = rank_to
! don't send if this processor hasn't been reassigned
        if (send_to.ne.rank_from .and. recv_from.ne.rank_from) then
! FIXME should we pack first?
           call mpi_isend(M_compact_temp,ncolsmax*nrowsmax,mpi_whatever_real,send_to,0,&
                &comm_from,request,ierror)
           call mpi_recv (M_compact_local,ncolsmax*nrowsmax,mpi_whatever_real,recv_from,0,&
                &comm_from,status,ierror)
           call MPI_WAIT (request, status, ierror)
        end if

! Deallocate Arrays
! ===========================================================================
        deallocate(M_compact_temp)

! Format Statements
! ===========================================================================
 
        return
        end



! Reorders the distributed segments of a matrix when the MPI communicator is changed
! from comm_from to comm_to.  Processor k gets the segment from processor l under the
! old ranks, where k is the old rank and l is the new rank.
        subroutine reorder_matrix (ncolsmax, nrowsmax, &
     &                              numM_local, listM_local, M_compact_local,&
     &                              comm_from, comm_to)

        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: comm_from, comm_to
        integer, intent (in) :: ncolsmax
        integer, intent (in) :: nrowsmax

        integer, intent (inout), dimension (nrowsmax) :: numM_local
        integer, intent (inout), dimension (ncolsmax, nrowsmax) :: listM_local

        real, intent (inout), dimension (ncolsmax, nrowsmax) :: M_compact_local

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror, request1, request2, request3
        integer rank_from,rank_to,recv_from,send_to
! MPI receive status
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer group_from,group_to
        integer, dimension (:), allocatable :: numM_temp
        integer, dimension (:,:), allocatable :: listM_temp
        real, dimension (:,:), allocatable :: M_compact_temp

! Allocate Arrays
! ============================================================================
 
! Procedure
! ============================================================================

        allocate(numM_temp(nrowsmax),listM_temp(ncolsmax,nrowsmax),M_compact_temp(ncolsmax,nrowsmax))
        numM_temp = numM_local
        listM_temp = listM_local
        M_compact_temp = M_compact_local

        call MPI_COMM_RANK (comm_from, rank_from, ierror)
        call MPI_COMM_RANK (comm_to, rank_to, ierror)
        call mpi_comm_group(comm_from,group_from,ierror)
        call mpi_comm_group(comm_to,group_to,ierror)
        call mpi_group_translate_ranks(group_to,1,rank_from,&
             &group_from,send_to,ierror)
        recv_from = rank_to
! don't send if this processor hasn't been reassigned
        if (send_to.ne.rank_from .and. recv_from.ne.rank_from) then
! FIXME should we pack first?
           call mpi_isend(numM_temp,nrowsmax,mpi_integer,send_to,0,&
                &comm_from,request1,ierror)
           call mpi_isend(listM_temp,ncolsmax*nrowsmax,mpi_integer,send_to,0,&
                &comm_from,request2,ierror)
           call mpi_isend(M_compact_temp,ncolsmax*nrowsmax,mpi_whatever_real,send_to,0,&
                &comm_from,request3,ierror)
           call mpi_recv (numM_local,nrowsmax,mpi_integer,recv_from,0,&
                &comm_from,status,ierror)
           call mpi_recv (listM_local,ncolsmax*nrowsmax,mpi_integer,recv_from,0,&
                &comm_from,status,ierror)
           call mpi_recv (M_compact_local,ncolsmax*nrowsmax,mpi_whatever_real,recv_from,0,&
                &comm_from,status,ierror)
           call MPI_WAIT (request1, status, ierror)
           call MPI_WAIT (request2, status, ierror)
           call MPI_WAIT (request3, status, ierror)
        end if

! Deallocate Arrays
! ===========================================================================
        deallocate(numM_temp,listM_temp,M_compact_temp)

! Format Statements
! ===========================================================================
 
        return
        end



! Reorders the distributed segments of the index sets of a matrix when the MPI communicator is changed
! from comm_from to comm_to.  Processor k gets the segment from processor l under the
! old ranks, where k is the old rank and l is the new rank.
        subroutine reorder_matrix_indexonly (ncolsmax, nrowsmax, &
     &                              numM_local, listM_local, &
     &                              comm_from, comm_to)

        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: comm_from, comm_to
        integer, intent (in) :: ncolsmax
        integer, intent (in) :: nrowsmax

        integer, intent (inout), dimension (nrowsmax) :: numM_local
        integer, intent (inout), dimension (ncolsmax, nrowsmax) :: listM_local

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror, request1, request2
        integer rank_from,rank_to,recv_from,send_to
! MPI receive status
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer group_from,group_to
        integer, dimension (:), allocatable :: numM_temp
        integer, dimension (:,:), allocatable :: listM_temp

! Allocate Arrays
! ============================================================================
 
! Procedure
! ============================================================================

        allocate(numM_temp(nrowsmax),listM_temp(ncolsmax,nrowsmax))
        numM_temp = numM_local
        listM_temp = listM_local

        call MPI_COMM_RANK (comm_from, rank_from, ierror)
        call MPI_COMM_RANK (comm_to, rank_to, ierror)
        call mpi_comm_group(comm_from,group_from,ierror)
        call mpi_comm_group(comm_to,group_to,ierror)
        call mpi_group_translate_ranks(group_to,1,rank_from,&
             &group_from,send_to,ierror)
        recv_from = rank_to
! don't send if this processor hasn't been reassigned
        if (send_to.ne.rank_from .and. recv_from.ne.rank_from) then
! FIXME should we pack first?
           call mpi_isend(numM_temp,nrowsmax,mpi_integer,send_to,0,&
                &comm_from,request1,ierror)
           call mpi_isend(listM_temp,ncolsmax*nrowsmax,mpi_integer,send_to,0,&
                &comm_from,request2,ierror)
           call mpi_recv (numM_local,nrowsmax,mpi_integer,recv_from,0,&
                &comm_from,status,ierror)
           call mpi_recv (listM_local,ncolsmax*nrowsmax,mpi_integer,recv_from,0,&
                &comm_from,status,ierror)
           call MPI_WAIT (request1, status, ierror)
           call MPI_WAIT (request2, status, ierror)
        end if

! Deallocate Arrays
! ===========================================================================
        deallocate(numM_temp,listM_temp)

! Format Statements
! ===========================================================================
 
        return
        end



! Reorders a distributed integer array when the MPI communicator is changed
! from comm_from to comm_to.  Processor k gets the segment from processor l under the
! old ranks, where k is the old rank and l is the new rank.
        subroutine reorder_int_array (array_seg, array_len, comm_from, comm_to)

        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: array_len, comm_from, comm_to
        integer, intent (inout), dimension (array_len) :: array_seg

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror, request
        integer rank_from,rank_to,recv_from,send_to
! MPI receive status
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer group_from,group_to
        integer, dimension (:), allocatable :: itemp

! Allocate Arrays
! ============================================================================
 
! Procedure
! ============================================================================

        allocate(itemp(array_len))
        itemp = array_seg

        call MPI_COMM_RANK (comm_from, rank_from, ierror)
        call MPI_COMM_RANK (comm_to, rank_to, ierror)
        call mpi_comm_group(comm_from,group_from,ierror)
        call mpi_comm_group(comm_to,group_to,ierror)
        call mpi_group_translate_ranks(group_to,1,rank_from,&
             &group_from,send_to,ierror)
        recv_from = rank_to
        if (send_to.ne.rank_from .and. recv_from.ne.rank_from) then
           call mpi_isend(itemp,array_len,mpi_integer,send_to,0,&
                &comm_from,request,ierror)
           call mpi_recv (array_seg,array_len,mpi_integer,recv_from,0,&
                &comm_from,status,ierror)
           call MPI_WAIT (request, status, ierror)
        end if

! Deallocate Arrays
! ===========================================================================

        deallocate(itemp)

! Format Statements
! ===========================================================================

        return
        end
