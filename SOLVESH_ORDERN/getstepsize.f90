! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
! University of Utah - James P. Lewis, Chair
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
! Ohio State University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! getstepsize.f90
! Program Description
! ===========================================================================
! Determine the stepsize using a line search.
! ===========================================================================
! Original Order-N compilation by Spencer Shellman
! Code rewritten by:
! James P. Lewis 
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

       recursive subroutine getstepsize (myrank, iteration, reset_iter, &
           & nprows, ipstart, ncrows, icstart, N, eta, &
           & nGmax, CXmax, HXmax, FXmax, &
           & numCX, listCX, numFX, listFX, numHX, listHX, &
           & numct, listct, ct_compact, &
           & ebs_local, stepsize, stop_now, xuse)
        use interactions
        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'

! The parameters represent variables found in eandg.
! Here ebs_local and stepsize are input parameters, as well as stepsize being a output parameter.
! ebs_local indicates the energy value at the current point, while stepsize is used to choose
! nearby points to evaluate energy values and interpolate.  Generally stepsize will be one of
! the proposed stepsizes returned by the previous call to getstepsize.
! This stepsize should be chosen carefully; if it is too large the interpolation becomes unstable.
! The increased precision of this stepsize function comes as a result of computing the energy
! at four nearby points (the value of stepsize determines how near) and interpolating to obtain
! a fourth-degree polynomial.
! The first iteration performs a line search instead of using the 4th-degree polynomial method.
! If reset_iter = 0 then a line search is performed.  This indicates that the energy did not
! decrease.  If the leading coefficient of the polynomial is zero then a line search is done.
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: myrank, iteration, reset_iter
        integer, intent (in) :: icstart
        integer, intent (in) :: ipstart
        integer, intent (in) :: ncrows
        integer, intent (in) :: nprows
        integer, intent (in) :: nGmax
        integer, intent (in) :: CXmax, HXmax, FXmax
        integer, intent (in) :: N         ! total number of electrons in system
        real, intent (in) :: eta, ebs_local

        integer, dimension (ncrowsmax), intent (in) :: numFX
        integer, dimension (FXmax, nprowsmax), intent (in) :: listFX

        integer, dimension (ncrowsmax), intent (in) :: numHX
        integer, dimension (HXmax, ncrowsmax), intent (in) :: listHX

        integer, dimension (ncrowsmax), intent (in) :: numCX
        integer, dimension (CXmax, nprowsmax), intent (in) :: listCX

        integer, dimension (ncrowsmax), intent (in) :: numct
        integer, dimension (nctmax, ncrowsmax), intent (in) :: listct
        real, dimension (nctmax, ncrowsmax), intent (in) :: ct_compact

! use the X matrix?
        logical, intent (in) :: xuse

! Output
        real, intent (inout) :: stepsize
        logical, intent (out) :: stop_now

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================

        integer nGtmax
        integer ierror
        integer imu, inu, jmu, jnu
        integer index, indexH, indexS
        integer chooseA(8), chooseB(8)
        integer nodiv, nomod, nbdiv, nbmod
        integer qr_retcode

! ****************************************************************************
! Step directions (Rt) transpose matrix declaration - local pieces.
! ****************************************************************************
        integer, dimension (:), allocatable :: numRGt_local       ! RG transpose
        integer, dimension (:, :), allocatable :: listRGt_local   ! RG transpose
        real, dimension (:, :), allocatable :: RGt_compact_local  ! RG transpose

        integer, dimension (:, :), allocatable :: numcp_local
        integer, dimension (:, :, :), allocatable :: listcp_local
        real, dimension (:, :, :), allocatable :: cp_compact_local

        integer, dimension (:, :), allocatable :: numcpt_local
        integer, dimension (:, :, :), allocatable :: listcpt_local
        real, dimension (:, :, :), allocatable :: cpt_compact_local

        real, dimension (:, :, :), allocatable :: CXp_compact_local

        real, dimension (:, :, :, :), allocatable :: FX_FsX_compact_local
        real, dimension (:, :, :, :), allocatable :: HX_SX_compact_local

!$ volatile numHX,listHX,HX_SX_compact_local

        real multA, multB
        double precision eta_factor      ! value mult'ed by eta for each polynomial coeff.
        double precision ebs_temp(4), ebs(4), ebs1, ebs2, eta_f1, eta_f2
        real prev_step, next_step, prev_ebs
        logical, dimension (:, :), allocatable, save :: HX_SX_bit_matrix, RGt_bit_matrix, &
             & CX_bit_matrix
        logical larger_step

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! QR algorithm variables
        double precision ebs_roots(3)
        double precision qr_detil
        double precision, dimension (3) :: ai_roots
        double precision, dimension (3) :: coeff
        double precision, dimension (4) :: coeff_new
        integer,          dimension (3) :: qr_work_area
        double precision, dimension (3, 3) :: qr_work_matrix
        double precision, dimension (3) :: r_roots
        double precision, dimension (4,4) :: matA
        double precision, dimension (4) :: vecB


! Code begins
! ===========================================================================

        stop_now = .false.

! sizes of local sections
        nodiv = norbitals / nactualprocs
        nomod = mod(norbitals,nactualprocs)
        nbdiv = nbands / nactualprocs
        nbmod = mod(nbands,nactualprocs)

! Allocate some arrays here.
        if (iteration.eq.0) then
           if (allocated(HX_SX_bit_matrix)) deallocate(HX_SX_bit_matrix)
           allocate (HX_SX_bit_matrix (nactualprocs, nactualprocs))
           if (allocated(RGt_bit_matrix)) deallocate(RGt_bit_matrix)
           allocate (RGt_bit_matrix (nactualprocs, nactualprocs))
           if (xuse) then
           if (allocated(CX_bit_matrix)) deallocate(CX_bit_matrix)
           allocate (CX_bit_matrix (nactualprocs, nactualprocs))
           end if
        end if

        nGtmax = 2*nFtmax
        allocate (numRGt_local (ncrowsmax))
        allocate (listRGt_local (nGtmax, ncrowsmax))
        allocate (RGt_compact_local (nGtmax, ncrowsmax))

! Compute the transpose of the step direction matrix pieces.
! Calculate RG transpose for the local RG matrix contained on this processor.
!       write (*,*) ' Transpose the RG matrix, myrank = ', myrank
        call build_transpose (nactualprocs, myrank, RGt_bit_matrix, iteration.eq.0,   &
             &                        nGmax, nGtmax, nprows, nprowsmax, ncrows,      &
             &                        ncrowsmax, numRG_local, listRG_local,            &
             &                        RG_compact_local, ipstart, numRGt_local,         &
             &                        listRGt_local, RGt_compact_local, icstart,       &
             &                        norbitals, nbands)

        if (reset_iter .gt. 0) then
! 4th-degree polynomial method; assumes that a stepsize is available.
           allocate (numcp_local (nprowsmax, 4))
           allocate (listcp_local (ncmax, nprowsmax, 4))
           allocate (cp_compact_local (ncmax, nprowsmax, 4))

           allocate (numcpt_local (ncrowsmax, 4))
           allocate (listcpt_local (nctmax, ncrowsmax, 4))
           allocate (cpt_compact_local (nctmax, ncrowsmax, 4))

           if (xuse) then
              allocate (CXp_compact_local (CXmax, nprowsmax, 4))
           end if

           allocate (FX_FsX_compact_local (FXmax, nprowsmax, 2, 4))
           allocate (HX_SX_compact_local (HXmax, ncrowsmax, 2, 4))

           do jmu = 1,4
              multA = 1.0d0
              if (jmu.eq.1) then
                 multB = stepsize
              elseif (jmu.eq.2) then
                 multB = 2.0d0*stepsize
              elseif (jmu.eq.3) then
                 multB = -stepsize
              elseif (jmu.eq.4) then
                 multB = -2.0d0*stepsize
              end if
              call sparse_add (nprows, ncmax, nGmax, ncmax, nprowsmax, nprowsmax, &
                   &                    nprowsmax, numc_local, listc_local,              &
                   &                    c_compact_local, multA, ipstart,                  &
                   &                    numRG_local,                          &
                   &                    listRG_local,                       &
                   &                    RG_compact_local, multB, ipstart,   &
                   &                    numcp_local, listcp_local, cp_compact_local(:,:,jmu))

              call sparse_add (ncrows, nctmax, nGtmax, nctmax, ncrowsmax, ncrowsmax, &
                   &                    ncrowsmax, numct, listct,              &
                   &                    ct_compact, multA, icstart,                  &
                   &                    numRGt_local,                          &
                   &                    listRGt_local,                       &
                   &                    RGt_compact_local, multB, icstart,   &
                   &                    numcpt_local, listcpt_local, cpt_compact_local(:,:,jmu))
           end do

           if (xuse) then
              chooseA(1) = 1
              chooseB(1) = 1
              chooseA(2) = 2
              chooseB(2) = 1
              chooseA(3) = 3
              chooseB(3) = 1
              chooseA(4) = 4
              chooseB(4) = 1

              call sendrecv (myrank, 1, nbdiv, nbmod, CX_bit_matrix,       &
                   &                 iteration.eq.0,                                       &
                   &                 numcp_local, listcp_local, cp_compact_local, nprows,  &
                   &                 ncmax, nprowsmax, ipstart, numX_local,              &
                   &                 listX_local, X_compact_local, ncrows, Xmax,&
                   &                 ncrowsmax, icstart, numCX, listCX,      &
                   &                 CXp_compact_local, chooseA, chooseB, 4,     &
                   &                 CXmax, nprowsmax)
           end if

           chooseA(1) = 1
           chooseA(2) = 2
           chooseB(1) = 1
           chooseB(2) = 1
           chooseA(3) = 1
           chooseA(4) = 2
           chooseB(3) = 2
           chooseB(4) = 2
           chooseA(5) = 1
           chooseA(6) = 2
           chooseB(5) = 3
           chooseB(6) = 3
           chooseA(7) = 1
           chooseA(8) = 2
           chooseB(7) = 4
           chooseB(8) = 4

           if (xuse) then
              call sendrecv (myrank, 1, nodiv, nomod, hs_bit_matrix, .false., &
                   & numh, listh, h_s_compact, nprows, nhmax, nprowsmax, ipstart, &
                   & numCX, listCX, CXp_compact_local, nprows, CXmax, nprowsmax, ipstart, &
                   & numFX, listFX, FX_FsX_compact_local, chooseA, chooseB, 8, FXmax, nprowsmax)
           else
              call sendrecv (myrank, 1, nodiv, nomod, hs_bit_matrix, .false., &
                   & numh, listh, h_s_compact, nprows, nhmax, nprowsmax, ipstart, &
                   & numcp_local, listcp_local, cp_compact_local, nprows, ncmax, nprowsmax, ipstart, &
                   & numFX, listFX, FX_FsX_compact_local, chooseA, chooseB, 8, FXmax, nprowsmax)
           end if
           chooseA(1) = 1
           chooseA(2) = 1
           chooseB(1) = 1
           chooseB(2) = 2
           chooseA(3) = 2
           chooseA(4) = 2
           chooseB(3) = 3
           chooseB(4) = 4
           chooseA(5) = 3
           chooseA(6) = 3
           chooseB(5) = 5
           chooseB(6) = 6
           chooseA(7) = 4
           chooseA(8) = 4
           chooseB(7) = 7
           chooseB(8) = 8

           call sendrecv (myrank, 1, nodiv, nomod, HX_SX_bit_matrix, iteration.eq.0, &
                & numcpt_local, listcpt_local, cpt_compact_local, ncrows, nctmax, ncrowsmax, icstart, &
                & numFX, listFX, FX_FsX_compact_local, nprows, FXmax, nprowsmax, ipstart, &
                & numHX, listHX, HX_SX_compact_local, chooseA, chooseB, 8, HXmax, &
                & ncrowsmax)

           do jmu = 1,4
! Add positive and negative values separately for more stability.
              indexH = 1
              indexS = 2
              ebs1 = 0.0d0
              ebs2 = 0.0d0
              eta_f1 = N
              eta_f2 = 0.0d0
!$omp parallel do private(index) reduction(+:ebs1,ebs2,eta_f1,eta_f2)
              do imu = 1, ncrows
                 do inu = 1, numHX(imu)
                    index = listHX(inu,imu)
                    if (index .eq. imu + icstart - 1) then
                       if (HX_SX_compact_local(inu,imu,indexH,jmu) .ge. 0) then
                          ebs1 = ebs1 + 4.0d0*HX_SX_compact_local(inu,imu,indexH,jmu)
                       else
                          ebs2 = ebs2 + 4.0d0*HX_SX_compact_local(inu,imu,indexH,jmu)
                       end if
                       if (HX_SX_compact_local(inu,imu,indexS,jmu) .le. 0) then
                          eta_f1 = eta_f1 + (-4.0d0)*HX_SX_compact_local(inu,imu,indexS,jmu)
                       else
                          eta_f2 = eta_f2 + (-4.0d0)*HX_SX_compact_local(inu,imu,indexS,jmu)
                       end if
                    end if
                    if (HX_SX_compact_local(inu,imu,indexH,jmu)*HX_SX_compact_local(inu,imu,indexS,jmu) .le. 0) then
                       ebs1 = ebs1                                                &
                            &     + (-2.0d0)*HX_SX_compact_local(inu,imu,indexH,jmu)*HX_SX_compact_local(inu,imu,indexS,jmu)
                    else
                       ebs2 = ebs2                                                &
                            &     + (-2.0d0)*HX_SX_compact_local(inu,imu,indexH,jmu)*HX_SX_compact_local(inu,imu,indexS,jmu)
                    end if
                    eta_f1 = eta_f1                                            &
                         &     + 2.0d0*HX_SX_compact_local(inu,imu,indexS,jmu)*HX_SX_compact_local(inu,imu,indexS,jmu)
                 end do
              end do
              eta_factor = eta_f1 + eta_f2
              ebs_temp(jmu) = ebs1 + ebs2 + eta * eta_factor
           end do
           call MPI_ALLREDUCE (ebs_temp, ebs, 4, mpi_whatever_double,       &
                &                      MPI_SUM, MPI_BTN_WORLD, ierror)

! Compute the coefficients of the interpolating polynomial.
           matA(1,1) = stepsize
           matA(1,2) = matA(1,1) * stepsize
           matA(1,3) = matA(1,2) * stepsize
           matA(1,4) = matA(1,3) * stepsize

           matA(2,1) = matA(1,1) * 2.0d0
           matA(2,2) = matA(1,2) * 4.0d0
           matA(2,3) = matA(1,3) * 8.0d0
           matA(2,4) = matA(1,4) * 16.0d0

           matA(3,1) = -matA(1,1)
           matA(3,2) = matA(1,2)
           matA(3,3) = -matA(1,3)
           matA(3,4) = matA(1,4)

           matA(4,1) = -matA(2,1)
           matA(4,2) = matA(2,2)
           matA(4,3) = -matA(2,3)
           matA(4,4) = matA(2,4)

           vecB = ebs - ebs_local

           call NGAUSS(4,matA,4,vecB,coeff_new)

! If the leading coefficient is zero then do a line search instead.
           if (coeff_new(4) .eq. 0.0d0) then
              if (myrank.eq.0) write (*,*) 'Warning: Leading coefficient is zero, executing line search'
              call getstepsize (myrank, iteration, 0, &
                   & nprows, ipstart, ncrows, icstart, N, eta, &
                   & nGmax, CXmax, HXmax, FXmax, &
                   & numCX, listCX, numFX, listFX, numHX, listHX, &
                   & numct, listct, ct_compact, &
                   & ebs_local, stepsize, stop_now, xuse)
           else
! 3rd-degree derivative
              coeff(1) = 0.25d0*coeff_new(1)/coeff_new(4)
              coeff(2) = 0.5d0*coeff_new(2)/coeff_new(4)
              coeff(3) = 0.75d0*coeff_new(3)/coeff_new(4)

! FIXME is there a better way to compute cubic roots?
              call qr_algeq_solver (3, coeff, epsilon(r_roots), r_roots, ai_roots,&
                   &                         qr_detil, qr_work_matrix, qr_work_area,       &
                   &                         qr_retcode)
              if (qr_retcode.ne.0) then
                 stop_now = .true.
                 return
              end if

! We choose the real root that yields the greatest decrease.
              if (ai_roots(1).ne.0 .or. ai_roots(2).ne.0) then
! only one real root
                 do jnu = 1, 3
                    if (ai_roots(jnu).eq.0) then
                       jmu = jnu
                       exit
                    end if
                 end do
              else
! compute polynomial at roots
                 do jnu = 1,3
                    ebs_roots(jnu) = (((coeff_new(4)*r_roots(jnu) + coeff_new(3)) &
                         & *r_roots(jnu) + coeff_new(2))*r_roots(jnu) &
                         & + coeff_new(1))*r_roots(jnu)
                 end do
                 jmu = 1
                 do jnu = 2, 3
                    if ((r_roots(jmu).ge.0.and.r_roots(jnu).lt.0) .or. ebs_roots(jnu).lt.ebs_roots(jmu)) jmu = jnu
                 end do
              end if
              stepsize = r_roots(jmu)
! The stepsize should be negative, otherwise execute a line search
              if (stepsize .ge. 0.0d0) then
                 if (myrank.eq.0) write (*,*) 'Warning: stepsize is nonnegative, executing line search'
                 call getstepsize (myrank, iteration, 0, &
                      & nprows, ipstart, ncrows, icstart, N, eta, &
                      & nGmax, CXmax, HXmax, FXmax, &
                      & numCX, listCX, numFX, listFX, numHX, listHX, &
                      & numct, listct, ct_compact, &
                      & ebs_local, stepsize, stop_now, xuse)
              end if
           end if
        else
! line search; no stepsize is available, we must compute one

! FIXME It would be more efficient to pass to sendrecv a list (i1,...,in) and R, and let it
! perform the multiplication with c+i1*R,...,c+in*R, rather than have to send c+i1*R,...,c+in*R
! to the various processors.
           allocate (numcp_local (nprowsmax, 1))
           allocate (listcp_local (ncmax, nprowsmax, 1))
           allocate (cp_compact_local (ncmax, nprowsmax, 1))

           allocate (numcpt_local (ncrowsmax, 1))
           allocate (listcpt_local (nctmax, ncrowsmax, 1))
           allocate (cpt_compact_local (nctmax, ncrowsmax, 1))

           if (xuse) then
              allocate (CXp_compact_local (CXmax, nprowsmax, 1))
           end if

           allocate (FX_FsX_compact_local (FXmax, nprowsmax, 2, 1))
           allocate (HX_SX_compact_local (HXmax, ncrowsmax, 2, 1))

! Perform a line search.

           prev_step = 0.0d0
           prev_ebs = ebs_local
           next_step = -0.001d0 !min_stepsize
           larger_step = .true.
           do jmu = 1,32

              multA = 1.0d0
              multB = next_step
              call sparse_add (nprows, ncmax, nGmax, ncmax, nprowsmax, nprowsmax, &
                   &                    nprowsmax, numc_local, listc_local,              &
                   &                    c_compact_local, multA, ipstart,                  &
                   &                    numRG_local,                          &
                   &                    listRG_local,                       &
                   &                    RG_compact_local, multB, ipstart,   &
                   &                    numcp_local, listcp_local, cp_compact_local)

              call sparse_add (ncrows, nctmax, nGtmax, nctmax, ncrowsmax, ncrowsmax, &
                   &                    ncrowsmax, numct, listct,              &
                   &                    ct_compact, multA, icstart,                  &
                   &                    numRGt_local,                          &
                   &                    listRGt_local,                       &
                   &                    RGt_compact_local, multB, icstart,   &
                   &                    numcpt_local, listcpt_local, cpt_compact_local)

              if (xuse) then
                 chooseA(1) = 1
                 chooseB(1) = 1

                 call sendrecv (myrank, 1, nbdiv, nbmod, CX_bit_matrix,       &
                      &                 iteration.eq.0,                                       &
                      &                 numcp_local, listcp_local, cp_compact_local, nprows,  &
                      &                 ncmax, nprowsmax, ipstart, numX_local,              &
                      &                 listX_local, X_compact_local, ncrows, Xmax,&
                      &                 ncrowsmax, icstart, numCX, listCX,      &
                      &                 CXp_compact_local, chooseA, chooseB, 1,     &
                      &                 CXmax, nprowsmax)
              end if

              chooseA(1) = 1
              chooseA(2) = 2
              chooseB(1) = 1
              chooseB(2) = 1

              if (xuse) then
                 call sendrecv (myrank, 1, nodiv, nomod, hs_bit_matrix, .false., &
                      & numh, listh, h_s_compact, nprows, nhmax, nprowsmax, ipstart, &
                      & numCX, listCX, CXp_compact_local, nprows, CXmax, nprowsmax, ipstart, &
                      & numFX, listFX, FX_FsX_compact_local, chooseA, chooseB, 2, FXmax, nprowsmax)
              else
                 call sendrecv (myrank, 1, nodiv, nomod, hs_bit_matrix, .false., &
                      & numh, listh, h_s_compact, nprows, nhmax, nprowsmax, ipstart, &
                      & numcp_local, listcp_local, cp_compact_local, nprows, ncmax, nprowsmax, ipstart, &
                      & numFX, listFX, FX_FsX_compact_local, chooseA, chooseB, 2, FXmax, nprowsmax)
              end if
              chooseA(1) = 1
              chooseA(2) = 1
              chooseB(1) = 1
              chooseB(2) = 2

              call sendrecv (myrank, 1, nodiv, nomod, HX_SX_bit_matrix, iteration.eq.0, &
                   & numcpt_local, listcpt_local, cpt_compact_local, ncrows, nctmax, ncrowsmax, icstart, &
                   & numFX, listFX, FX_FsX_compact_local, nprows, FXmax, nprowsmax, ipstart, &
                   & numHX, listHX, HX_SX_compact_local, chooseA, chooseB, 2, HXmax, &
                   & ncrowsmax)

! Add positive and negative values separately for more stability.
              indexH = 1
              indexS = 2
              ebs1 = 0.0d0
              ebs2 = 0.0d0
              eta_f1 = N
              eta_f2 = 0.0d0
!$omp parallel do private(index) reduction(+:ebs1,ebs2,eta_f1,eta_f2)
              do imu = 1, ncrows
                 do inu = 1, numHX(imu)
                    index = listHX(inu,imu)
                    if (index .eq. imu + icstart - 1) then
                       if (HX_SX_compact_local(inu,imu,indexH,1) .ge. 0) then
                          ebs1 = ebs1 + 4.0d0*HX_SX_compact_local(inu,imu,indexH,1)
                       else
                          ebs2 = ebs2 + 4.0d0*HX_SX_compact_local(inu,imu,indexH,1)
                       end if
                       if (HX_SX_compact_local(inu,imu,indexS,1) .le. 0) then
                          eta_f1 = eta_f1 + (-4.0d0)*HX_SX_compact_local(inu,imu,indexS,1)
                       else
                          eta_f2 = eta_f2 + (-4.0d0)*HX_SX_compact_local(inu,imu,indexS,1)
                       end if
                    end if
                    if (HX_SX_compact_local(inu,imu,indexH,1)*HX_SX_compact_local(inu,imu,indexS,1) .le. 0) then
                       ebs1 = ebs1                                                &
                            &     + (-2.0d0)*HX_SX_compact_local(inu,imu,indexH,1)*HX_SX_compact_local(inu,imu,indexS,1)
                    else
                       ebs2 = ebs2                                                &
                            &     + (-2.0d0)*HX_SX_compact_local(inu,imu,indexH,1)*HX_SX_compact_local(inu,imu,indexS,1)
                    end if
                    eta_f1 = eta_f1                                            &
                         &     + 2.0d0*HX_SX_compact_local(inu,imu,indexS,1)*HX_SX_compact_local(inu,imu,indexS,1)
                 end do
              end do
              eta_factor = eta_f1 + eta_f2
              ebs_temp(1) = ebs1 + ebs2 + eta * eta_factor

              call MPI_ALLREDUCE (ebs_temp, ebs, 1, mpi_whatever_double,       &
                   &                      MPI_SUM, MPI_BTN_WORLD, ierror)

              if (ebs(1) .le. prev_ebs) then
                 prev_ebs = ebs(1)
                 prev_step = next_step
                 if (larger_step) then
                    next_step = next_step * 2.0d0
                 else
                    next_step = next_step / 2.0d0
                 end if
              else if (prev_step.eq.0.0d0) then
                 prev_ebs = ebs(1)
                 prev_step = next_step
                 larger_step = .false.
                 next_step = next_step / 2.0d0
              else
                 stepsize = prev_step
                 exit
              end if
              stepsize = next_step
           end do
        end if
        deallocate (FX_FsX_compact_local, HX_SX_compact_local)
        deallocate (numcp_local, listcp_local, cp_compact_local)
        deallocate (numcpt_local, listcpt_local, cpt_compact_local)
        if (xuse) then
           deallocate (CXp_compact_local)
        end if
        deallocate (numRGt_local, listRGt_local, RGt_compact_local)

      end subroutine getstepsize


      SUBROUTINE NGAUSS(N,A,IA,B,X)   
      DOUBLE PRECISION A,B,X
      DIMENSION A(IA,N),B(N),X(N)     
      DO 4 K = 1,N-1
        DO 3 I = K+1,N      
          XMULT = A(I,K)/A(K,K)       
          DO 2 J = K+1,N    
            A(I,J) = A(I,J) - XMULT*A(K,J)      
   2      CONTINUE
          A(I,K) = XMULT
          B(I) = B(I) - XMULT*B(K)    
   3    CONTINUE  
   4  CONTINUE    
      X(N) = B(N)/A(N,N)
      DO 6 I = N-1,1,-1     
        SUM = B(I)
        DO 5 J = I+1,N      
          SUM = SUM - A(I,J)*X(J) 
   5    CONTINUE  
        X(I) = SUM/A(I,I)   
   6  CONTINUE    
      RETURN
      END 
