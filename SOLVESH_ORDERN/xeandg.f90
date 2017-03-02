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

! xeandg.f90
! Program Description
! ===========================================================================
!       This routine optimizes the matrix X over the functional 
! min Tr[(-S)(2X-XSX)] as in equation (20) in Yang 1997.
!
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
        subroutine xeandg (ncrows, icstart, iteration, &
             & nSmax, numSlocal, listSlocal, S_compactlocal, &
             & trace_local, gnorm, stop_now)
        use interactions
        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ncrows, icstart, iteration, nSmax
        integer, intent (in), dimension (ncrowsmax) :: numSlocal
        integer, intent (in), dimension (nSmax, ncrowsmax) :: listSlocal
        real, intent (in), dimension (nSmax, ncrowsmax) :: S_compactlocal

! Output
        real, intent (out) :: trace_local, gnorm
        logical, intent (out) :: stop_now
  
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ichoosesize
        integer ierror
        integer imu, inu, jmu, index
        integer nodiv, nomod, nbdiv, nbmod
        integer myrank
        integer, save :: SXmax, SSXmax, XGmax
        integer XGmax_local, SXmax_local, SSXmax_local
        integer, save :: reset_iter
        integer, dimension (2) :: ichooseA
        integer, dimension (2) :: ichooseB

        double precision trace1, trace2
        double precision, dimension (2) :: trace_temp, trace_temp2
        double precision c1, c2

        double precision, dimension (2) :: prod(2)
        double precision, dimension (2) :: prodt(2)

        real multA 
        real multB
        real, save :: stepsize, trace_old
        real gnorm_local

! Holds the X matrix for the next step.
        real, dimension (:, :), allocatable, save :: altX_compact_local

        logical need_list

! Bit matrix indicating which blocks of a product matrix are nonzero.
        logical, dimension (:, :), allocatable, save :: SX_bit_matrix, SSX_bit_matrix

! ****************************************************************************
! Temporary storage vectors. 
! ****************************************************************************
        integer, dimension (:), allocatable :: numTg1_local
        integer, dimension (:), allocatable :: numTg2_local
        integer, dimension (:, :), allocatable :: listTg1_local
        integer, dimension (:, :), allocatable :: listTg2_local

        real, dimension (:, :), allocatable, target :: Tg1_compact_local, Tg2_compact_local

! ****************************************************************************
! SX matrix declaration - local pieces. Defined by SX = S * X. 
! ****************************************************************************
        integer, dimension (:), allocatable, save :: numSX_local
        integer, dimension (:, :), allocatable, save :: listSX_local

        real, dimension (:, :), allocatable :: SX_compact_local

! ****************************************************************************
! SXp matrix declaration - local pieces.
! ****************************************************************************
        real, dimension (:, :, :), allocatable :: SXp_compact_local

! ****************************************************************************
! SSX matrix declaration - local pieces. Defined by SSX = S^2 * X. 
! ****************************************************************************
        integer, dimension (:), allocatable, save :: numSSX_local
        integer, dimension (:, :), allocatable, save :: listSSX_local

        real, dimension (:, :), allocatable :: SSX_compact_local

! ****************************************************************************
! Xp matrix declaration - local pieces. X +/- stepsize * G.
! ****************************************************************************
        integer, dimension (:), allocatable :: numXp_local
        integer, dimension (:, :), allocatable :: listXp_local

        real, dimension (:, :, :), allocatable :: Xp_compact_local

! ****************************************************************************
! XG matrix declaration - local pieces.
! ****************************************************************************
        integer, dimension (:), allocatable, save :: numXG_local
        integer, dimension (:, :), allocatable, save :: listXG_local

        real, dimension (:, :), allocatable, save :: XG_compact_local

! ****************************************************************************
! R matrix declaration - local pieces.
! ****************************************************************************
        integer, dimension (:), allocatable, save :: numR_local
        integer, dimension (:, :), allocatable, save :: listR_local

        real, dimension (:, :), allocatable, save :: R_compact_local

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ============================================================================
 
! Procedure
! ============================================================================
        stop_now = .false.

! Initialize the variable need_list
! On the first iteration this is set to true, so that the index sets of the 
! various matrices may be constructed.  Once constructed they do not change 
! (theoretically).
        need_list = .false.
        if (iteration .eq. 0) then
         need_list = .true.
         if (allocated(altX_compact_local)) deallocate(altX_compact_local)
         allocate (altX_compact_local (Xmax, ncrowsmax))
        else
! Here we set X to its new value computed during the last call.
! FIXME inefficient copy
         X_compact_local = altX_compact_local
        end if
! If this is the first iteration, reset the CG accumulator.
        if (iteration .eq. 0) reset_iter = 0

! Find out which processor this is.
        call MPI_COMM_RANK (MPI_BTN_WORLD, myrank, ierror)

! sizes of local sections
        nodiv = norbitals/nactualprocs
        nomod = mod(norbitals,nactualprocs)
        nbdiv = nbands/nactualprocs
        nbmod = mod(nbands,nactualprocs)

! ****************************************************************************
!                            SX   M A T R I X
! ****************************************************************************
! Allocate appropriate arrays.
!       write (*,*) ' Building SX matrix, myrank = ', myrank
        if (iteration .eq. 0) then
         call set_maxdimension (nactualprocs, myrank, ncrows, nSmax, Xmax, &
     &                          ncrowsmax, ncrows, ncrowsmax, numSlocal,   &
     &                          listSlocal, icstart, numX_local,           &
     &                          listX_local, icstart, 1, SXmax_local)
         call MPI_ALLREDUCE (SXmax_local, SXmax, 1, MPI_INTEGER, MPI_MAX,&
     &                       MPI_BTN_WORLD, ierror)
         if (allocated (numSX_local)) deallocate (numSX_local, listSX_local)
         allocate (numSX_local (ncrowsmax))
         allocate (listSX_local (SXmax, ncrowsmax))
         if (allocated(SX_bit_matrix)) deallocate(SX_bit_matrix)
         allocate (SX_bit_matrix (nactualprocs, nactualprocs))
        end if
        allocate (SX_compact_local (SXmax, ncrowsmax))

! Here we multiply S by X to get SX.
        ichoosesize = 1
        ichooseA(1) = 1
        ichooseB(1) = 1
        call sendrecv (myrank, iteration, nbdiv, nbmod, SX_bit_matrix,  &
     &                 iteration.eq.0, numSlocal, listSlocal, S_compactlocal,  &
     &                 ncrows, nSmax, ncrowsmax, icstart, numX_local,       &
     &                 listX_local, X_compact_local, ncrows, Xmax,       &
     &                 ncrowsmax, icstart, numSX_local, listSX_local,      &
     &                 SX_compact_local, ichooseA, ichooseB,            &
     &                 ichoosesize, SXmax, ncrowsmax)

! ****************************************************************************
!                            SSX   M A T R I X
! ****************************************************************************
! Allocate appropriate arrays.
!       write (*,*) ' Building SSX matrix, myrank = ', myrank
        if (iteration .eq. 0) then
         call set_maxdimension (nactualprocs, myrank, ncrows, nSmax, SXmax, &
     &                          ncrowsmax, ncrows, ncrowsmax, numSlocal,   &
     &                          listSlocal, icstart, numSX_local,           &
     &                          listSX_local, icstart, 1, SSXmax_local)

         call MPI_ALLREDUCE (SSXmax_local, SSXmax, 1, MPI_INTEGER, MPI_MAX,&
     &                       MPI_BTN_WORLD, ierror)
         if (allocated (numSSX_local)) deallocate (numSSX_local, listSSX_local)
         allocate (numSSX_local (ncrowsmax))
         allocate (listSSX_local (SSXmax, ncrowsmax))
         if (allocated(SSX_bit_matrix)) deallocate(SSX_bit_matrix)
         allocate (SSX_bit_matrix (nactualprocs, nactualprocs))
        end if
        allocate (SSX_compact_local (SSXmax, ncrowsmax))

! Here we multiply S by SX to get SSX.
        ichoosesize = 1
        ichooseA(1) = 1
        ichooseB(1) = 1
        call sendrecv (myrank, iteration, nbdiv, nbmod, SSX_bit_matrix,  &
     &                 iteration.eq.0, numSlocal, listSlocal, S_compactlocal,  &
     &                 ncrows, nSmax, ncrowsmax, icstart, numSX_local,       &
     &                 listSX_local, SX_compact_local, ncrows, SXmax,       &
     &                 ncrowsmax, icstart, numSSX_local, listSSX_local,      &
     &                 SSX_compact_local, ichooseA, ichooseB,            &
     &                 ichoosesize, SSXmax, ncrowsmax)

! Compute the trace value at the current point, using the functional 
! formula. We add positive and negative values separately for greater stability.
        trace1 = 0.0d0
        trace2 = 0.0d0
!$omp parallel do private(index) reduction(+:trace1,trace2)
        do imu = 1, ncrows
         do inu = 1, numSX_local(imu)
          index = listSX_local(inu,imu)
          if (index .eq. imu + icstart - 1) then
           if (SX_compact_local(inu,imu) .le. 0) then
            trace1 = trace1 - 2.0d0*SX_compact_local(inu,imu)
           else
            trace2 = trace2 - 2.0d0*SX_compact_local(inu,imu)
           end if
          end if
          trace1 = trace1 + SX_compact_local(inu,imu) * SX_compact_local(inu,imu)
         end do
        end do
        trace_temp(1) = trace1 + trace2

! Compute trace value at this point.
        call MPI_ALLREDUCE (trace_temp, trace_temp2, 1,              &
     &                      mpi_whatever_double, MPI_SUM, MPI_BTN_WORLD,    &
     &                      ierror)

        trace_local = trace_temp2(1)

! Reset CG (take a SD step) if the trace value has increased.
        if (reset_iter .gt. 0 .and. trace_local .ge. trace_old) then
         if (myrank.eq.0) write (*,*)                                        &
     &    ' Warning: energy has not decreased, taking a steepest descent step! '
         reset_iter = 0
        end if

! ****************************************************************************
!                     C A L C U L A T E   G R A D I E N T S
! ****************************************************************************
! First Determine the maximum value for the parameter XGmax.
! Allocate appropriate arrays.
!       write (*,*) ' Building gradient matrices, myrank = ', myrank
        if (iteration .eq. 0) then
         XGmax = nSmax + SSXmax
! Allocate the gradient matrix
         if (allocated (numXG_local))                                         &
     &    deallocate (numXG_local, listXG_local, XG_compact_local)
         allocate (numXG_local (ncrowsmax))
         allocate (listXG_local (XGmax, ncrowsmax))
         allocate (XG_compact_local (XGmax, ncrowsmax))
! Allocate the CG step direction, which accumulates from one iteration to 
! the next
         if (allocated (numR_local))                                        &
     &    deallocate (numR_local, listR_local, R_compact_local)
         allocate (numR_local (ncrowsmax))
         allocate (listR_local (XGmax, ncrowsmax))
         allocate (R_compact_local (XGmax, ncrowsmax))
        end if

        allocate (numTg1_local (ncrowsmax))
        allocate (listTg1_local (XGmax, ncrowsmax))
        allocate (Tg1_compact_local (XGmax, ncrowsmax))
        allocate (numTg2_local (ncrowsmax))
        allocate (listTg2_local (XGmax, ncrowsmax))
        allocate (Tg2_compact_local (XGmax, ncrowsmax))

! Add terms of the gradient expression
        multA = -2.0d0
        multB = 2.0d0
        call sparse_add (ncrows, nSmax, SSXmax, XGmax, ncrowsmax, ncrowsmax,  &
     &                   ncrowsmax, numSlocal, listSlocal,             &
     &                   S_compactlocal, multA, icstart, numSSX_local,    &
     &                   listSSX_local, SSX_compact_local, multB, icstart,   &
     &                   numTg1_local, listTg1_local, Tg1_compact_local)

! constrain the gradient by the X matrix index set
        call sparse_mask (ncrows, XGmax, Xmax, XGmax, ncrowsmax, ncrowsmax,  &
     &                   ncrowsmax, numTg1_local, listTg1_local,             &
     &                   Tg1_compact_local, icstart, numX_local,    &
     &                   listX_local, icstart,   &
     &                   numTg2_local, listTg2_local, Tg2_compact_local)

        if (reset_iter .gt. 0) then

! Saves the difference between the old and new gradients.
         multA = 1.0d0
         multB = -1.0d0
         call sparse_add (ncrows, XGmax, XGmax, XGmax, ncrowsmax, ncrowsmax, &
     &                    ncrowsmax, numTg2_local, listTg2_local,            &
     &                    Tg2_compact_local, multA, icstart, numXG_local,     &
     &                    listXG_local, XG_compact_local, multB, icstart,      &
     &                    numTg1_local, listTg1_local, Tg1_compact_local)

         call sparse_norm2 (ncrows, XGmax, ncrowsmax, numXG_local,            &
              &                      listXG_local, XG_compact_local, 1, .false., prodt(2))

        end if

! Saves the new gradient for use by the next iteration.
        numXG_local = numTg2_local
        listXG_local = listTg2_local
        XG_compact_local = Tg2_compact_local

        if (reset_iter .gt. 0) then
! Conjugate gradient (Polak-Ribiere; see Nocedal/Wright, Numerical
! Optimization, pp. 120-122)
         call sparse_vecprod (ncrows, XGmax, XGmax, ncrowsmax, ncrowsmax,    &
     &                        numTg1_local, listTg1_local, Tg1_compact_local,&
     &                        1, numTg2_local, listTg2_local,                &
     &                        Tg2_compact_local, 1, .false., prodt(1))

         call MPI_ALLREDUCE (prodt, prod, 2, mpi_whatever_double, MPI_SUM,  &
     &                       MPI_BTN_WORLD, ierror)

! FIXME make sure elsewhere that XG is nonzero so the division will be safe
         multA = 1.0d0
         multB = prod(1) / prod(2)
         call sparse_add (ncrows, XGmax, XGmax, XGmax, ncrowsmax, ncrowsmax, &
     &                    ncrowsmax, numTg2_local, listTg2_local,            &
     &                    Tg2_compact_local, multA, icstart, numR_local,    &
     &                    listR_local, R_compact_local, multB, icstart,    &
     &                    numTg1_local, listTg1_local, Tg1_compact_local)

! FIXME inefficient copy
         multA = 1.0d0
         call sparse_copy (XGmax, XGmax, ncrows, ncrowsmax, ncrows,          &
     &                     ncrowsmax, numTg1_local, listTg1_local,           &
     &                     Tg1_compact_local, icstart, icstart, 0, .false.,  &
     &                     .false., multA, numR_local, listR_local,        &
     &                     R_compact_local)
        else

! Initialize the CG search direction with the current gradient.
! FIXME inefficient copy
         numR_local = numTg2_local
         listR_local = listTg2_local
         R_compact_local = Tg2_compact_local
        end if

! Calculate the RMS of the gradient
        gnorm_local = 0.0d0
!$omp parallel do reduction(+:gnorm_local)
        do imu = 1, ncrows
         do inu = 1, numXG_local(imu)
          gnorm_local = gnorm_local + XG_compact_local(inu,imu)**2
         end do
        end do

        call MPI_ALLREDUCE (gnorm_local, gnorm, 1, mpi_whatever_real,        &
     &                      MPI_SUM, MPI_BTN_WORLD, ierror)
        gnorm = sqrt(gnorm / (nbands*nbands))

        deallocate (numTg1_local, listTg1_local, Tg1_compact_local)
        deallocate (numTg2_local, listTg2_local, Tg2_compact_local)
        deallocate (SX_compact_local, SSX_compact_local)

! ****************************************************************************
!                   D E T E R M I N E   S T E P   S I Z E
! ****************************************************************************

! Computes the optimal step length by computing the value at two nearby points
! and interpolating.
        if (iteration.eq.0) then
           stepsize = -min_stepsize
        else
           stepsize = abs(stepsize)
        end if

        allocate (numXp_local(ncrowsmax))
        allocate (listXp_local(Xmax, ncrowsmax))
        allocate (Xp_compact_local(Xmax, ncrowsmax, 2))

        multA = 1.0d0
        multB = stepsize
        call sparse_add (ncrows, Xmax, XGmax, Xmax, ncrowsmax, ncrowsmax, &
     &                    ncrowsmax, numX_local, listX_local,                &
     &                    X_compact_local, multA, icstart,                   &
     &                    numR_local, listR_local,                           &
     &                    R_compact_local, multB, icstart,   &
     &                    numXp_local, listXp_local, Xp_compact_local(:,:,1))

        multA = 1.0d0
        multB = -stepsize
        call sparse_add (ncrows, Xmax, XGmax, Xmax, ncrowsmax, ncrowsmax, &
     &                    ncrowsmax, numX_local, listX_local,                &
     &                    X_compact_local, multA, icstart,                   &
     &                    numR_local, listR_local,                           &
     &                    R_compact_local, multB, icstart,   &
     &                    numXp_local, listXp_local, Xp_compact_local(:,:,2))

        allocate (SXp_compact_local (SXmax, ncrowsmax, 2))

! Here we multiply S by Xp to get SXp.
        ichoosesize = 2
        ichooseA(1) = 1
        ichooseB(1) = 1
        ichooseA(2) = 1
        ichooseB(2) = 2
        call sendrecv (myrank, 1, nbdiv, nbmod, SX_bit_matrix,  &
     &                 iteration.eq.0, numSlocal, listSlocal, S_compactlocal,  &
     &                 ncrows, nSmax, ncrowsmax, icstart, numXp_local,       &
     &                 listXp_local, Xp_compact_local, ncrows, Xmax,       &
     &                 ncrowsmax, icstart, numSX_local, listSX_local,      &
     &                 SXp_compact_local, ichooseA, ichooseB,            &
     &                 ichoosesize, SXmax, ncrowsmax)

        do jmu = 1, 2
           trace1 = 0.0d0
           trace2 = 0.0d0
!$omp parallel do private(index) reduction(+:trace1,trace2)
           do imu = 1, ncrows
              do inu = 1, numSX_local(imu)
                 index = listSX_local(inu,imu)
                 if (index .eq. imu + icstart - 1) then
                    if (SXp_compact_local(inu,imu,jmu) .le. 0) then
                       trace1 = trace1 - 2.0d0*SXp_compact_local(inu,imu,jmu)
                    else
                       trace2 = trace2 - 2.0d0*SXp_compact_local(inu,imu,jmu)
                    end if
                 end if
                 trace1 = trace1 + SXp_compact_local(inu,imu,jmu) * SXp_compact_local(inu,imu,jmu)
              end do
           end do
           trace_temp(jmu) = trace1 + trace2
        end do

! Compute trace value at this point.
        call MPI_ALLREDUCE (trace_temp, trace_temp2, 2,              &
     &                      mpi_whatever_double, MPI_SUM, MPI_BTN_WORLD,    &
     &                      ierror)

!        if (myrank.eq.0) write (*,*) 'trace_local=',trace_local,'trace1=',trace_temp2(1),'trace-1=',trace_temp2(2)
        c1 = (trace_temp2(1)-trace_temp2(2)) / (2*stepsize)
        c2 = (trace_temp2(2)+trace_temp2(1)-2*trace_local) / (2*stepsize*stepsize)
!        if (myrank.eq.0) write (*,*) 'c1=',c1,'c2=',c2
        stepsize = -c1/(2*c2)
! If the stepsize is positive make it negative.
! This seems to work better in practice.
! If the function is positive definite then the stepsize
! should, theoretically, always be negative, yielding a
! descent direction.
        stepsize = -abs(stepsize)
! Make sure abs(stepsize) is no less than the minimum, to ensure progress.
        if (stepsize .gt. min_stepsize) stepsize = min_stepsize
!        if (myrank.eq.0) write (*,*) 'stepsize=',stepsize
        deallocate (SXp_compact_local)
        deallocate (numXp_local, listXp_local, Xp_compact_local)

        allocate (numTg1_local (ncrowsmax))
        allocate (listTg1_local (Xmax, ncrowsmax))

! Saves the next X matrix.  The X matrix is not changed until the next iteration.
! We assume here that the gradient has the same or less sparsity as X.
        multA = 1.0d0
        multB = stepsize
        call sparse_add (ncrows, Xmax, XGmax, Xmax, ncrowsmax, ncrowsmax, &
     &                    ncrowsmax, numX_local, listX_local,                &
     &                    X_compact_local, multA, icstart,                   &
     &                    numR_local, listR_local,                           &
     &                    R_compact_local, multB, icstart,   &
     &                    numTg1_local, listTg1_local, altX_compact_local)

        deallocate (numTg1_local, listTg1_local)

! Save the trace value to compare it later.
        trace_old = trace_local

! Update the CG iteration index.
        reset_iter = reset_iter + 1


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
	 
        return
        end subroutine xeandg

