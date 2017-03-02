! copyright info:
!
! @Copyright 2005
! Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
!
! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio State University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! bfgs.f90
! Program Description
! ===========================================================================
!
! FUNCTION bfgs
! it call l-bfgs-b routine till
!	1)the convergence criteria defined in cg.optional are reached (f(i,j)<fmax for all i,j) - istatus=1
!   2)the algorithm is not able to find the smaller energy in spite of f>fmax (common in case the forces are not precise) - istatus=1
!	3)the maximum number of steps exceeds (istatus = 3)
!   4)an error in l-bfgs-b routine occure - program is stopped with an error (it could never happen)
!
! output:
!
! ===========================================================================
! Code written by:
! Z.Chromcova
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! email: chrom@fzu.cz
!
! ===========================================================================
! functions and subroutine declaration

!
! Program Declaration
! ===========================================================================

Subroutine bfgs(istatus)
    use configuration, only: natoms, ratom
    use energy, only: etot
    use forces, only: ftot
    use optimization, only: cg_maxstep, cg_minint, force_tol, freeParamsCount, mask, icg2md
    Implicit None
    !return value
    integer, intent(inout) :: istatus
    !local parameters:
    !parameters used in routine l-bfgs-b
    character(60) :: task, csave
    logical :: lsave(4)
    integer :: n, m, iprint, nbd(freeParamsCount), iwa(3 * freeParamsCount)
    integer :: isave(44)
    double precision :: f, factr, pgtol
    double precision :: dsave(29)
    double precision, allocatable, save :: wa(:)
    double precision :: l(freeParamsCount)
    double precision :: u(freeParamsCount)
    double precision :: grad(freeParamsCount)
    double precision :: xact(freeParamsCount)

    real :: ratomOld(3, natoms), forceOld(3, natoms)
    double precision :: Eold, gabs2
    integer :: iter, ind, nuphill
    integer :: withForce
    integer :: MAX_ITER

    !the local function delclaration
    logical, external :: minimumReached
    external :: writeOutRFE, getEnergyAndForces, ratom2xact, xact2ratom, ftot2grad, writeAnswerBas
    integer :: i
    integer :: nwrong = 0
    withForce = 0
    MAX_ITER = cg_maxstep
    
    iter = 0
    iprint = -1
    ! We suppress both code-supplied stopping tests because the
    ! user is providing his own stopping criteria. pgtol = maximum force on one component of g
    factr = 0.0
    pgtol = 0.0
    ! We now specify nbd which defines the bounds on the variables:
    !                l   specifies the lower bounds,
    !                u   specifies the upper bounds.
    !     nbd(i)=0 if x(i) is unbounded,
    !            1 if x(i) has only a lower bound,
    !            2 if x(i) has both lower and upper bounds, and
    !            3 if x(i) has only an upper bound.

    nbd(1:freeParamsCount) = 0
    l(1:freeParamsCount) = 0
    u(1:freeParamsCount) = 0

    n = freeParamsCount
    m = 7 ! m is the number of limited memory corrections 3<=m<=20
    ! inicialization
    if (.not.allocated(wa)) allocate(wa(2 * m * n + 4 * n + 12 * m * m + 12 * m))
    iwa = 0
    wa = 0.d0
    isave = 0
    dsave = 0.d0
    lsave = .false.

    call ratom2xact(ratom, xact, natoms)
    call getEnergyAndForces(gabs2, iter)
    call writeOutRFE(ratom, ftot, etot, 'coordinates and forces        ', iter)
    ratomOld = ratom
    forceOld = ftot
    Eold = etot
    gabs2 = 0
    !     We start the iteration by initializing task.
    task = 'START'
    nuphill = 0
    
    !------- the beginning of the loop ----------
    do iter = 1, MAX_ITER 
        if (withForce /= 1)then
            f = etot
        else
            !this part is used in case that Etot does not change more but the forces are still big
            f = gabs2 !		f=sqrt(gabs2) ! it find better minima but it is quite time consuming
        endif

        call ftot2grad(ftot, grad, natoms)
        call setulb(n, m, xact, l, u, nbd, f, grad, factr, pgtol, wa, iwa, task, iprint, &
        csave, lsave, isave, dsave)
        call xact2ratom(ratom, xact, natoms)

        if (task(1:2) == 'FG') then
            !the routine ask for the next function value and gradient
            call getEnergyAndForces(gabs2, iter)
        elseif (task(1:5) == 'NEW_X' .or. task(1:5) == 'START') then
            !the minimization routine has returned with a new iterate, or the minimization was restarted
            call getEnergyAndForces(gabs2, iter)

            if (etot <= Eold)then
                call writeOutRFE(ratom, ftot, etot, 'coordinates and forces        ', iter)
                call writeAnswerBas(ratom, natoms)
                if (minimumReached(ratom, ftot, natoms, etot, Eold, iter))then
                    istatus = 1
                    write(*,*)'++ The user-defined criteria of convergence filled. Fmax < Ftol'
                    deallocate(wa)
                    return
                endif
                Eold = etot
                ratomOld = ratom
                forceOld = ftot
                !				nwrong = 0
            elseif (.false. .and. withForce == 1) then
                !the minimization according forces only (for case the forces are not precise)
                write(*, *) '++ BFGS: FORCES DECREASE, BUT ENERGY DOES NOT CHANGE MORE. Restating...'
                task = 'STOP'
                call setulb(n, m, xact, l, u, nbd, f, grad, factr, pgtol, wa, iwa, task, iprint, &
                csave, lsave, isave, dsave)
                task = 'START'
                withForce = 2
                nwrong = 0
            elseif(withForce == 1 )then
                istatus = 2
                deallocate(wa)
                write(*,*)'++ BFGS MINIMIZATION FAILED - energy increases'
!go back to the best position
                ratom = ratomOld
                ftot = forceOld
                etot = Eold
                return
            endif
        elseif (task(1:5) == 'ERROR')then
            write(*, *) task
            stop
        elseif (task(1:11) == 'CONVERGENCE')then
            write(*, *) '++ BFGS: ', task
            task = 'STOP'
            call setulb(n, m, xact, l, u, nbd, f, grad, factr, pgtol, wa, iwa, task, iprint, &
            csave, lsave, isave, dsave)
            nwrong = nwrong + 1
            if (withForce == 2)then
                istatus = 2
                if (etot > Eold)then
                    ratom = ratomOld
                    ftot = forceOld
                    etot = Eold
                endif
                write(*, *) '++ MINIMIZATIN STOPPED in spite of the user defined critera is not full-filled'
                deallocate(wa)
                return
            elseif ((nwrong >= cg_minint) .and. withForce == 0) then
!NOT USED NOW (see l. this+15)
                !restarting the minimization with function f(ratom)=abs(f)+etot (the Etot does not change more, now we try to minimize forces only)
!                withForce = 1
!                write(*, *) '++ RESTARTING THE MINIMIZATION WITH FUNCTI0N (ftot*ftot) or abs(ftot)'
!            elseif
                if(icg2md>0)then
                    write(*, *) '++ MINIMIZATION FAILS. PROCEEDS WITH ',icg2md,' STEPS OF FIRE MINIMIZATION'
                    call FIRE_loop_ex(icg2md,iter)
                    write(*, *) '++ RESTARTING THE MINIMIZATION WITH FUNCTI0N (etot) AFTER 100 FIRE STEPS'
                else
                    write(*, *) '++ RESTARTING THE MINIMIZATION WITH FUNCTI0N (etot) LASTTIME'
                    ratom = ratomOld
                    ftot = forceOld
                    etot = Eold
                endif
                withForce = 2
            else
                ratom = ratomOld
                ftot = forceOld
                etot = Eold
            endif
            task = 'START'
        elseif (task(1:30) == 'ABNORMAL_TERMINATION_IN_LNSRCH') then
            task = 'STOP'
            call setulb(n, m, xact, l, u, nbd, f, grad, factr, pgtol, wa, iwa, task, iprint, &
            csave, lsave, isave, dsave)
            write(*, *) '++ MINIMIZATION STOPPED - ABNORMAL TERMINATION. The user defined critera is not full-filled'
            deallocate(wa)
            return
        else
            write(*, *) '++ ERROR: Unknown task in minimization routine BFGS'
            write(*, *) task
            stop
        endif
    enddo
    istatus = 3
    ratom = ratomOld
    ftot = forceOld
    etot = Eold
    deallocate(wa)
end subroutine bfgs

subroutine ftot2grad(ftot, grad, natoms)
    use optimization, only:mask
    use optimization, only: mask, freeParamsCount
    implicit none
    integer :: natoms
    real :: ftot(3, natoms)
    real :: grad(freeParamsCount)
    integer :: ind, j, i
    i = 1

    do ind = 1, natoms
        do j = 1, 3
            if (mask(j, ind) == 1) then
                grad(i) = -ftot(j, ind) !nN to eV/nm
                i = i + 1
            endif
        enddo
    enddo
end subroutine

subroutine writeAnswerBas(ratom, natoms)
    use charges, only:nzx
    use interactions, only: imass
    use configuration, only: shifter, ishiftO
    use options, only: iclassicMD
    implicit none
    integer :: iatom
    integer, intent(in) :: natoms
    real, intent(in) :: ratom(3, natoms)

    if (iclassicMD == 0)then
        open (unit = 86, file = 'answer.bas', status = 'unknown')
        write (86, *) natoms
        do iatom = 1, natoms
            write(*, 701) nzx(imass(iatom)), ratom(:, iatom)
            if (ishiftO == 1) then
                write (86, 700) nzx(imass(iatom)), ratom(:, iatom) - shifter(:)
            else
                write (86, 700) nzx(imass(iatom)), ratom(:, iatom)
            endif
        enddo
        close(unit = 86)
    endif
    700 format (2x, i2, 3(2x, f18.8))
    701 format ('++ ', 2x, i2, 3(2x, f18.8))
end subroutine

subroutine ratom2xact(ratom, xact, natoms)
    use optimization, only: mask, freeParamsCount
    implicit none
    integer :: natoms
    real :: ratom(3, natoms)
    real :: xact(freeParamsCount)
    integer :: ind, j, i

    i = 1
    do ind = 1, natoms
        do j = 1, 3
            if (mask(j, ind) == 1) then
                xact(i) = ratom(j, ind)
                i = i + 1
            endif
        enddo
    enddo
end subroutine

subroutine xact2ratom(ratom, xact, natoms)
    use optimization, only: mask, freeParamsCount
    implicit none
    integer :: natoms
    real :: ratom(3, natoms)
    real :: xact(freeParamsCount)
    integer :: ind, j, i
    i = 1

    do ind = 1, natoms
        do j = 1, 3
            if (mask(j, ind) == 1) then
                ratom(j, ind) = xact(i)
                i = i + 1
            endif
        enddo
    enddo
end subroutine

subroutine getEnergyAndForces(fabs, iter)
    use options, only: iclassicMD, icluster, iimage, iquench
    use forces, only: ftot
    use configuration, only:natoms
    implicit none
    real, intent(inout) :: fabs
    integer, intent(in) :: iter
    integer :: i

    if (iimage .ge. 1) call imaged(icluster, iimage, iter, 1)

    ! ***************************************************************************
    ! junkermeier: this little section will incrementally increase/decrease
    !              the temperature over the course of a calculation.
    !THIS PART IS IMPRACTICAL, iquench=-5
    !         if (iendtemp == 1 .and. iquench == 0) then
    !            T_wantPrev = T_want
    !            T_want = T_initial + T_increment*itime_step
    !            call resetNHC(natoms,T_want,T_wantPrev)
    !         end if
    ! ***************************************************************************

    if (iclassicMD == 0)then
        ! perform SCF LOOP
        call scf_loop(iter)
        ! optionally perform post-processing (DOS etc.)
        call postscf()
        ! calculate the total energy
        call getenergy(iter)
    end if
    call getforces()
    fabs = 0
    do i = 1, natoms
        fabs = fabs + ftot(1, i) * ftot(1, i) + ftot(2, i) * ftot(2, i) + ftot(3, i) * ftot(3, i)
    enddo
end subroutine getEnergyAndForces

function minimumReached(ratom, force, natoms, energy, energyOld, iter)
    use optimization, only: mask, force_tol, energy_tol

    implicit none
    real, intent(in) :: energy, energyOld
    real, dimension(3, natoms) :: force, ratom
    integer, intent(in) :: natoms, iter

    interface
        function getMaxForce(idx, force, natoms, maxForce)
            real :: force(3, natoms), maxForce, getMaxForce
            integer :: natoms, idx
        end function
    end interface

    logical :: minimumReached
    integer :: i, j
    real :: maxForce
    character(30) :: msg, emptymsg = '                              '

    minimumReached = .false.
    maxForce = 0
    do i = 1, natoms
        maxForce = getMaxForce(i, force, natoms, maxForce)
    enddo

    if (maxForce < force_tol)then
        minimumReached = .true.
    elseif (abs(energy - energyOld) < energy_tol)then
        minimumReached = .false.
    endif

    if (minimumReached)then
        write(*, *) '============================================================'
        write(*, *) '= MINIMUM REACHED ROUTINE RETURNS TRUE: maxForce<force_tol ='
        write(*, *) '=  maxForce=', maxForce, 'force_tol=', force_tol
        write(*, *) '============================================================'
        msg = emptymsg
        msg = 'ratom forces:'
        call writeOutRFE(ratom, force, energy, msg, iter)
    endif
end function

function getMaxForce(idx, force, natoms, maxForce)
    use optimization, only: mask
    implicit none
    real :: force(3, natoms), maxForce, fProj, getMaxForce
    integer :: i, natoms, idx

    fProj = mask(1, idx) * force(1, idx) * force(1, idx) + &
    mask(2, idx) * force(2, idx) * force(2, idx) + &
    mask(3, idx) * force(3, idx) * force(3, idx)

    if (fProj > maxForce * maxForce)then
        getMaxForce = sqrt(fProj)
    else
        getMaxForce = maxForce
    endif

    return
end function

subroutine writeOutRFE(ratom, force, val, msg, cg_iter)
    use configuration, only: natoms, symbol, shifter, ishiftO, basisfile
    use interactions, only: wrtout, imass, nssh
    use charges, only:nzx
    use outputs, only: iwrtxyz
    use options, only:iqout, iquench
    use charges, only:Qin
    implicit none

    interface
        function getMaxForce(idx, force, natoms, maxForce)
            real :: force(3, natoms), maxForce, getMaxForce
            integer :: natoms, idx
        end function
    end interface

    real, intent(in) :: val
    real, dimension(3, natoms), intent(inout) :: force
    real, dimension(3, natoms), intent(inout) :: ratom
    character(30), intent(in) :: msg
    integer, intent(in) :: cg_iter
    real :: maxForce

    integer :: i, ix

    maxForce = 0
    if (wrtout) then
        write(*, *) '++ ========================================='
        write(*, *) '++', msg, '++'
        write(*, *) '++ Etot= ', val, '[eV]'
        do i = 1, natoms
            write(*, 36) i, ratom(:, i), force(:, i), sqrt(force(1, i)**2 + force(2, i)**2 + force(3, i)**2)
            maxForce = getMaxForce(i, force, natoms, maxForce)
        enddo
    else
        do i = 1, natoms
            maxForce = getMaxForce(i, force, natoms, maxForce)
        enddo
    endif

    if (iquench == -5)then
        write (*, 501) cg_iter, val, maxForce
    elseif (iquench == -6)then
        write (*, 503) cg_iter, val, maxForce
    else
        write (*, 505) cg_iter, val, maxForce
    endif

    ! writeout XYZ
    if (iwrtxyz .eq. 1)then
        open (unit = 87, file = 'answer.xyz', status = 'unknown', position = 'append')
        write (87, *) natoms
        if (iquench == -5)then
            write (87, 500) cg_iter, val, maxForce
        elseif (iquench == -6)then
            write (87, 502) cg_iter, val, maxForce
        else
            write (87, 504) cg_iter, val, maxForce
        endif

        do i = 1, natoms
            if (ishiftO .eq. 1) then
                write (87, 800) symbol(i), ratom(:, i) - shifter(:)
            else
                write (87, 800) symbol(i), ratom(:, i)
            endif
        end do
        close(unit = 87)
    endif

    ! writeout CHARGES
    open (unit = 87, file = 'CHARGES.min', status = 'unknown')
    write (87, 600) natoms, basisfile, iqout
    do i = 1, natoms
        write(87, 601) (Qin(ix, i), ix = 1, nssh(imass(i)))
    end do
    close(87)
    500 format ('## BFGS step no.', i6, ', E=', f20.8, ' fmax=', e18.4)
    501 format ('+++ BFGS step no.', i6, ', E=', f20.8, ' fmax=', e18.4)
    502 format ('## FIRE step no.', i6, ', E=', f20.8, ' fmax=', e18.4)
    503 format ('+++ FIRE step no.', i6, ', E=', f20.8, ' fmax=', e18.4)
    504 format ('## CG step no.', i6, ', E=', f20.8, ' fmax=', e18.4)
    505 format ('+++ CG step no.', i6, ', E=', f20.8, ' fmax=', e18.4)
    600 format (2x, i5, 2x, a40, 2x, i2)
    601 format (2x, 10f14.8)
    700 format (2x, i2, 3(2x, f18.8))
    800 format (2x, a2, 3(2x, f12.6))
    36 format(' ++ atom= ', i4, 3f12.5, ', f= ', 3e12.2, ', f=', e12.2)
end subroutine

subroutine writeOutRE(ratom, val, msg)
    use configuration, only: natoms
    use interactions, only: wrtout
    implicit none

    real, intent(in) :: val
    real, dimension(3, natoms), intent(inout) :: ratom
    character(30), intent(in) :: msg
    integer :: i

    if (wrtout) then
        write(*, *) '++ ========================================='
        write(*, *) msg
        write(*, *) '++ Etot= ', val, '[eV]'
        do i = 1, natoms
            write(*, 37) i, ratom(:, i)
        enddo
    endif
    37 format(' ++ atom= ', i4, 3e12.2)
end subroutine

