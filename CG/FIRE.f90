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


! FIRE.f90
! Program Description
! ===========================================================================
!
! minimization algorithm PRL97(2006)17021 (method of blind skier)
!
! output:
!
! ===========================================================================
! Code written by (it is ispirated by FIRE routine in LAMPS):
! Z.Chromcova
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! email: chrom@fzu.cz
! ===========================================================================
! functions and subroutine declaration

subroutine iterate_FIRE(ntimestep, restart)
    use configuration, only: natoms, ratom, vatom, xmass, ratom
    use interactions, only: imass
    use energy, only: etot
    use forces, only: ftot
    use MD,only: dt
    use optimization, only: cg_drmax, force_tol, mask
    implicit none
    integer :: ntimestep
    logical, intent(inout) :: restart
    ! CONSTANTS:
    integer, parameter :: DELAYSTEP = 5
    real, parameter :: DT_GROW = 1.5
    real, parameter :: DT_SHRINK = 0.5
    real, parameter :: ALPHA0 = 0.1
    real, parameter :: ALPHA_SHRINK = 0.99
    real, save :: TMAX = 50
    real, save :: TMIN = 1.0e-6
    real, save :: DMAX
    ! VARIABLES
    real :: vmax, vdotf, vdotfall, vdotv, fdotf
    real :: scale2
    real :: dtv, dtfm
    real, save :: alpha = ALPHA0
    real, save :: dtmax
    integer, save :: last_negative, nrestarts=0
    real, allocatable, save :: v(:,:)
    logical, save :: v_initialized = .false.
    real, save :: ddt = 0.5
    integer :: i,j
    external :: getEF

    if (.not.v_initialized)then
        allocate(v(3, natoms))
        last_negative = ntimestep
        v_initialized = .true.
        ddt = dt
        dtmax = ddt*TMAX
        DMAX = cg_drmax
        v(:, :) = 0
    endif

    if (restart)then
        last_negative = ntimestep
        nrestarts = nrestarts + 1
        alpha = ALPHA0
!may be the time step is wrong.
!        if(dt*(DT_SHRINK)**nrestarts > TMIN ) then
!            ddt=dt*(DT_SHRINK)**nrestarts
!        else
!            ddt = TMIN
!        endif
        do i=1,natoms
            do j=1,3
                v(j,i) = ftot(j,i)*ddt*mask(j,i)
            enddo
        enddo
        restart = .false.
        write(*, *) '++ ====================restart E no. ',nrestarts,' ===================='
        write(*,*)'++ ddt=',ddt,'alpha=',alpha
    else
        nrestarts = 0
    endif
    !// vdotfall = v dot f

    vdotf = sum(v * ftot)

    if (vdotf > 0.0) then
        vdotv = sum(v * ftot)
        fdotf = sum(ftot*ftot)
 
        if (fdotf == 0.0) then
            scale2 = 0.0
        else
            scale2 = alpha * sqrt(vdotv/fdotf)
        endif

        do i=1,natoms
            do j=1,3
                if(mask(j,i)==1)then
                    v(j,i) = (1.0 - alpha) * v(j,i) + scale2 * ftot(j,i)
                else
                    v(j,i)=0.0
                endif
            enddo
        enddo

        if (ntimestep - last_negative > DELAYSTEP) then
            ddt = MIN(ddt * DT_GROW, dtmax)
            alpha = alpha * ALPHA_SHRINK
        endif
    else
        last_negative = ntimestep
        if(ddt * DT_SHRINK >TMIN)then
            ddt = ddt * DT_SHRINK
        else
            ddt = TMIN
        endif
        alpha = ALPHA0
        v = 0
    endif

    !// no particle moves further than dmax
    dtv = ddt

    do i = 1, natoms
        vmax = MAX(abs(v(1, i)), abs(v(2, i)))
        vmax = MAX(vmax, abs(v(3, i)))
        if (dtv * vmax > dmax) dtv = dmax/vmax
    enddo
    !// Euler integration step

    do i = 1, natoms
        do j=1,3
            dtfm = dtv / xmass(imass(i))
            if(mask(j,i)==1)then
                ratom(:, i) = ratom(:, i) + dtv * v(:, i)
                v(:, i) = v(:, i) + dtfm * ftot(:, i)
            endif
        enddo
    enddo

end subroutine

subroutine getEF(e)
    use configuration, only: natoms
    use energy, only: etot
    use forces, only: ftot
    implicit none
    external :: getEnergyAndForces
    integer :: itmp
    real :: dtmp, e

    itmp = 0
    dtmp = 0
    call getEnergyAndForces(dtmp, itmp)
    e = etot
end subroutine

! Program Declaration
! ===========================================================================
subroutine FIRE_loop()
    call FIRE_loop_ex(0,0)
end subroutine

subroutine FIRE_loop_ex(nsteps,act_iter)
    use configuration, only: natoms, ratom, shifter, ishiftO
    use forces, only: ftot
    use energy, only: etot
    use optimization, only: cg_maxstep
    use charges, only:nzx
    use interactions, only: imass
    use options, only: iquench, iclassicMD
    use optimization, only: cg_maxstep, cg_minint !cgopt.optional input file
    implicit none
    integer :: nsteps,act_iter
    integer :: iatom
    integer :: istatus, nrestart
    real :: eclas, eMin, fMin(3,natoms), ratomMin(3,natoms), eOld, ecurrent, f(3, natoms)
    logical :: restart
    integer :: iter, ndownhillsteps
    !=====    EXTERNAL ROUTINES def in BFGS.f90:   ======
    !====================================================
    external :: writeOutRE, writeOutRFE, writeAnswerBas
    logical, external :: minimumReached
    real :: dist, fabs
    !====================================================
    write (*, *) '++   '
    write (*, *) '++  ======================================================= '
    write (*, *) '++              Begin FIRE minimization. '
    write (*, *) '++  ======================================================= '
    write (*, *) '++   '
    eOld = 100.0
    call getEF(ecurrent, f)
    eMin = etot
    ndownhillsteps = 0

    restart = .false.
    nrestart = 0
    do iter = 1 + act_iter, CG_MAXSTEP -1
        call iterate_FIRE(iter, restart)
        eOld = etot

        call getEF(ecurrent)
        if (etot > eOld)then
            restart = .true.
        else
            nrestart = 0
            call writeOutRFE(ratom, ftot, etot, 'coordinates and forces        ', iter)
            call writeAnswerBas(ratom, natoms)

            if(etot<eMin)then
                ndownhillsteps = ndownhillsteps +1
                eMin = etot
                fMin = ftot
                ratomMin = ratom
            endif
        endif
        !checking if the minimization could be stopped
        if (minimumReached(ratom, ftot, natoms, etot, eOld, iter))then
            write(*, *) '++ =============================='
            write(*, *) '++ ==     MINIMUM REACHED      =='
            write(*, *) '++ =============================='
            write(*, '(a,i5.1,a)') '++ ===== niter=',iter,' ======'
            exit
        endif
        if(nsteps>0 .and. nsteps <= ndownhillsteps)then
            write(*,*)'++ BACK TO BFGS MINIMIZATON AFTER ',nsteps,' DOWNHILL STEPS'
            exit
        endif
    enddo

    if (iter == CG_MAXSTEP)then
        write(*, *) '++ =============================='
        write(*, *) '++ ==   MINIMIZATION STOPPED   =='
        write(*, *) '++ =============================='
        write(*, '(a,i5.1,a,a)') '++ The minimization in the progress but you want to stop minimization after', CG_MAXSTEP, &
                    'iteration.', ' You can change this number in configuration file.'
    endif

    etot = eMin
    ftot = fMin
    ratom = ratomMin
    ! WRITE OUT OUTPUTS
    call writeOutRFE(ratom, ftot, etot, 'coordinates and forces        ', iter)
    call writeAnswerBas(ratom, natoms)

    if (iclassicMD == 2)then
        !debug
        eclas = etot
        call scf_loop(istatus)
        ! optionally perform post-processing (DOS etc.)
        call postscf()
        ! calculate the total energy
        call getenergy(istatus)
        write(*, *) '++ ==== ENERGY class ===='
        write(*, *) '++ etot=', eclas
    endif

    open (unit = 86, file = 'answer.bas', status = 'unknown')
    write (86, *) natoms
    write(*, *) '++ ==== ENERGY ===='
    write(*, *) '++ etot=', etot

    write(*, *) '++ ==== COORDINATES: ===='
    do iatom = 1, natoms
        write(*, 701) nzx(imass(iatom)), ratom(:, iatom)
        if (ishiftO == 1) then
            write (86, 700) nzx(imass(iatom)), ratom(:, iatom) - shifter(:)
        else
            write (86, 700) nzx(imass(iatom)), ratom(:, iatom)
        endif
    enddo
    close(unit = 86)

    700 format (2x, i2, 3(2x, f18.8))
    701 format ('++ ', 2x, i2, 3(2x, f18.8))
    return
end subroutine

