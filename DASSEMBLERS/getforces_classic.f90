! copyright info:
!
!                             @Copyright 2010
!                           Fireball Committee
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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! getforces_classic.f90
! Program Description
! ===========================================================================
! definiton of subroutine which call the right potential. This subroutine is
! called in ../DASSEMBLERS/getforces.f90 in case iclassicMD = 1 . The relevat
! empirical potentials going to be defined in files cMD_POTENTIAL_NAME.f90.
! In this file is defined Lennard-Jones potential.
! ===========================================================================
! Code written by:
! Zdenka Chromcova
! Institute of Physics of  the  AS CR,  v. v. i.
! Cukrovarnicka 10
! CZ-162 00 Praha 6
! Czech Republic
! email: chrom@fzu.cz
! webpage with description: http://nanosurf.fzu.cz/wiki/doku.php?id=classical_md
! ===========================================================================


! this prodcedure is called in module getforces if option iclassicMD = 1
subroutine getforces_classic()
    use forces, only: ftotnew, ftotold, ftot
    use energy, only: etot, etotper, etotnew, etotold
    use configuration, only: natoms, ratom
    use classicMD,only: distance2, Potential,getRGLforce,getLJforce
    implicit none

    interface
        ! debug - comapre numerical and analitical derivatives
        subroutine num_derivation(dx, ftot)
            use configuration, only: natoms
            real, intent(out) :: ftot(3, natoms)
            real, intent(in) :: dx
        end subroutine
    end interface

    real :: fnumder(3, natoms) !deleteme
    real :: e, distance(natoms, natoms) !, numftot(3, natoms)
    real :: dx

    fnumder = 0.0
    etotold = etotnew
    etot = 0.0

    distance(:,:) = -1

!	call num_derivation(1.0e-5,numftot)
    call getClassForce(ratom, ftot, e, distance)

!GLOBAL VARIABLES ARE SET HERE
!	ftot=numftot
    ftotold = ftotnew
    ftotnew = ftot
    etot = e
    etotper = e/natoms
    etotnew = etotper

!========= NUM DERIVATIVE =========
!    dx = 0.000001

!    call num_derivation(dx, ftot)
!    call getNumGrad(fqm, nparam, gradNum, dx)
!    call getNumGrad(fqm, nparam, gradNum2, -dx)
!    gradNum = (gradNum+gradNum2)/2

!    if (maxval(ftot - ftotnew)/maxval(ftotnew) > 0.05)then
!        write(*, *) 'numericka a analyticka derivace nesedi, ftot-fnum = '
!        write(*, *) ftot - ftotnew
!        write(*, *) '==============numerical======================='
!        write(*, *) ftot
!        write(*, *) '================analytical====================='
!        write(*, *) ftotnew
!        stop
!    endif
!    write(*,*)'==============OK==============='
!    etot=e
!    ftot = ftotnew
!========= NUM DERIVATIVE =========
end subroutine getforces_classic


!this should in individual file i.e. cMD_Lennard-Jones.dat
! ======= definition of local functions =======
subroutine getLJforce(ratom, f, e, distance)
    use classicMD, only: potential, distance2
    use neighbor_map, only: neigh_classic, neighn_classic, neigh_b_classic
    use configuration, only: xl, natoms, nspecies
	use interactions, only: imass
	
    implicit none
    !    type(PotentialParam),intent(in), dimension(nspecies,nspecies) :: potential
    real, intent(out) :: f(3, natoms)
    double precision, intent(inout) :: e
    real, intent(in) :: ratom(3, natoms)
    real, intent(in) :: distance(natoms, natoms)
    integer i, ineigh, k
    real :: eps, ro, cutoff, r1
    real, dimension(3) :: rneigh, ftmp

    e = 0.0
    do k = 1, natoms
        f(:, k) = (/0.0, 0.0, 0.0/)
        do i = 1, neighn_classic(k)
            ineigh = neigh_classic(i, k)
            rneigh = ratom(:, ineigh) + xl(:, neigh_b_classic(i, k))
            r1 = distance2(rneigh, ratom(:, k))
            cutoff = potential(imass(ineigh), imass(k)) % cutoff

            if (sum(xl(:, neigh_b_classic(i, k)) * xl(:, neigh_b_classic(i, k))) /= 0 &
                .and. r1 < 2 * potential(imass(ineigh), imass(k)) % cutoff)then
                eps = potential(imass(ineigh), imass(k)) % params(1)
                ro = potential(imass(ineigh), imass(k)) % params(2)
                !                eps =0.5*(potential(imass(ineigh), imass(ineigh)) % params(1)+potential(imass(k), imass(k)) % params(1))
                !                ro = 0.5*(potential(imass(ineigh), imass(ineigh)) % params(2)+potential(imass(k), imass(k)) % params(2))

                e = e + eps * ((ro * ro/r1)**6 - 2 * (ro * ro/r1)**3)
                ftmp = 12 * eps * ((ro**6)/(r1**4) - (ro**12)/(r1**7))*(ratom(:, k) - rneigh)
                f(:, k) = f(:, k) - 2 * ftmp * 1.602 !units: 1.602 - convert eV/nm to nN
            endif
        enddo
    enddo
end subroutine getLJforce

subroutine num_derivation(dx, ftotn)
    use configuration, only: natoms, ratom
    implicit none

    real, intent(out) :: ftotn(3, natoms)
    real, intent(in) :: dx
    real, dimension(3, natoms) :: rat2
    real :: distancet(natoms, natoms)
    double precision :: f(3, natoms), edx, e
    integer :: i, j

    rat2 = ratom

    do i = 1, natoms
        do j = 1, 3
            rat2(j, i) = ratom(j, i) + dx
            distancet(:,:) = -1.0
            call getClassForce(rat2, f, edx, distancet)
            rat2(j, i) = ratom(j, i) - dx
            distancet(:,:) = -1.0
            call getClassForce(rat2, f, e, distancet)
            ftotn(j, i) = -1.602 * (edx - e)/(2 * dx) !units eV/nm -> nN
            rat2(j, i) = ratom(j, i)
        enddo
    enddo
end subroutine

subroutine getClassForce(ratom, ftot, e, distance)
    use configuration, only: natoms
    use classicMD, only: Potential, getRGLforce, getLJforce, getforce_vdw
	implicit none
    real, dimension(3, natoms), intent(in) :: ratom
    real, dimension(3, natoms), intent(inout) :: ftot
    real, dimension(3, natoms) :: ftotn
    double precision, intent(inout) :: e
    real, dimension(natoms, natoms), intent(inout) :: distance

    e = 0.0
    ftot = 0.0
    select case (Potential(1, 1) % ptype)
    case (1)
        call getLJforce(ratom, ftot, e, distance)
    case (2)
        call getRGLforce(ratom, ftot, e, distance)
    case (3)
        call getTersoffforce(ratom, ftot, e, distance)
    case (4)
	call getforce_vdw(ratom, ftot, e)
!        call getSWforce(ratom, ftot, e, distance)
    case default
        write(*, *) 'ERROR - getforces_classic: I did not found the empirical potential ', Potential(1, 1) % ptype
        stop
    end select
end subroutine

function atom(i, ineigh, ratom, N)
    use neighbor_map, only: neigh_classic, neighn_classic, neigh_b_classic
    use configuration, only: xl
    real, dimension(3, N), intent(in) :: ratom
    real :: atom(3)
    integer, intent(in) :: i, ineigh, N
    if (neigh_classic(ineigh, i) == 0)then
        write(*, *) 'ERROR: atom(i,ineigh,ratom,N) - ratom(:,0), neigh_classic(ineigh, i) == 0!'
        if (neighn_classic(i) < ineigh)then
            write(*, *) 'neigh_classic(inegh=', ineigh, ',i=', i, ')=', neigh_classic(ineigh, i)
            write(*, *) 'neighn_classic=', neighn_classic(i), '<ineigh=', ineigh
        endif
        stop
    endif
    atom = ratom(:, neigh_classic(ineigh, i)) + xl(:, neigh_b_classic(ineigh, i))
end function atom

function distanceNeigh(i, ineigh, ratom, N)
    use classicMD, only:atom, distance2
    implicit none
    integer, intent(in) :: N, i, ineigh
    real, intent(in) :: ratom(3, N)
    real :: distanceNeigh, Rneigh(3)

    Rneigh = atom(i, ineigh, ratom, N)

    distanceNeigh = sqrt(distance2(ratom(:, i), Rneigh(:)))

end function distanceNeigh

function distance2(a, b)
    real :: distance2
    real, intent(in) :: a(3), b(3)
    distance2 = (a(1) - b(1))*(a(1) - b(1))+(a(2) - b(2))*(a(2) - b(2))+(a(3) - b(3))*(a(3) - b(3))
end function distance2

