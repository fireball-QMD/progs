! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
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


! get_vdw.f90
! Program Description
! ===========================================================================
! 	This routine calculates the vdW energy and forces on the atoms in the system.
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
!
! modified by P. Jelinek (June 2009)
! modified by Z. Chromcova (April 2013)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine getforce_vdw(ratom, ftot, etot)
    use configuration, only: natoms, nspecies, xl
!    use constants_fireball
!    use forces
!    use interactions
!    use neighbor_map
    use configuration, only: natoms
    use interactions, only: imass
    use classicMD, only: PotentialParam, potential
    use neighbor_map, only: neigh_classic, neighn_classic, neigh_b_classic
    implicit none

    ! Argument Declaration and Description
    ! ===========================================================================
    ! Input

    ! Output

    ! Local Parameters and Data Declaration
    ! ===========================================================================

    ! Local Variable Declaration and Description
    ! ===========================================================================
    real, intent(inout) :: etot, ftot(3,natoms)
    real, intent(in) :: ratom(3,natoms)
    integer :: iatom, in1, in2
    integer :: ineigh, jatom, matom, mbeta

    real :: C6factor, distance, factor, dfactor, alpha, Rfactor, vdw_piece

    real, dimension (3) :: eta
    real, dimension (3) :: r1, r2

    ! Allocate Arrays
    ! ===========================================================================

    ! Procedure
    ! ===========================================================================
    ! Loop over the neighbors of each iatom.
    etot = 0.0d0
    ftot = 0.0d0
    do iatom = 1, natoms
        matom = iatom !neigh_vdw_self(iatom)
        r1(:) = ratom(:, iatom)
        in1 = imass(iatom)
        do ineigh = 1, neighn_classic(iatom) ! <==== loop over i's neighbors
            mbeta = neigh_b_classic(ineigh, iatom)
            jatom = neigh_classic(ineigh, iatom)
            r2(:) = ratom(:, jatom) + xl(:, mbeta) !ratom(:, neigh_classic(ineigh, i)) + xl(:, neigh_b_classic(ineigh, i))
            in2 = imass(jatom)

            if ((iatom .ne. jatom) .or. (sum(xl(:, mbeta)*xl(:, mbeta))/=0)) then
                distance = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2 +(r2(3) - r1(3))**2)
                ! JEL-UNITS
                !           distance = distance/abohr
                if (distance .lt. 1.0d-03) then
                    write (*, *) ' WARNING! The distance between atoms in the vdw '
                    write (*, *) ' routine is nearly zero.  What is going on here? '
                    write (*, *) ' iatom, jatom = ', iatom, jatom
                    write (*, *) ' ineigh, matom = ', ineigh, matom
		    write (*, *) ' xl(:,mbeta) = ', xl(:, mbeta)
                    write (*, *) ' Sorry, we must stop!  '
                    stop
                end if

                ! Scale the van der Waals interactions according to this factor.
                ! Basically, closer to the nucleus this term is zero, but asymptotically
                ! allows the interactions to be 1/R^6.  See Elstner et al. J. Chem. Phys.
                ! v. 114 p. 5149 (2001)
!                Rfactor = (R0(in1)**3 + R0(in2)**3)/(R0(in1)**2 + R0(in2)**2)
                 Rfactor = ((potential(in1,in1)%params(3))**3 + (potential(in2,in2)%params(3))**3)/&
                            ((potential(in1,in1)%params(3))**2 + (potential(in2,in2)%params(3))**2)
                
!                C6factor = 2 * C6(in1) * C6(in2) * p_alpha(in1) * p_alpha(in2) &
                 C6factor = 2 * potential(in1,in1)%params(1) * potential(in2,in2)%params(1) *&
                            potential(in1,in1)%params(2) * potential(in2,in2)%params(2) &
!                           & /(C6(in1) * p_alpha(in1)**2 + C6(in2) * p_alpha(in2)**2)
                            & /(potential(in1,in1)%params(1) * (potential(in1,in1)%params(2))**2 + &
                            potential(in2,in2)%params(1) * (potential(in2,in2)%params(2))**2)
                ! dumping term
                alpha = -3.0d0 * (distance/Rfactor)**7
                factor = (1 - exp(alpha))**4
                vdw_piece = -factor * C6factor/distance**6
                ! assemble the energy term
                etot = etot + vdw_piece

                ! assemble the forces
                ! the derivative of the dumping term
                dfactor = (-4.0d0 * exp(alpha) + 12.0d0 * (exp(alpha)**2) &
                & -12.0d0 * (exp(alpha)**3) + 4.0d0 * (exp(alpha)**4)) * 7.0d0 * alpha/distance
                eta(:) = (r2(:) - r1(:))/distance

                ftot(:, iatom) = ftot(:, iatom) - (6.0d0 * vdw_piece * eta(:)/distance &
                & +dfactor * C6factor/distance**6 * eta(:))

                ftot(:, jatom) = ftot(:, jatom) + (6.0d0 * vdw_piece * eta(:)/distance &
                & +dfactor * C6factor/distance**6 * eta(:))
            end if
        end do
    end do
    ! Deallocate Arrays
    ! ===========================================================================

    ! Format Statements
    ! ===========================================================================

    return
end

