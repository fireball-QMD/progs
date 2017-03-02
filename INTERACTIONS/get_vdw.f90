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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

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
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine get_vdw ()
        use configuration
        use constants_fireball
        use forces
        use interactions
        use neighbor_map
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1, in2
        integer ineigh
        integer jatom
        integer matom
        integer mbeta

        real C6factor
        real distance
        real factor
        real dfactor
        real alpha
        real Rfactor
        real vdw_piece

        real, dimension (3) :: eta
        real, dimension (3) :: r1, r2

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Loop over the neighbors of each iatom.
        vdw = 0.0d0
        fvdw = 0.0d0
        do iatom = 1, natoms
         matom = neigh_vdw_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
         do ineigh = 1, neighn_vdw(iatom)       ! <==== loop over i's neighbors
          mbeta = neigh_b_vdw(ineigh,iatom)
          jatom = neigh_j_vdw(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

          if (ineigh .ne. matom) then
           distance = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2           &
     &                     + (r2(3) - r1(3))**2)
! JEL-UNITS
!           distance = distance/abohr
           if (distance .lt. 1.0d-03) then
            write (*,*) ' WARNING! The distance between atoms in the vdw '
            write (*,*) ' routine is nearly zero.  What is going on here? '
            write (*,*) ' iatom, jatom = ', iatom, jatom
            write (*,*) ' ineigh, matom = ', ineigh, matom
            write (*,*) ' Sorry, we must stop!  '
            stop
           end if

! Scale the van der Waals interactions according to this factor.
! Basically, closer to the nucleus this term is zero, but asymptotically
! allows the interactions to be 1/R^6.  See Elstner et al. J. Chem. Phys.
! v. 114 p. 5149 (2001)
           Rfactor = (R0(in1)**3 + R0(in2)**3)/(R0(in1)**2 + R0(in2)**2)
           C6factor = 2*C6(in1)*C6(in2)*p_alpha(in1)*p_alpha(in2)            &
     &               /(C6(in1)*p_alpha(in1)**2 + C6(in2)*p_alpha(in2)**2)
! dumping term
           alpha = -3.0d0*(distance/Rfactor)**7
           factor = (1 - exp(alpha))**4
           vdw_piece = - factor*C6factor/distance**6
! assemble the energy term
           vdw = vdw + vdw_piece

! assemble the forces
! the derivative of the dumping term
           dfactor = (-4.0d0*exp(alpha) + 12.0d0*(exp(alpha)**2)    &
     &      - 12.0d0*(exp(alpha)**3) + 4.0d0*(exp(alpha)**4))*7.0d0*alpha/distance
           eta(:) = (r2(:) - r1(:))/distance

           fvdw(:,iatom) = fvdw(:,iatom) + 6.0d0*vdw_piece*eta(:)/distance &
     &      	+ dfactor*C6factor/distance**6*eta(:)
           fvdw(:,jatom) = fvdw(:,jatom) - 6.0d0*vdw_piece*eta(:)/distance &
     &          - dfactor*C6factor/distance**6*eta(:)

          end if
         end do
        end do


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end
