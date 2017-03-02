! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

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

! tdbc.f90
! Program Description
! ===========================================================================
!       This routine sets the boundary conditions on the edge, where ions move
! the new wf-coefficients are evaluated to avoid instability
!    C(t+dt) = S(t+dt)^(-1/2)*S(t)^(1/2)*C(t)
!
! ===========================================================================
! Code rewritten by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine set_tdbc ()

        use tdse
        use dimensions
        use interactions
        use kpoints

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Local Parameters and Data Declaration
! ===========================================================================

       integer  ielec
       integer  imu
       integer  inu
       integer  ikpoint

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================
! NOTE: s12psi has been calculated in the subroutine tddenmat for
! former ionic configuration

!        do ikpoint = 1, nkpoints
!         do ielec = 1, nelec
!          do imu = 1, norbitals
!           do inu = 1, norbitals
!            psi(imu,ielec,ikpoint) = psi(imu,ielec,ikpoint) +              &
!     &         sm12(imu,inu,ikpoint)*s12psi(inu,ielec,ikpoint)
!           end do ! inu
!          end do ! imu
!         end do ! ielec
!       end do ! ikpoint


! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine set_tdbc

