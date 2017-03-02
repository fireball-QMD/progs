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

! postete.f90
! Program Description
! ===========================================================================
!      This subroutine assemble density matrix, energy & forces after
! time-dependent evoution of wavefunction with fixed ionic positions
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
        subroutine postete (itime_step)
        use energy
        use tdse


        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
       integer, intent (in) :: itime_step

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

! assemble the density matrix & band structure energy
        call tddenmat ()

! calculate energy
        write (*,*) ' Assemble total energy'
        call getenergy (itime_step)

! calc forces
        call getforces ()

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine postete

