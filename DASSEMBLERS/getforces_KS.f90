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

! forces_KS.f90
! Program Description
! ===========================================================================
!       This routine assemble forces for Kohn-Sham
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getforces_KS ()

        use options
        use outputs
        use mpi_main
        use configuration
        use forces

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================

! ============================================================================
!                              Dassemble_2c
! ============================================================================
          write (*,*) '  '
          write (*,*) '  '
          write (*,*) ' ***************************************************** '
          write (*,*) ' Dassemble two-center force contributions. '
          call Dassemble_2c (nprocs, iordern, igauss)

          write (*,*) ' Dassemble two-center PP force contributions. '
          call Dassemble_2c_PP (nprocs, iordern)

! ===========================================================================
!                               Dassemble_3c
! ===========================================================================
          write (*,*) '  '
          write (*,*) ' Dassemble three-center force contributions. '
          call Dassemble_3c (nprocs, iordern, igauss)

          write (*,*) ' Dassemble three-center PP force contributions. '
          call Dassemble_3c_PP (nprocs, iordern)

          write (*,*) ' ***************************************************** '

! ============================================================================
!                                assemble_F
! ============================================================================
! Call assemble_F: This program assembles all the forces we have calculated
! in Dassemble_2c and Dassemble_3c, and assemble_usr.
          write (*,*) '  '
          write (*,*) '  '
          write (*,*) ' ***************************************************** '
          write (*,*) ' Assemble all force contributions. '
          call assemble_F (natoms, itheory, itheory_xc, igauss, ivdw,       &
     &     iharmonic, ibias, iwrtfpieces)

! Reassign forces for tolerance testing.
          ftotold = ftotnew
          ftotnew = ftot


! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine getforces_KS

