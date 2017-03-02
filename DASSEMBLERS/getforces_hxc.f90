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

! forces_hxc.f90
! Program Description
! ===========================================================================
!       This routine assemble forces for Horsfield
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getforces_hxc ()

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

! Call the exchange-correlation interactions based on method chosen
! (i.e. itheory_xc).
          write (*,*) ' Dassemble Horsfield exchange-correlation forces. '
          call Dassemble_hxc_2c (nprocs, iordern, itheory, igauss)


          if (itheory .eq. 1) then
           write (*,*) ' Dassemble two-center DOGS force contributions. '
           call Dassemble_ca_2c (nprocs, iordern)
          endif

! ===========================================================================
!                               Dassemble_3c
! ===========================================================================
          write (*,*) '  '
          write (*,*) ' Dassemble three-center force contributions. '
          call Dassemble_3c (nprocs, iordern, igauss)

          write (*,*) ' Dassemble three-center PP force contributions. '
          call Dassemble_3c_PP (nprocs, iordern)

          if (itheory .eq. 1) then
           write (*,*) ' Dassemble three-center DOGS force contributions. '
           call Dassemble_ca_3c (nprocs, iordern, igauss)
           write (*,*) ' Dassemble three-center long-range contributions. '
           call Dassemble_lr (nprocs, iordern)
          end if

! Call the exchange-correlation interactions based on method chosen
! (i.e. itheory_xc).

          write (*,*)'Dassemble Horsfield exchange-correlation forces.'
          call Dassemble_hxc_3c (nprocs, iordern, itheory, igauss)
!

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
        end subroutine getforces_hxc

