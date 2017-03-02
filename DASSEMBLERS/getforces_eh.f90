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

! forces_eh.f90
! Program Description
! ===========================================================================
!       This routine assemble forces for Extended Hubbard
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getforces_eh ()
 
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
          if (itheory_xc .eq. 1) then
           write (*,*) ' Dassemble on-site SNXC force contributions. '
           call Dassemble_snxc_on (nprocs, iordern)
           call Dassemble_snxc_2c (nprocs, iordern)
          end if
          if (itheory_xc .eq. 2) then
           write (*,*) ' Dassemble on-site OSLXC force contributions. '
           call Dassemble_olsxc_on (nprocs, iordern)
           call Dassemble_olsxc_2c (nprocs, iordern)
          end if
          
          write (*,*) ' Dassemble two-center extended-hubbard contributions. '
          call Dassemble_eh_2c (nprocs, my_proc, iordern)        

! ===========================================================================
!                               Dassemble_3c
! ===========================================================================
          write (*,*) '  '
          write (*,*) ' Dassemble three-center force contributions. '
          call Dassemble_3c (nprocs, iordern, igauss)

          write (*,*) ' Dassemble three-center PP force contributions. '
          call Dassemble_3c_PP (nprocs, iordern)


! Call the exchange-correlation interactions based on method chosen
          if (itheory_xc .eq. 1) then
           write (*,*) ' Dassemble off-site SN exchange-correlation forces. '
           call Dassemble_snxc_3c (nprocs, iordern, igauss)
          else if(itheory_xc .eq. 2) then
           write (*,*) ' Dassemble off-site OLS exchange-correlation forces. '
           call Dassemble_olsxc_3c (nprocs, iordern, igauss)
          end if
          
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
     &     iharmonic, iwrtfpieces)

! Reassign forces for tolerance testing. 
          ftotold = ftotnew
          ftotnew = ftot


! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine getforces_eh
 
