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


! forces_mcweda.f90
! Program Description
! ===========================================================================
!       This routine assemble forces for McWeda
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getforces_mcweda ()

        use options
        use outputs
        use mpi_main
        use configuration
        use forces
        !use qmmm_module, only : qmmm_struct

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
          !write (*,*) '  '
          !write (*,*) '  '
          !write (*,*) ' ***************************************************** '
          !write (*,*) ' Dassemble two-center force contributions. '
          call Dassemble_2c (nprocs, iordern, igauss)

          !write (*,*) ' Dassemble two-center PP force contributions. '
          call Dassemble_2c_PP (nprocs, iordern)
!bias
           if (ibias .eq. 1) then
             call Dassemble_bias (nprocs, iordern)
           endif



! Call the exchange-correlation interactions based on method chosen
! (i.e. itheory_xc).
          if (itheory_xc .eq. 1) then
           !write (*,*) ' Dassemble on-site SNXC force contributions. '
           if (itheory .eq. 1) then
             call Dassemble_ca_snxc_on (nprocs, iordern)
             call Dassemble_ca_snxc_2c (nprocs, iordern)
           else
             call Dassemble_snxc_on (nprocs, iordern)
             call Dassemble_snxc_2c (nprocs, iordern)
           endif
          end if
          if (itheory_xc .eq. 2) then
           !write (*,*) ' Dassemble on-site OSLXC force contributions. '
           if (itheory .eq. 1) then
             call Dassemble_ca_olsxc_on (nprocs, iordern)
             call Dassemble_ca_olsxc_2c (nprocs, iordern)
           else
             call Dassemble_olsxc_on (nprocs, iordern)
             call Dassemble_olsxc_2c (nprocs, iordern)
           endif
          end if

          if (itheory .eq. 1) then
           !write (*,*) ' Dassemble two-center DOGS force contributions. '
           if (idipole .eq. 0) call Dassemble_ca_2c (nprocs, iordern)
           if (idipole .eq. 1) call Dassemble_ca_2c_dip (nprocs, iordern)
          endif

! ===========================================================================
!                               Dassemble_3c
! ===========================================================================
          !write (*,*) '  '
          !write (*,*) ' Dassemble three-center force contributions. '
          call Dassemble_3c (nprocs, iordern, igauss)

          !write (*,*) ' Dassemble three-center PP force contributions. '
          call Dassemble_3c_PP (nprocs, iordern)

          if (itheory .eq. 1) then
           !write (*,*) ' Dassemble three-center DOGS force contributions. '
           if (idipole .eq. 0) call Dassemble_ca_3c (nprocs, iordern, igauss)
           if (idipole .eq. 1) call Dassemble_ca_3c_dip (nprocs, iordern, igauss)
           !write (*,*) ' Dassemble three-center long-range contributions. '
           if (idipole .eq. 0) call Dassemble_lr (nprocs, iordern)
           if (idipole .eq. 1) call Dassemble_lr_dip (nprocs, iordern)
           if (iqmmm .eq. 1) then
             !write (*,*) ' Dassemble three-center qm/mm contributions. '
             if (idipole .eq. 0) call Dassemble_qmmm (nprocs, iordern)
             if (idipole .eq. 1) call Dassemble_qmmm_dip (nprocs, iordern)
           else
             flrew_qmmm = 0.0d0
!             qmmm_struct%dxyzcl = 0.0d0
           end if
          end if

! Call the exchange-correlation interactions based on method chosen
          if (itheory_xc .eq. 1) then
           !write (*,*) ' Dassemble off-site SN exchange-correlation forces. '
           if (itheory .eq. 1) then
! Include gaussians
            call Dassemble_ca_snxc_3c (nprocs, iordern, igauss)
           else
            call Dassemble_snxc_3c (nprocs, iordern, igauss)
           endif
          else if(itheory_xc .eq. 2) then
           !write (*,*) ' Dassemble off-site OLS exchange-correlation forces. '
           if (itheory .eq. 1) then
            call Dassemble_ca_olsxc_3c (nprocs, iordern, igauss)
           else
            call Dassemble_olsxc_3c (nprocs, iordern, igauss)
           endif
          end if

          !write (*,*) ' ***************************************************** '

! ============================================================================
!                                assemble_F
! ============================================================================
! Call assemble_F: This program assembles all the forces we have calculated
! in Dassemble_2c and Dassemble_3c, and assemble_usr.
          !write (*,*) '  '
          !write (*,*) '  '
          !write (*,*) ' ***************************************************** '
          !write (*,*) ' Assemble all force contributions. '
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
        end subroutine getforces_mcweda

