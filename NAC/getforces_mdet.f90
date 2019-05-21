! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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


! forces_mdet.f90
! Program Description
! ===========================================================================
!       This routine assemble forces for MDET/McWeda
!
! JOM : adapted from getforces_mcweda to also calculate 
! the gradient of the Hamiltonian
! G < mu | H | nu > and overlap contributions for the calculation
! of the nonadiabatic couplings
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getforces_mdet ()

        use options
        use outputs
        use mpi_main
        use configuration
        use forces
        use nonadiabatic

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
          !write (*,*) ' In getforces_mdet '
          if (itheory_xc .ne. 2) then
          !write (*,*) ' must stop, only ready for itheory_xc = 2 '
          stop
          end if
!         if (itheory .ne. 0) then
!         !write (*,*) ' must stop, only ready for itheory = 0 (Harris) '
!         stop
!         end if
          if (icluster .ne. 1) then
          !write (*,*) ' must stop, only ready for icluster = 1 '
          stop
          end if
! JOM-add allocate nonadiabatic variables
!         !write(*,*) ' JOM : later allocate+reallocate if doing MD'
          call allocate_nac(natoms)
! JOM-warning : later allocate+reallocate if doing MD (changing number
! of neighbors)
!         call reallocate_nac(natoms)
! JOM-end allocate nonadiabatic variables
!--------------------------------------------------------------------
          !write (*,*) '  '
          !write (*,*) 'assemble gover = < Grad mu | nu > contributions.'
          call assemble_G_S (nprocs, iordern)
!--------------------------------------------------------------------
          !write (*,*) ' ********************************************** '
          !write (*,*) ' Dassemble two-center force contributions. '
          call Dassemble_2c_mdet (nprocs, iordern, igauss)

          if (itheory .eq. 1) then
           !write (*,*) ' Dassemble two-center DOGS force contributions.'
           if (idipole .eq. 1) then
            call Dassemble_ca_2c_mdet_dip (nprocs, iordern)
           else
            call Dassemble_ca_2c_mdet (nprocs, iordern)
           end if !end if idipole .eq. 1
          endif
!--------------------------------------------------------------------

          !write (*,*) ' Dassemble two-center PP force contributions. '
          call Dassemble_2c_PP_mdet (nprocs, iordern)
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
             call Dassemble_ca_olsxc_on_mdet (nprocs, iordern)
             call Dassemble_ca_olsxc_2c_mdet (nprocs, iordern)
           else
             call Dassemble_olsxc_on_mdet (nprocs, iordern)
             call Dassemble_olsxc_2c_mdet (nprocs, iordern)
           endif
          end if

! ===========================================================================
!                               Dassemble_3c
! ===========================================================================
          !write (*,*) '  '
          !write (*,*) ' Dassemble three-center force contributions. '
          call Dassemble_3c_mdet (nprocs, iordern, igauss)

          !write (*,*) ' Dassemble three-center PP force contributions. '
          call Dassemble_3c_PP_mdet (nprocs, iordern)

          if (itheory .eq. 1) then
           !write (*,*) ' Dassemble three-center DOGS force. '
           if (idipole .eq. 1) then
            call Dassemble_ca_3c_mdet_dip (nprocs, iordern,igauss)
            call Dassemble_lr_mdet_dip (nprocs, iordern)
           else !else idipole .eq. 1
            call Dassemble_ca_3c_mdet (nprocs, iordern, igauss)
            call Dassemble_lr_mdet (nprocs, iordern)
           !write (*,*) ' Dassemble three-center long-range '
           end if !end if idipole .eq. 1 
           if (iqmmm .eq. 1) then
             !write (*,*) ' Dassemble three-center qm/mm contributions. '
           if (idipole .eq. 1) then
                   !Dassemble_qmmm_mdet_dip
             call Dassemble_qmmm_mdet_dip (nprocs, iordern)
           else ! else idipole .eq. 1
             call Dassemble_qmmm_mdet (nprocs, iordern)
           end if ! end if idipole .eq. 1
           else !else if iqmmm .eq. 1
             flrew_qmmm = 0.0d0
           end if ! end if iqmmm .eq. 1
          end if !itheory

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
            call Dassemble_ca_olsxc_3c_mdet (nprocs, iordern, igauss)
           else
            call Dassemble_olsxc_3c_mdet (nprocs, iordern, igauss)
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
          !write (*,*) ' ********************************************* '
          !write (*,*) ' Assemble all force contributions. '
          call assemble_F (natoms, itheory, itheory_xc, igauss, ivdw,       &
     &     iharmonic, ibias, iwrtfpieces)

! Reassign forces for tolerance testing.
          ftotold = ftotnew
          ftotnew = ftot
!
!JOM-add : NAC!
! ============================================================================
!                                NAC
! ============================================================================
! Call nacouplings.f90 : Calculates the non-adiabatic couplings between 
! Kohn-Sham states
          call nacouplings()
! Instead of calling nacouplings, now we are going to calculate them using
! the numerical derivative
!         call delta_t_ks()
          call deallocate_nac(natoms)
!JOM-end

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine getforces_mdet

