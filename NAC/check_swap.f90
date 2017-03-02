! copyright info:
!
!                             @Copyright 2009
!                FAST (Fireball Atomic Simulation Techniques)
! West Virginia University - James P. Lewis, Chair
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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! overlap_sign.f90
! Program Description
! ===========================================================================
! Calculate non-adiabatic coupling(d_{jk}  contribution V.d_{jk}
! using Kohn-Sham states at different
! time steps, and compares with the equivalent contribution obtained
! directly using the non-adiabtic couplings calculated in
! nacouplings.f90
!
! ===========================================================================
! Code written by Enrique Abad Gonzalez
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine check_swap (itime_step,Kscf)


        use configuration
        use nonadiabatic
        use density
        use kpoints
        use interactions
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer itime_step
        integer Kscf


! Local Parameters and Data Declaration
! ===========================================================================
       real, parameter :: hbar = 0.65822d0


! Local Variable Declaration and Description
! ===========================================================================
        integer it
        integer imu, inu, jmu,jnu, iorbital
        integer in1, in2, in3
        integer iatom, jatom
        integer isorp, interaction
        integer ix
        integer ia
        integer ik, ij
        integer ikpoint
!        integer max_inu, iband, jband, ioccupy_oldimu, ioccupy_oldinu
!        real foccupy_oldimu, foccupy_oldinu, max_sumb, check_otherstates
        real max_sumb, check_otherstates
        real y, rcutoff_i, rcutoff_j, range
        real diff
        real delta
        integer, dimension (nele) :: banda_swap, max_inu
        integer, dimension (norbitals) :: ioccupy_swap
        real, dimension (norbitals) :: foccupy_swap
        real, dimension (3) :: r1, r2, r21, sighat
        real, dimension (3,3,3) :: deps
        real, dimension (3,3) :: eps
        real, dimension (numorb_max, numorb_max) :: sx
        real, dimension (3,numorb_max, numorb_max) :: spx
        real, dimension (norbitals, norbitals) :: s
        real, dimension (nele, nele) :: suma
!        real, dimension (nele, nele) :: sumb
        complex aim
        complex a0
        complex a1



! ===========================================================================
        aim = cmplx(0.0d0, 1.0d0)
        a0 = cmplx(0.0d0, 0.0d0)
        a1 = cmplx(1.0d0, 0.0d0)

! ===========================================================================
! Calculate non-adiabatic couplings using Kohn-Sham states at different
! time steps
! (1) Define overlap matrix between orbitals mu at time "t" and orbitals
! nu at time "t+dt"
!print *, 'empezamos'
        if ( (itime_step .eq. 1) .and. (Kscf .eq. 1) ) then
         allocate (sumb(nele,nele))
         allocate (ratom_old(3,natoms))
         allocate (bbnkre_old(norbitals,norbitals,nkpoints))
         allocate (blowre_old(norbitals,norbitals,nkpoints))
         ratom_old = ratom
         bbnkre_old = bbnkre
         blowre_old = blowre
        end if
        call overlap_numeric()

    do ikpoint = 1, nkpoints
      do imu = 1, nele
        max_inu(imu) = 0
        max_sumb = 0.0d0
        check_otherstates = 0.0d0
        do inu = 1, nele
          check_otherstates = check_otherstates + (sumb(imu,inu)*sumb(imu,inu))
          if ( abs(sumb(imu,inu)) .gt. max_sumb ) then
            max_sumb = abs(sumb(imu,inu))
            max_inu(imu) = inu
          end if
	end do
        if ( abs(check_otherstates - 1.0 ) .gt. 0.3 ) then
          write (*,*) 'Looks like there is a swap between one state and another one not considered in mdet.input'
          write (*,*) 'We are going to stop, and encourage you to consider more states in mdet.input'
          STOP
        end if
        if ( max_inu(imu) .ne. imu ) then
          write (*,*) 'The most similar state to',imu,'is',max_inu(imu),'not itself.'
          write (*,*) 'Probably state swapping has ocurred; we will swap them in map_ks'
          write (*,*) 'and change their occupancy accordingly'
        end if
      end do

          ! swapping... Better a separate subroutine?
      do imu = 1, nele
         banda_swap(imu) = map_ks(imu)
         ioccupy_swap(imu) = ioccupy_na(map_ks(imu),ikpoint)
         foccupy_swap(imu) = foccupy_na(map_ks(imu),ikpoint)
         print *, 'imu,map_ks,iocc,focc',imu,map_ks(imu),ioccupy_na(map_ks(imu),ikpoint),foccupy_na(map_ks(imu),ikpoint)
      end do
      do imu = 1, nele
         map_ks(imu) = banda_swap(max_inu(imu))
         ioccupy_na(map_ks(imu),ikpoint) = ioccupy_swap(imu)
         foccupy_na(map_ks(imu),ikpoint) = foccupy_swap(imu)
         print *, 'imu,max_inu,map_ks,iocc,foccNEW',imu,max_inu(imu),map_ks(imu),ioccupy_na(map_ks(imu),ikpoint),foccupy_na(map_ks(imu),ikpoint)
      end do
     end do ! end do ikpoints




! Check if this overlap matrix is not unit matrix and change in consequently
! NOW DONE IN OVERLAP_SIGN!!
!      do imu = 1, nele
!            if ( sumb(imu,max_inu(imu)) .lt. -0.1 ) then
!              write (*,*) 'The arbitrary sign of wavefunctions of this and previous time step are different'
!              write (*,*) 'We will change the sign in order to have the same sign in all time steps!!'
!              do iorbital = 1, norbitals
!                bbnkre(iorbital,map_ks(imu),ikpoint) = - bbnkre(iorbital,map_ks(imu),ikpoint)
!              end do
!            end if
!      end do
!     end do ! end do ikpoints



! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format ('c_na',2i4,f8.4,2(f7.3))
300     format ('NAC-SUMS',2i4,2f8.4)
301     format ('S(t,tprime)',2i4,1f8.4)
400     format ('S',4f7.3)


        return
        end subroutine check_swap

