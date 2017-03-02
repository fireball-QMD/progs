subroutine coords_forces_charges_to_amber(escf)

use interactions
use configuration
use forces
use qmmm_module
use charges
use energy
use options

implicit none

integer :: iatom, issh, in1
real , intent(out) :: escf


qmmm_struct%dxyzqm = -ftot*23.061d0
escf = etot*23.061d0
qmmm_struct%qm_vel = vatom

        do iatom = 1, natoms
         qmmm_struct%Qneutral_TOT(iatom) = 0
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
          qmmm_struct%Qneutral_TOT(iatom) = qmmm_struct%Qneutral_TOT(iatom) + Qneutral(issh,in1)
         end do
        end do
	
	if (qmmm_nml%fb_iqout .eq. 1) then
 	  qm2_struct%scf_mchg = -(QLowdin_TOT - qmmm_struct%Qneutral_TOT) ! qm2_struct%scf_mchg
	else if (qmmm_nml%fb_iqout .eq. 2) then
	  qm2_struct%scf_mchg = -(QMulliken_TOT - qmmm_struct%Qneutral_TOT) ! qm2_struct%scf_mchg
        else if (qmmm_nml%fb_iqout .eq. 3) then
	  qm2_struct%scf_mchg = -(QLowdin_TOT - qmmm_struct%Qneutral_TOT) ! qm2_struct%scf_mchg
        end if

!        if (iqmmm .eq. 0) then
!           qmmm_struct%dxyzcl = 0.0d0
!        endif

!        write (*,*) ' Atom Velocities from ratom File: '
!        write (*,200)
!        write (*,201)
!        write (*,200)
!        do iatom = 1, qmmm_struct%nquant_nlink
!         in1 = imass(iatom)
!         symbol(iatom) = symbolA(in1)
!         write (*,202) iatom, symbol(iatom), qmmm_struct%qm_vel(:,iatom), imass(iatom)
!        end do


        write (*,200)
        write (*,*) '  '



! Format Statements
! ===========================================================================
200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 5x,   &
     &              ' x ', 8x, ' y ', 8x, ' z ', 6x, ' Species # ')
202     format (3x, i5, 7x, a2, 3(2x,f9.3), 7x, i2)
203     format (2x, ' Atom # ', 2x, ' Type ', 5x,   &
     &              'fx ', 8x, 'fy ', 8x, 'fz ', 6x, ' Species # ')

204     format (3x, i5, 7x, a2, 2x, f9.3, 7x, i2)


end subroutine coords_forces_charges_to_amber
