subroutine coords_to_fireball

use configuration
use interactions
use qmmm_module

!vatom = (ratom-qmmm_struct%qm_coords)/0.0005
ratom = qmmm_struct%qm_coords
vatom = qmmm_struct%qm_vel



!        write(*,*) 'estamos en coords_to_fireball'
!       
!        write(*,*) 'vatom=', vatom
!        write(*,*) 'ratom=', ratom
!
!
!        write (*,*) '  '
!        write (*,*) '  '
!        write (*,*) ' Atom Coordinates from ratom File: '
!        write (*,200)
!        write (*,201)
!        write (*,200)
!        do iatom = 1, qmmm_struct%nquant_nlink
!         in1 = imass(iatom)
!         symbol(iatom) = symbolA(in1)
!         write (*,202) iatom, symbol(iatom), ratom(:,iatom), imass(iatom)
!        end do
!
!
!        write (*,*) '  '
!        write (*,*) '  '
!        write (*,*) ' Atom Coordinates from qmmm_struct%qm_coords File: '
!        write (*,200)
!        write (*,201)
!        write (*,200)
!        do iatom = 1, qmmm_struct%nquant_nlink
!         in1 = imass(iatom)
!         symbol(iatom) = symbolA(in1)
!         write (*,202) iatom, symbol(iatom), qmmm_struct%qm_coords(1,iatom), qmmm_struct%qm_coords(2,iatom), qmmm_struct%qm_coords(3,iatom), imass(iatom)
!        end do
!
!
!        write (*,*) ' Atom Velocities from vatom File: '
!        write (*,200)
!        write (*,201)
!        write (*,200)
!        do iatom = 1, qmmm_struct%nquant_nlink
!         in1 = imass(iatom)
!         symbol(iatom) = symbolA(in1)
!         write (*,202) iatom, symbol(iatom), vatom(:,iatom), imass(iatom)
!        end do
!
!
!        write (*,*) ' Atom Velocities from qmmm_struct%qm_vel File: '
!        write (*,200)
!        write (*,201)
!        write (*,200)
!        do iatom = 1, qmmm_struct%nquant_nlink
!         in1 = imass(iatom)
!         symbol(iatom) = symbolA(in1)
!         write (*,202) iatom, symbol(iatom), qmmm_struct%qm_vel(:,iatom), imass(iatom)
!        end do
!
!
!        write (*,200)
!        write (*,*) '  '
!


! Format Statements
! ===========================================================================
200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 5x,   &
     &              ' x ', 8x, ' y ', 8x, ' z ', 6x, ' Species # ')
202     format (3x, i5, 2x,f9.3)
203     format (i4)


end subroutine coords_to_fireball
