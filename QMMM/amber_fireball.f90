subroutine amber_fireball

use qmmm_module

implicit none

integer :: k 

!crea un fichero fireball.in para que fireball pueda leer las opciones
!estas opciones tienen que ser variables definidas en qmmm_module 
!en el type qmmm_namelist y tendran la forma qmmm_nml%()

open (unit = 16, file = 'fireball.in', status = 'unknown')

write (16,600) "&option"

if (qmmm_nml%fb_ihorsfield .eq. 1) then
write (16,600) "imcweda = 0"
write (16,600) "ihorsfield = 1"
end if
write (16,601) "qstate = ", -1*qmmm_nml%qmcharge
write (16,600) "icluster = 1"
write (16,602) "iqout = ", qmmm_nml%fb_iqout
write (16,603) "max_scf_iterations = ", qmmm_nml%itrmax
write (16,602) "iensemble = ", qmmm_nml%fb_iensemble
write (16,602) "imdet = ", qmmm_nml%fb_imdet
write (16,605) "nddt = ", qmmm_nml%fb_nddt
write (16,602) "iqmmm = ", qmmm_nml%qmmm_int
write (16,602) "ifixcharge =",qmmm_nml%fb_ifixcharge
write (16,600) 'iquench = 4'
write (16,605) "tempfe =",qmmm_nml%fb_tempfe
write (16,602) "idftd3 =",qmmm_nml%fb_idftd3
write (16,602) "idipole =",qmmm_nml%fb_idipole
write (16,600) "&end"
write (16,600) "&output"
write (16,600) "iwrtcharges = 1"!, qmmm_nml%wrtcharges
write (16,602) "iwrtpop = ", qmmm_nml%fb_iwrtpop
write (16,602) "iwrtvel = ", qmmm_nml%fb_iwrtvel
write (16,602) "iwrteigen = ", qmmm_nml%fb_iwrteigen
write (16,602) "iwrtefermi = ", qmmm_nml%fb_iwrtefermi
write (16,602) "iwrtdos = ", qmmm_nml%fb_iwrtdos
write (16,602) "iwrtewf =",qmmm_nml%fb_iwrtewf
write (16,602) "iwrtatom =",qmmm_nml%fb_iwrtatom
write (16,600) "&end"
if (qmmm_nml%fb_iwrtewf .eq. 1) then
write (16,600) "&mesh"
write (16,603) "iewform = ", qmmm_nml%fb_iewform
write (16,603) "npbands = ", qmmm_nml%fb_npbands
write (16,608) "pbands = ", qmmm_nml%fb_pbands
write (16,600) "&end"
end if


write (*,*)  'fireball.in creado'

! crea el fichero input.bas y asi fireball ya sabe que tipos de atomos hay en el sitema
! y coge unicamente estos del Fdata


open (unit = 17, file = 'input.bas', status = 'unknown')

write (17,*) qmmm_struct%nquant_nlink
do k = 1, qmmm_struct%nquant_nlink
        write (17,700) qmmm_struct%iqm_atomic_numbers(k), qmmm_struct%qm_coords(:,k)
end do

write (*,*)  'input.bas creado'



! Format Statements
! ===========================================================================
600     format (a)!(a,f8.5)
601     format (a,i2)
602     format (a,i1)
603     format (a,i3)
604     format (a,i1)
605     format (a,i4)
608     format (a,199(i4,','),i3)
700     format (i2, 3(2x,f8.4))

end subroutine
