! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold

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

 
! writeout_charges.f90
! Program Description
! ===========================================================================
!       Write out the charges for restart capabilities.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine writeout_charges (natoms, ifixcharge, iqout, iwrtcharges, iwrtdensity, basisfile, symbol,ab)
  use charges
  use outputs, only : iwrtdipole
  use density
  use interactions
  use neighbor_map
  use scf, only : scf_achieved, Kscf
  use MD, only : itime_step_g
  use configuration, only : ratom, xl
  implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
  integer, intent (in) :: natoms
  integer, intent (in) :: ifixcharge
  integer, intent (in) :: iqout
  integer, intent (in) :: iwrtcharges
  integer, intent (in) :: iwrtdensity
  integer, intent (in) :: ab !ab = 0: before mixing. ab=1: after mixing
  character (len=40) basisfile
  character (len=2), dimension (natoms) :: symbol

 
! Local Parameters and Data Declaration
! ===========================================================================
  real, parameter ::  Debye = 0.208194
! Local Variable Declaration and Description
! ===========================================================================
  integer iatom
  integer imu
  integer in1, in2
  integer ineigh
  integer inu
  integer issh
  integer jatom
  integer mbeta
  real    Qtot


 


! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! ****************************************************************************
!
! W R I T E    O U T    D E N S I T Y  
! ****************************************************************************
! Write out the density for the first atom.
  if (iwrtdensity .eq. 1) then
    do iatom = 1, natoms
      write (*,*) '  '
      write (*,*) ' Write out the density matrix for iatom = ', iatom
      write (*,*) ' Number of neighbors = ', neighn(iatom)
      in1 = imass(iatom)
      write (*,*) ' Number of orbitals on iatom, num_orb(in1) = ', num_orb(in1)
      do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass(jatom)
        write (*,*) '  '
        write (*,*) ' iatom, jatom = ', iatom, jatom
        write (*,*) ' Number of orbitals on jatom, num_orb(in2) = ', num_orb(in2)
        write (*,*) ' ----------------------------------------------'
        do imu = 1, num_orb(in1)
          write (*,400) (rho(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
        end do
      end do
    end do         
  end if ! iwrtdensity



  if (iwrtcharges .eq. 1) then 
    if (iqout .eq. 1 .or. iqout .eq. 3) then !Lowdin y NPA
      write (*,*) '  '
      write (*,*) '  '
      write (*,*) ' LOWDIN CHARGES (by shell): '
      write (*,500)
      write (*,501)
      write (*,500)
      do iatom = 1, natoms
        in1 = imass(iatom)
        write (*,502) iatom, symbol(iatom), nssh(in1),(Qin(issh,iatom), issh = 1, nssh(in1))
      end do
      write (*,500)
      write (*,*) '  '
      write (*,*) '  '
      write (*,*) ' Total Lowdin charges for each atom: '
      write (*,500)
      do iatom = 1, natoms
        write (*,503) iatom, symbol(iatom), QLowdin_TOT(iatom)
      end do
      write (*,500)
      write (*,*) '  '
      write (*,*) '  '
    end if ! Low NAP

    if ((iqout .eq. 2) .or. (iqout .eq. 4) .or. (iqout .eq. 7)) then !Mullk y Mulliken-dipole
      write (*,*) '  '
      write (*,*) '  '
      write (*,*) ' MULLIKEN CHARGES (by shell): '
      write (*,500)
      write (*,501)
      write (*,500)
      do iatom = 1, natoms
       in1 = imass(iatom)
       write (*,502) iatom, symbol(iatom), nssh(in1), (Qin(issh,iatom), issh = 1, nssh(in1))
      end do

      write (*,500)
      write (*,*) '  '
      write (*,*) '  '
      write (*,*) ' Total Mulliken charges for each atom: '
      write (*,500)
      do iatom = 1, natoms
       write (*,503) iatom, symbol(iatom), QMulliken_TOT(iatom)
      end do
      write (*,500)
      write (*,*) '  '
      write (*,*) '  '
    end if !end if (iqout .eq. 2 y 4)

    write (*,500)
    write (*,*) '  '
  end if  !end if (iwrtcharges .eq. 1)

! ****************************************************************************
! W R I T E    O U T    T E M P O R A L    S E R I E S   O F    C H A R G E S
! ****************************************************************************

  if (iwrtcharges .eq. 2) then
    if ( scf_achieved ) then
      open(unit = 333, file = 'CHARGES.xyz', status = 'unknown', position = 'append')
      open(unit = 334, file = 'CHARGES_series', status = 'unknown', position = 'append')
      write(333,*) natoms
      write(333,*) 'Time step = ', itime_step_g
      do iatom = 1, natoms
        Qtot=0
        in1 = imass(iatom)
        do issh = 1,nssh(in1)
          Qtot = Qtot+Qin(issh,iatom)
        end do
        write(333,601) (Qin(issh,iatom), issh = 1, nssh(in1)), Qtot
        write(334,400,advance="no") Qtot
      end do !end do iatom = 1,natoms
      close(334)
      close(333)
    end if !end if ( scf_achieved )
  end if !end if (iwrtcharges .eq. 2)


  if (iwrtcharges .eq. 3) then !!! cambiar el formato a n C xyz charge
    if (ab .eq. 1) then
      open(unit = 333, file = 'CHARGES.dat', status = 'unknown', position = 'append')
    else ! ab .eq. 0
       open(unit = 333, file = 'CHARGES_no_mix.dat', status = 'unknown', position = 'append')
    end if ! if ab .eq. 1
    write(333,*) natoms
    write(333,*) 'Time step = ', itime_step_g, 'Kscf = ', Kscf
    do iatom = 1, natoms
      Qtot=0 
      in1 = imass(iatom)
      do issh = 1,nssh(in1)
        Qtot = Qtot+Qout(issh,iatom)
      end do
      !write(333,601) (Qout(issh,iatom), issh = 1, nssh(in1)), Qtot, Qtot-Q0_TOT(iatom)
    end do !end do iatom = 1,natoms
    close(333)     
  end if !end if (iwrtcharges .eq. 3)
  
    
! ****************************************************************************
! W R I T E    O U T    C H A R G E S    F I L E
! ****************************************************************************
! Open the file CHARGES which contain the Lowdin charges for restart purposes.
  if (ifixcharge .eq. 0 .and. wrtout) then
    open (unit = 21, file = 'CHARGES', status = 'unknown')
    write (21,600) natoms, basisfile, iqout
    do iatom = 1, natoms
      in1 = imass(iatom)
      write (21,601) (Qin(issh,iatom), issh = 1, nssh(in1))
    end do
    close (unit = 21)
  end if

! ****************************************************************************
! W R I T E   T H E   D I P O L E   (For testing purposes)
! ****************************************************************************
  if (iwrtdipole .gt. 0) then
  
    if (iqout .ne. 7) then
      call Dipole_proyection()
    end if
 
     open( unit = 667, file = 'Charges_and_Dipoles', status = 'unknown')
     write(667,*)   '+++++++++++++++++++ N E W   S T E P ++++++++++++++++++' 
     write(667,444) 'dip_TOT',dipTot_x/Debye, dipTot_y/Debye, dipTot_z/Debye, dipTot_tot/Debye
     write(667,444) 'dip_Qin',dipQin_x/Debye, dipQin_y/Debye, dipQin_z/Debye, dipQin_tot/Debye
     write(667,444) 'dip_Qout',dipQout_x/Debye, dipQout_y/Debye, dipQout_z/Debye, dipQout_tot/Debye
     write(667,444) 'dip_Int', dipIntra_x/Debye, dipIntra_y/Debye, dipIntra_z/Debye, dipIntra_tot/Debye
     write(667,444) 'dip_res',dipRes_x/Debye, dipRes_y/Debye,dipRes_z/Debye, dipRes_tot/Debye
     close(667)
     
  end if !end if (iwrtdipole .gt. 0)

! END OF WRITING THE DIPOLE
!******************************************************************************
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
400     format (9f9.4)
500     format (70('='))
501     format (2x, ' Atom # ', 2x, ' Type ', 2x, ' Shells ', 1x,' Charges ')
502     format (3x, i5, 7x, a2, 5x, i2, 4x, 8(1x, f5.2))
503     format (3x, i5, 7x, a2, 4x, f10.4)
600     format (2x, i5, 2x, a40, 2x, i2)
601     format (2x, 10f14.8)
444     format (a7,4f10.4)
445     format (a2,6f10.4)

        return
        end

subroutine load_partial_charges()
  use dimensions
  use interactions
  use charges
  use configuration
  integer issh
  integer in1
  integer iatom
  do iatom = 1, natoms
    in1 = imass(iatom)
    Q_partial(iatom) = Q0_TOT(iatom)
    do  issh = 1, nssh(in1)
      Q_partial(iatom) = Q_partial(iatom) - Qout(issh,iatom)
    end do
  ! write(*,*)'Q_partial',iatom,Q_partial(iatom)
  end do
end subroutine load_partial_charges



