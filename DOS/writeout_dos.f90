! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

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
! C.Gonzalez Pascual, UAM, Spain

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


! writeout_dos.f90
! Program Description
! ===========================================================================
! This routine writes the calculates and writes the DOS of the atoms we want
! (from natom_beg to natom_end). If iwrttip writes the DOS (in a special format
! in the tip_e_str.inp we need in the STM calculation
! ===========================================================================
! Code written by:
! Cesar Gonzalez Pascual
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madird SPAIN
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
! modified by P. Jelinek
! ==========================================================================

        subroutine writeout_dos( )

        use dimensions
        use interactions
        use neighbor_map
        use module_dos
        use configuration
        use charges

        implicit none

! Argument Declaration and Description
! Input

! Output

! Local Parameters and Data Declaration
! ===========================================================================
       real, parameter :: pi=3.141592653589793230d0

! Local Variable Declaration and Description
! ==========================================================================
        character(80)   :: titf1
        character(80)   :: titf2
        character(3)    :: tit1
        integer         :: iatom
        integer         :: jatom
        integer         :: iorb
        integer         :: jorb
        integer         :: norbi
        integer         :: norbj
        integer         :: ie
        integer         :: f1
        integer         :: f2
        integer         :: in1
        integer         :: in2
        integer         :: fall
        integer         :: nsteps_tip
        real            :: energy
        real            :: ener_tip_bottom
        real            :: dos_tot
        real            :: dos_min
        real            :: dos_atoms
        real            :: charge_atoms
        real, dimension(numorb_max)  :: dos_orb
        real, dimension(natom_beg:natom_end,numorb_max)  :: dos_orb_aux
        real, dimension(natom_beg:natom_end)  :: charge_atom
        real, dimension(natom_beg:natom_end,numorb_max)  :: charge_orb

! Procedure
! ===========================================================================
! Initialize some things
! Open the dos files 
      do iatom = natom_beg,natom_end         !CGP 17-V-2002
        write (tit1,'(i3.3)') iatom          !CGP 17-V-2002
        titf1 = 'dens_'//tit1//'.dat'        !CGP 17-V-2002
        titf2 = 'char_'//tit1//'.dat'        !CGP 17-V-2002
        f1 = 10 + iatom                      !CGP 17-V-2002
        open (unit=f1, file=titf1)           !CGP 17-V-2002
        f2 = 10 + iatom + natom_end          !CGP 17-V-2002
        open (unit=f2, file=titf2)           !CGP 17-V-2002
      end do 
      if (iwrttip .ge. 1) open(9,file='tip_e_str.inp',status='unknown')
      
        fall = f2+1
        open (unit=fall, file='dens_TOT.dat')

      dos_orb = 0.0d0
      dos_orb_aux = 0.0d0
      nsteps_tip = 0
      charge_orb = 0.0d0
      norbi = 0
      energy = ener_beg

! Loop over energy
      do ie = 1, nener
! increase energy step counter for tip
        if(iwrttip .ge. 1) then
          if (energy .gt. ener_min .and. energy .lt. ener_max) nsteps_tip = nsteps_tip + 1
          if (nsteps_tip .eq. 1) ener_tip_bottom = energy - efermi
        end if ! end if (iwrttip. ge. 1)
        norbi = 0
        dos_atoms = 0.0d0
        charge_atoms = 0.0d0
! Loop over atoms
        do iatom = natom_beg, natom_end
          in1 = imass(iatom)
          f1 = 10 + iatom
          f2 = 10 + iatom + natom_end
          dos_orb = 0.0d0
          do iorb = 1, num_orb(in1)
            norbi = norbi + 1
            dos_orb(iorb) = -1.0d0/pi * imag(green(norbi,norbi,ie))
            dos_min = min(dos_orb(iorb), dos_orb_aux(iatom,iorb))
            charge_orb(iatom,iorb) = charge_orb(iatom,iorb) +          &
              ener_step*(dos_min + abs(dos_orb(iorb)-dos_orb_aux(iatom,iorb))/2.0d0)
            if(iwrttip .ge. 1) then
              if(energy .gt. ener_min .and. energy .lt. ener_max) then
                norbj= 0
                do jatom = natom_beg,natom_end
                  in2 = imass(jatom)
                  do jorb = 1,num_orb(in2)
                    norbj = norbj + 1
                    write(9,*) imag(green(norbi,norbj,ie)),             &
     &                          real(green(norbi,norbj,ie))
                    write(9,*) -imag(green(norbi,norbj,ie)),            &
     &                          real(green(norbi,norbj,ie))
                    write(9,*) -1.0d0/pi*imag(green(norbi,norbj,ie))
                  end do
                end do
              end if ! end if (energy gt.ener_min) 
            end if ! end if (iwrttip. ge. 1)
            dos_orb_aux(iatom,iorb) = dos_orb(iorb)
          end do             
          dos_tot = sum(dos_orb(:))
          charge_atom(iatom) = sum(charge_orb(iatom,:))
          dos_atoms = dos_atoms + dos_tot
          charge_atoms = charge_atoms + charge_atom(iatom)
          write (f1,'(20F10.4)') energy,dos_orb(1:num_orb(in1)),          &
     &                            dos_tot,charge_atom(iatom)        
          write (f2,'(20F10.4)') energy,charge_orb(iatom,1:num_orb(in1)), &
     &                            charge_atom(iatom)    
        end do
        write (fall,'(20F10.4)') energy,dos_atoms,charge_atoms
        energy = energy + ener_step
      end do

! close dos files
      close (unit=fall)
      do iatom = natom_beg,natom_end                   ! CGP 17-V-2002
         close (unit=10+iatom)           ! CGP 17-V-2002
         close (unit=10+iatom+natom_end) ! CGP 17-V-2002
      end do     

! write down tip_g_str.inp file which contains structural information about tip
      if(iwrttip .ge. 1) then 

	open(10,file='tip_g_str.inp',status='unknown')
! write range of atoms
	write (10,*) natom_end-natom_beg+1, natom_end-natom_beg+1
! atomic coordinates
        do iatom = natom_beg, natom_end
          in1 = imass(iatom)
          if (ishiftO .eq. 1) then
	   write (10,100) ratom(:,iatom) - shifter(:), in1, num_orb(in1)
          else
	   write (10,100) ratom(:,iatom), in1, num_orb(in1)
          endif
        end do
! number of shells
        write (10,*) nssh(1:nspecies)
! l for each shell
        do in1 = 1, nspecies
         write (10,*) lssh (1:nssh(in1),in1)
        end do
! energy range
        write (10,200) ener_tip_bottom, ener_step*nsteps_tip, nsteps_tip         
! open files
        close (9)
        close (10)
        stop
       endif
! Format Statements
! ===========================================================================
100     format (3f16.9, 2x, 2i2) 
200     format (2f16.6, i5) 

      return
  
      end subroutine writeout_dos
