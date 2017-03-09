! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
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


! readpressure.f90
! Program Description
! ===========================================================================
!       This routine reads in the information from the pressure.optional
! file. This information is needed for doing constant pressure runs.
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
        subroutine readpressure (iPRESSURE, Pext, Wcell, icell_dynamics,     &
     &                           iquenchs, iquenchl, betaq, alphaq,          &
     &                           volume_file, lkinetic_file, skinetic_file,  &
     &                           latvec_file)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: iPRESSURE
 
! Output
        integer, intent (out) :: icell_dynamics
        integer, intent (out) :: iquenchl
        integer, intent (out) :: iquenchs
 
        real, intent (out) :: alphaq
        real, intent (out) :: betaq
        real, intent (out) :: Pext
        real, intent (out) :: Wcell
 
        character (len=30), intent (out) :: volume_file
        character (len=30), intent (out) :: lkinetic_file
        character (len=30), intent (out) :: skinetic_file
        character (len=30), intent (out) :: latvec_file
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
        if (iPRESSURE .eq. 1) then
         write (*,*) '  '
         write (*,100)
         write (*,*) ' You are doing constant pressure. '
         write (*,*) ' Read information from pressure.optional file. '
 
         open (unit = 79, file = 'pressure.optional', status = 'old')
 
         write (*,*) '  '
         write (*,*) ' Insert pressure in GPa: '
         read (79,*) Pext
         write (*,*) ' Pext = ', Pext, ' (GPa) '
         Pext = Pext/160.0d0
         write (*,*) ' Pext = ', Pext, ' (eV/A**3) '
         write (*,*) '  '
         write (*,*) ' Insert the mass of the cell (e.g. 1000). '
         read (79,*) Wcell
         write (*,*) ' The mass of the cell is = ', Wcell
 
! Three types of dynamics.
! icell_dynamics =         1               2              3
!                Parrinello-Rahman  Projected force   Wentzcovitch
 
! Parrinello Rahman dynamics is described in PRL 45, 1196 (1980).
! Projected cell dynamics is where we compute de/dh and use only that part
! of it projected along the directions of h1, h2, and h3. This way we maintain
! the box shape (actually its angles), but the size and relative length change.
! Thus a box can go to tetragonal or orthorhombic, but not hexagonal.
! Wentzcovitch dynamics means Eq. 8 of Phys. Rev. B 44, 2358 (1991).
! It is similar to PR dynamics, but you multiply a the final part of
! the "force" by (SAB(daggger)*SAB)^-1. This uses a reference SAB.
         write (*,*) '  '
         write (*,*) ' The variable icell_dynamics tells us what kind of cell '
         write (*,*) ' dynamics we want. '
         write (*,*) ' icell_dynamics = 1 Parrinello Rahman dynamics '
         write (*,*) ' icell_dynamics = 2 Projected force dynamics '
         write (*,*) ' icell_dynamics = 3 Wentzcovitch dynamics '
         write (*,*) ' Default is 3 '
         write (*,*) ' Insert icell_dynamics: '
         read (79,*) icell_dynamics
         write(*,*) ' icell_dynamics = ', icell_dynamics
         if (icell_dynamics .lt. 1 .or. icell_dynamics .gt. 3) then
           write (*,*) ' bad icell_dynamics, check pressure.optional '
           stop
         end if
 
! You may want to quench your atoms and not to quench your lattice:
         write (*,*) '  '
         write (*,*) ' Insert scaled velocity iquenchs (yes/no -> 1/0):'
         read (79,*) iquenchs
         write (*,*) ' iquenchs = ',iquenchs
         write (*,*) '  '
         write (*,*) ' Insert the lattice iquenchl (yes/no -> 1/0): '
         read (79,*) iquenchl
         write (*,*) ' iquenchl = ',iquenchl
 
! When the scaled temperature reachres a maximum
! sdot(1,ix,iy) = betaq*sdot(1,ix,i)
         write (*,*) '  '
         write (*,*) ' When the scaled temperature reachres a maximum, '
         write (*,*) ' then sdot(1,ix,iatom) = betaq*sdot(1,ix,iatom). '
         write (*,*) ' Insert betaq (0.6 is a reasonable choice): '
         read (79,*) betaq
         write (*,*) ' betaq = ', betaq
 
! To quench the box in the PRESSURE run we use a slightly different technique.
! When the lattice temperature reaches a maximum
! hdot(1,ix,iy) = alphaq*hdot(1,ix,iy)
         write (*,*) '  '
         write (*,*) ' When the lattice temperature reachres a maximum'
         write (*,*) ' hdot(1,ix,iy) = alphaq*hdot(1,ix,iy)'
         write (*,*) ' Insert alphaq (0.4 is a reasonable choice): '
         read (79,*) alphaq
         write (*,*) ' alphaq = ', alphaq
 
! Get the file names.
         write (*,*) '  '
         write (*,*) ' Insert filenames: '
         write (*,*) ' (volume_file, skinetic_file, '
         write (*,*) '  lkinetic_file, latvec_file) '
 
         read (79,200) volume_file
         write (*,*) ' volume_file = ', volume_file
         read (79,200) skinetic_file
         write (*,*) ' skinetic_file = ', skinetic_file
         read (79,200) lkinetic_file
         write (*,*) ' lkinetic_file = ', lkinetic_file
         read (79,200) latvec_file
         write (*,*) ' latvec_file = ', latvec_file
         write (*,100)
 
         close (unit = 79)
        end if
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (a30)
 
        return
        end
