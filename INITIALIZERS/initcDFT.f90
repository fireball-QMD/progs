! Copyright info:
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


! initcDFT.f90
! Program Description
! ===========================================================================
! ===========================================================================
! Code written by:
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine initcDFT ()

   use scf
   use interactions
   use charges 
   use MD
   use options
   use kpoints

   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Output

! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
      real diffq

! Format Statements
      write (*,*) '  ---- Initialize cDFT -----'
      open (unit = 33, file = 'cDFT.optional', status = 'old')
      read (33,*) id_hole
      read (33,*) id_elec
      read (33,*) occup_elec
      close (unit = 33)
      occup_hole = 1.0d0 - occup_elec
      write (*,200) id_hole, occup_hole
      write (*,300) id_elec, occup_elec

! allocate arrays
      allocate (wf_hole(norbitals, nkpoints))
      allocate (wf_elec(norbitals, nkpoints))

      cDFT_active = .false.
      
! only 2 MD steps (ground & excited state) & no forces
!      nstepf = 3
!      iforce = 0

! ===========================================================================
200     format (2x, ' Hole @ state : ',i8,' with occupancy ',f8.4)
300     format (2x, ' Electron @ state : ',i8,' with occupancy ',f8.4)

   return
 end subroutine initcDFT

