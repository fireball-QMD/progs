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


! formc_compact.f90
! Program Description
! ===========================================================================
!       This subroutine builds the Localized Wave Functions in a compact
! format, centered on ATOMS (searching for atoms within a cutoff radius 
! rcutoff), and assigns a RANDOM initial guess to start the CG minimization.
!
! Criterion to build LWF's: 
! 1) Method of Kim et al: use more localized orbitals than occupied orbitals.
!    We assign the minimum number of orbitals so that there is place for more 
!    electrons than those in the system; for instance:
!      H:        1 LWF
!      C,Si:     3 LWF's
!      N:        3 LWF's
!      O:        4 LWF's
!      ...
! 2) Method of Ordejon et al.: number of localized orbitals equal to number 
!    of occupied orbitals. Only available when each atom has an EVEN number of
!    electrons.

! ===========================================================================
! This code was originally provided by Pablo Ordejon.
! Code rewritten by:
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
        subroutine formc_compact (natoms, nspecies, ioptionlwf, ratom,       & 
     &                            icrows, isendstart)
        use charges
        use constants_fireball
        use interactions
        use neighbor_map
        use ordern
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icrows
        integer, intent (in) :: isendstart
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies

        integer, intent (in) :: ioptionlwf ! : Build LWF's according to:
                                           ! 0 = Read blindly from disk
                                           ! 1 = Functional of Kim et al.
                                           ! 2 = Functional of Ordejon-Mauri

        real, intent (in), dimension (3, natoms) :: ratom

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer icount
        integer imu
        integer in1, in2
        integer index
        integer indexa, indexb
        integer indexi
        integer ineigh
        integer inu
        integer jatom
        integer mmu
        integer numloc

        integer, dimension (:), allocatable :: indexloc

        real cvalue
        real distance
        real factor
        real xnorm

        logical skip_it
 
! Allocate Arrays
! ===========================================================================
        allocate (indexloc(norbitals))
 
! Procedure
! ===========================================================================
! Initialize some arrays to zero.
        index = 0
        numc_local = 0.0d0
        listc_local = 0.0d0
        c_compact_local = 0.0d0

! Loop over all the atoms in the central cell.
        do iatom = 1, natoms
         in1 = imass(iatom)

! Determine how many LWF's depending on the atomic species (or the number of
! electrons).
         if (ioptionlwf .eq. 1) then
          indexi = nelectron(iatom)/2
          if (nzx(in1) .eq. 1) indexi = 1
          if (nzx(in1) .eq. 7) indexi = indexi + 1           ! For HMX
          if (nzx(in1) .eq. 8) then
           icount = 0
           do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh,iatom)
            distance = sqrt((ratom(1,iatom) - ratom(1,jatom))**2             &
     &                    + (ratom(2,iatom) - ratom(2,jatom))**2             &
     &                    + (ratom(3,iatom) - ratom(3,jatom))**2)
            if (distance .lt. 1.5d0) icount = icount + 1
           end do
           if (icount .le. 2) indexi = indexi - 1
          end if
         else if (ioptionlwf .eq. 2) then
          if ((nelectron(iatom)/2)*2 .ne. nelectron(iatom)) then
           write (*,*) ' Wrong Order-N functional option in formc_compact. '
           write (*,*) ' You can only use the functional of Ordejon-Mauri '
           write (*,*) ' for atoms with an even number of electrons. '
           stop
          end if
          indexi = nelectron(iatom)/2
         else
          write (*,*) ' Wrong functional option in formc_compact '
          stop
         end if

! Loop over LWF's centered on iatom.
	 do indexb = 1, indexi
	  index = index + 1

! Clear list of atoms considered within localization range.
          indexloc = 0
          numloc = 0

! Initialize stuff
          xnorm = 0.0d0
 
! Loop over the neighbors of iatom within rcutoff.
          do ineigh = 1, neighn(iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)

! Limit LWF's to those within established cutoff. 
           distance = sqrt((ratom(1,iatom) - ratom(1,jatom))**2              &
     &                   + (ratom(2,iatom) - ratom(2,jatom))**2              &
     &                   + (ratom(3,iatom) - ratom(3,jatom))**2)       
           if (distance .lt. rcutoff_lwf*abohr) then 

! Check if jatom has already been included in current lwf.
            skip_it = .false.
            do indexa = 1, numloc
             if (jatom .eq. indexloc(indexa)) skip_it = .true. 
            end do
            if (.not. skip_it) then
             numloc = numloc + 1
             indexloc(numloc) = jatom

! Loop over orbitals of jatom
             do imu = 1, num_orb(in2)
	      mmu = imu + degelec(jatom)
              if (mmu .ge. isendstart .and. mmu .le. icrows) then
               mmu = mmu - isendstart + 1
	       numc_local(mmu) = numc_local(mmu) + 1
	       listc_local(numc_local(mmu),mmu) = index

! Assign random guess for orbitals in iatom.
               if (jatom .eq. iatom) then
                call initguess (natoms, nspecies, iatom, imu, in1, cvalue)
                c_compact_local(numc_local(mmu),mmu) = cvalue
                xnorm = xnorm + cvalue**2
               end if
              end if
             end do
            end if
           end if
          end do
          if (xnorm .lt. 1.0d-4) xnorm = 1.0d0

! Normalize LWF's. Normalize to one if functions are expected to be occupied, 
! 0.5 if half occupied and 0.1 if empty)
          factor = 1.0d0
          if (ioptionlwf .eq. 1) then
           if (2*(nelectron(iatom)/2) .eq. nelectron(iatom)                  &
     &         .and. indexb .eq. indexi) factor = sqrt(0.1d0)
           if (2*(nelectron(iatom)/2) .ne. nelectron(iatom)                  &
     &         .and. indexb .eq. indexi) factor = sqrt(0.5d0)
          end if
          if (iatom .lt. natoms) then 
           do imu = 1 + degelec(iatom), degelec(iatom + 1)
            if (imu .ge. isendstart .and. imu .le. icrows) then
             mmu = imu - isendstart + 1
             do inu = 1, numc_local(mmu)
              if (listc_local(inu,mmu) .eq. index) then
               c_compact_local(inu,mmu) =                                    &
     &          c_compact_local(inu,mmu)*factor/sqrt(xnorm)
              end if
             end do
            end if
           end do
          else
           do imu = 1 + degelec(iatom), norbitals
            if (imu .ge. isendstart .and. imu .le. icrows) then
             mmu = imu - isendstart + 1
             do inu = 1, numc_local(mmu)
              if (listc_local(inu,mmu) .eq. index) then
               c_compact_local(inu,mmu) =                                    &
     &          c_compact_local(inu,mmu)*factor/sqrt(xnorm)
              end if
             end do
            end if
           end do
          end if
         end do
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (indexloc)
 
! Format Statements
! ===========================================================================
 
        return
        end
