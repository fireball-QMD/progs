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


! set_dimensions.f90
! Program Description
! ===========================================================================
!
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
        subroutine set_dimensions (natoms, ioptionlwf, ratom, ztot)
        use charges
        use constants_fireball
        use interactions
        use neighbor_map
        use ordern
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
!Input
        integer, intent (in) :: ioptionlwf
        integer, intent (in) :: natoms
 
        real, intent (in) :: ztot

        real, intent (in), dimension (3, natoms) :: ratom
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer icount
        integer imu
        integer in1, in2
        integer index1
        integer index2
        integer indexa
        integer indexb
        integer indexct
        integer indexFt
        integer indexi
        integer ineigh
        integer inu
        integer iplace
        integer issh
        integer jatom
        integer jneigh
        integer katom
        integer mmu
        integer nnu
        integer numloc

        integer, dimension (:), allocatable :: indexloc
        integer, dimension (:), allocatable :: indexlocF
        integer, dimension (:), allocatable :: numF
        integer, dimension (:), allocatable :: numT1_local
        integer, dimension (:), allocatable :: numT2_local
        integer, dimension (:, :), allocatable :: listT2_local
 
        real distance 
        real rhmax 
         
        logical skip_it

! Allocate Arrays
! ===========================================================================
        allocate (indexloc (norbitals))
        allocate (indexlocF (norbitals))
        allocate (numT1_local (norbitals))
        allocate (numT2_local (norbitals))
        allocate (listT2_local (norbitals, norbitals))
        allocate (numF (norbitals))
 
! Procedure
! ===========================================================================
! ***************************************************************************
!                                  nhmax
! ***************************************************************************
! Find the maximum number of non-zero elements for the sparse Hamiltonian
! and overlap matrices.
! First initialize the index counter to zero.
        indexloc = 0
        nhmax = 0
        numT2_local = 0
        listT2_local = 0

! Loop over all the atoms in the central cell.
        do iatom = 1, natoms
         in1 = imass(iatom)

! Loop over the neighbors of iatom
         index1 = 0
         index2 = 0
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)

          iplace = 1 + degelec(jatom)
          if (indexloc(iplace) .eq. 0) then
           indexloc(iplace) = index1 + 1

! Loop over orbitals of jatom
           do inu = 1, num_orb(in2) 
            index1 = index1 + 1

! Loop over orbitals of iatom
            do imu = 1, num_orb(in1)
             mmu = imu + degelec(iatom)
             listT2_local(index1,mmu) = inu + degelec(jatom)
            end do
           end do 
          else 
           index2 = indexloc(iplace)

! Loop over orbitals of jatom and iatom
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             index2 = index2 + 1
            end do
           end do
          end if
         end do
         if (index1 .gt. nhmax) nhmax = index1
         if (index2 .gt. nhmax) nhmax = index2

! Set counter of number of non-zero elements
         do imu = 1, num_orb(in1)
          mmu = imu + degelec(iatom)
          numT2_local(mmu) = index1
         end do

! Reinitialize control vector - index
         do imu = 1, index1
          iplace = 1 + degelec(iatom)
          inu = listT2_local(imu,iplace)
          indexloc(inu) = 0
         end do
        end do

 
! ***************************************************************************
! Find ncmax, nctmax, nFmax
! ***************************************************************************
! Find the maximum number of non-zero elements for the wavefunction 
! coefficients.
        index1 = 0
        indexloc = 0
        numF = 0
        numT1_local = 0
        ncmax = 0
        nctmax = 0
        nFmax = 0
        nFtmax = 0

! Loop over all the atoms in the central cell.
        do iatom = 1, natoms
         in1 = imass(iatom)

! Determine how many LWF's depending on the atomic species (or the number of 
! electrons).
         if (ioptionlwf .eq. 1) then
          indexi = nelectron(iatom)/2
          if (nzx(in1) .eq. 1) indexi = 1
          if (nzx(in1) .eq. 7) indexi = indexi + 1       ! For HMX
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
           write (*,*) ' Wrong Order-N functional option in set_dimensions. '
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
!        write (*,*) ' iatom, nzx(in1), indexi = ', iatom, nzx(in1), indexi
         do indexb = 1, indexi
          index1 = index1 + 1

! Clear list of atoms considered within localization range.
          indexloc = 0
          indexlocF = 0
          numloc = 0
          indexFt = 0 
          indexct = 0 

! Loop over the neighbors of iatom within rcutoff.
          do ineigh = 1, neighn(iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)

! Limit LWF's to those within established cutoff.  
           distance = sqrt((ratom(1,iatom) - ratom(1,jatom))**2             &
     &                   + (ratom(2,iatom) - ratom(2,jatom))**2             &
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
              numT1_local(mmu) = numT1_local(mmu) + 1
              if (numT1_local(mmu) .gt. ncmax) ncmax = numT1_local(mmu)
              indexct = indexct + 1

! Find out structure of F and Ft matrices.
              do inu = 1, numT2_local(mmu)
               nnu = listT2_local(inu,mmu)
               if (indexlocF(nnu) .eq. 0) then 
                indexlocF(nnu) = 1
                numF(nnu) = numF(nnu) + 1
                if (numF(nnu) .gt. nFmax) nFmax = numF(nnu)
                indexFt = indexFt + 1
               end if
              end do
             end do
            end if
           end if
          end do
          if (indexFt .gt. nFtmax) nFtmax = indexFt 
          if (indexct .gt. nctmax) nctmax = indexct
         end do
        end do


! ***************************************************************************
!                                  nbands
! ***************************************************************************
! Check that there is an even number of electrons in the system.  This is 
! required if the linear-scaling option is to be used.
        if (2*(int(ztot + 1.0d-3)/2) .ne. int(ztot + 1.0d-3)) then
         write (*,*) ' In formc_compact: Wrong total charge; odd charge: ', ztot
         write (*,*) ' Charge must be EVEN to use Order-N option! '
         stop
        end if

        nbands = index1
!       if (nbands .gt. sum(nelectron)/2) then
!        write (*,*) ' Number of LWFs (nbands) larger than half number of '
!        write (*,*) ' electrons. Something is wrong here. '
!        write (*,*) ' nbands = ', nbands
!        write (*,*) ' sum(nelectron)/2 = ', sum(nelectron)/2
!        stop
!       end if


! Deallocate Arrays
! ===========================================================================
        deallocate (indexloc)
        deallocate (indexlocF)
        deallocate (numT1_local)
        deallocate (numT2_local, listT2_local)
        deallocate (numF)
 
! Format Statements
! ===========================================================================
 
        return
        end
