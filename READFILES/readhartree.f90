! copyright info:
!
!                             @Copyright 2009
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
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


! readhartree.f90
! Program Description
! ===========================================================================
! This routine reads the hartree.input file in order to calculate
! the Hartree-Fock correction to the molecular gap
!

! ===========================================================================
! Code written by:
! Cesar Gonzalez Pascual & Enrique Abad Gonzalez
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madrid SPAIN
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
!===========================================================================
!
! Program Declaration
!===========================================================================
        subroutine readhartree (nspecies,natoms)

        use dimensions
        use interactions
        use neighbor_map
        use hartree_fock

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

       integer, intent(in) :: nspecies
       integer, intent(in) :: natoms
! Output

! Local Variable Declaration and Description
! ==========================================================================
       integer iatom
       integer in1,in2
       integer iorb,jorb
       integer imixoption
       integer icont
       real    Utemp, Jtemp

! Procedure
! ===========================================================================
! Initialize some things

! Opening the file
       open (unit = 121, file = 'hartree-fock.input')

! Read dos.input
       read (121,*) natomshf
! Read the number of atoms (the first ant the final atoms)
       read(121,*) natomhf_beg,natomhf_end
! Read if we want to calculate a media of the U's and J's of each atom (in a future of each orbital)
       read(121,*) imixoption, betha
       allocate(Uisigma(numorb_max, numorb_max, nspecies))
       allocate(Jijsigma(numorb_max, nspecies, numorb_max, nspecies))
       allocate(Jialpha(numorb_max, natomhf_beg:natomhf_end))
       allocate(nij(numorb_max, natoms, numorb_max, natoms))
!print *, 'leemos las Us'
       do in1 = 1, nspecies
        do iorb = 1, num_orb(in1)
!print *, 'leyendo... iorb,num_orb(in1),in1',iorb,num_orb(in1),in1
          read(121,*) Uisigma(iorb,1:num_orb(in1),in1)
!print *,    Uisigma(iorb,1:num_orb(in1),in1)
        end do
       end do
!print *, 'leemos las Js'
       do in1 = 1, nspecies
         do in2 = in1, nspecies
           do iorb = 1, num_orb(in1)
             read(121,*) Jijsigma(iorb,in1, 1:num_orb(in2),in2)
!print *, 'leyendo... iorb, in1, 1:num_orb(in2), in2', iorb, in1, num_orb(in2), in2
!print *, Jijsigma(iorb,in1, 1:num_orb(in2),in2)
             Jijsigma(1:num_orb(in2),in2, iorb,in1)= Jijsigma(iorb,in1, 1:num_orb(in2),in2)
           end do
         end do
       end do

       if ( imixoption .eq. 1 ) then
       do in1 = 1, nspecies
         Utemp = 0.0
	 icont = 0
           do iorb = 1, num_orb(in1)
	     do jorb = 1, num_orb(in1)
               Utemp = Utemp +  Uisigma(iorb,jorb,in1)
	       icont = icont + 1
             end do
	   end do
	   Utemp = Utemp/(real(icont))
!	   print *, 'U,especie',Utemp,in1
           do iorb = 1, num_orb(in1)
	     do jorb = 1, num_orb(in1)
               Uisigma(iorb,jorb,in1) = Utemp
             end do
	   end do
	end do

       do in1 = 1, nspecies
         do in2 = 1, nspecies
	 Jtemp = 0.0
	 icont = 0
           do iorb = 1, num_orb(in1)
	    do jorb = 1, num_orb(in2)
              Jtemp = Jtemp + Jijsigma(iorb,in1,jorb,in2)
	      icont = icont + 1
	    end do
	   end do
	   Jtemp = Jtemp/(real(icont))
!	   print *, 'J,especie1,especie2',Jtemp,in1,in2
	   do iorb = 1, num_orb(in1)
	      do jorb = 1, num_orb(in2)
	        Jijsigma(iorb,in1,jorb,in2) = Jtemp
	      end do
	   end do
	 end do
       end do



!print *, 'escribimos las Us'
       do in1 = 1, nspecies
        do iorb = 1, num_orb(in1)
!print *, 'escribiendo... iorb,num_orb(in1),in1',iorb,num_orb(in1),in1
!print *,    Uisigma(iorb,1:num_orb(in1),in1)
        end do
       end do
!print *, 'escribimos las Js'
       do in1 = 1, nspecies
         do in2 = in1, nspecies
           do iorb = 1, num_orb(in1)
!print *, 'escribiendo... iorb, in1, 1:num_orb(in2), in2', iorb, in1, num_orb(in2), in2
!print *, Jijsigma(iorb,in1, 1:num_orb(in2),in2)
           end do
         end do
       end do




       end if

       allocate(hf_mat(numorb_max, numorb_max, neigh_max, natoms))

       return
       end subroutine readhartree
