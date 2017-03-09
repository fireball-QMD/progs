! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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


! dos.f90
! Program Description
! ===========================================================================
! This routine calculates the Green function for the unit cell in real space
! uses the hamk calculated in SOLVES_DIAG/kspace.f90.
! ===========================================================================
!
! Code written by:
! Cesar Gonzalez Pascual
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madird SPAIN    
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
!
! modified by P. Jelinek 22/05/2006
! ==========================================================================
  subroutine dos(weightk, k_temp)

     use density
     use charges
     use dimensions
     use interactions
     use neighbor_map
     use module_dos
     use configuration

     implicit none

! Argument Declaration and Description
! Input
      real,intent(in)    :: weightk      ! weight of the k
      real, intent (in), dimension (3)         :: k_temp ! k vector

! Hamk: hamiltonian in k space for each k point
! hr_box: hamiltonian in real space in boxes for each atom with its neighbors

! Output

! Local Parameters and Data Declaration
! ===========================================================================
      real, parameter :: PI=3.14159265 

! Local Variable Declaration and Description
! ==========================================================================

      integer            :: norbi
      integer            :: norbj
      integer            :: norbig
      integer            :: norbjg
      integer            :: norb_beg
      integer            :: ie
      integer            :: in1
      integer            :: in2
      integer            :: iorb
      integer            :: jorb
      integer            :: iatom
      integer            :: jatom
      integer            :: katom
      integer            :: ineig
      integer            :: info
      integer            :: lwork

 
      complex            :: energy
      complex            :: a1
      complex            :: a0
 
! Matrix of identity
      real, dimension (numorb_max, numorb_max)          :: Ident

! A bunch of memory to be used in many ways
      complex, dimension(:,:),  allocatable             :: greenk
      complex, dimension(:),  allocatable               :: work
      integer, dimension(:),  allocatable               :: ipiv

! Procedure
! ===========================================================================

! Initialize some things
      write (*,*) 'Wellcome to dos subroutine ...  '

      a1 = (1.0d0, 0.0d0)
      a0 = (0.0d0, 0.0d0)
      lwork = norbitals*norbitals
      allocate(greenk(norbitals,norbitals))
      allocate(work(lwork))
      allocate(ipiv(norbitals))

      greenk = a0         
   
      Ident = a0
   
      do iorb = 1,numorb_max
       Ident(iorb,iorb) = a1 
      end do


! Add imaginary part to energy and shift Fermi level to zero
      energy = ener_beg*a1 + eta*(0.0d0, 1.0d0)
    
! Loop over energy range
      do ie = 1,nener
 
       greenk = 0.0d0
! Add energy step
       norbi = 0
! Loop over atoms
       do iatom = 1,natoms
        norbj = 0
        in1 = imass(iatom)
        if(ie .eq. 1) then
         if(iatom .eq. natom_beg) norb_beg = norbi 
        end if
        do jatom = 1,natoms
         in2 = imass(jatom)
         if(iatom .eq. jatom) then
          greenk(norbi+1:norbi+num_orb(in1),norbj+1:norbj+num_orb(in2)) = &
     &            energy*Ident(1:num_orb(in1),1:num_orb(in2)) -                &
     &            hamk(norbi+1:norbi+num_orb(in1),norbj+1:norbj+num_orb(in2))
         else
          greenk(norbi+1:norbi+num_orb(in1),norbj+1:norbj+num_orb(in2)) = &
     &            - hamk(norbi+1:norbi+num_orb(in1),norbj+1:norbj+num_orb(in2))
         end if
         norbj = norbj+num_orb(in2)
        enddo ! jatom
        norbi = norbi+num_orb(in1)
       enddo ! iatom

!!       call inv(greenk, greenk, norbitals, norbitals)
! LU matrix factorization of general matrix 
       call zgetrf ( norbitals, norbitals, greenk, norbitals, ipiv, info )
       if (info .ne. 0)  then
        write (*,*) ' ***  error in zgetrf  ' 
        write (*,*) ' ***  info = ',info
        stop
       endif
! matrix inversion of general matrix
       call zgetri ( norbitals , greenk, norbitals, ipiv, work, lwork, info)
       if (info .ne. 0) then
        write (*,*) ' ***  error in zgetri  ' 
        write (*,*) ' ***  info = ',info
        stop
       endif


       norbi = 0
       norbig = norb_beg
       do iatom = natom_beg,natom_end
        norbj = 0
        norbjg = norb_beg
        in1 = imass(iatom)
        do jatom = natom_beg,natom_end
         in2 = imass(jatom)
         green(norbi+1:norbi+num_orb(in1),norbj+1:norbj+num_orb(in2),ie) =  &
     &        green(norbi+1:norbi+num_orb(in1),norbj+1:norbj+num_orb(in2),ie)  &
     &        + cexp((0.0,1.0)*(dot_product(k_temp(1:3),(ratom(1:3,iatom)      &
     &        -ratom(1:3,jatom)))))*                                           &
     &        greenk(norbig+1:norbig+num_orb(in1),norbjg+1:norbjg+num_orb(in2))&
     &        *weightk
         norbj = norbj+num_orb(in2)
         norbjg = norbjg+num_orb(in2)
        end do ! jatom
        norbi = norbi+num_orb(in1)
        norbig = norbig+num_orb(in1)
       end do ! iatom
       energy = energy + ener_step*a1
      enddo ! ie

! Deallocates
      deallocate(greenk)
      deallocate(ipiv)
      deallocate(work)

! Return to fireball
      return

100    format (6(2x,f12.6))

 end subroutine dos
