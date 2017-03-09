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


! writeout_atom.f90
! Program Description
! ===========================================================================
! The soubroutine gives us the files Atomo_i for the STM calculation           
! ===========================================================================
! Code written by:
! Cesar Gonzalez Pascual (modified by P. Jelinek)
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madird SPAIN    
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
!
! ==========================================================================
!
! Program Declaration
! ==========================================================================
 subroutine writeout_atom ()  

     use configuration
     use module_dos
     use neighbor_map
     use interactions
     use kpoints
     use charges

     implicit none

! Argument Declaration and Description
! ===========================================================================
! input
     

! input and output

! output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
     real,allocatable                     :: dist_hop(:,:)
     real,dimension(3)                    :: vector
     real                                 :: distance
     integer                              :: n_tipos
     character(4)                         :: ind_aux
     integer                              :: iatom
     integer                              :: jatom
     integer                              :: iorb
     integer                              :: in1,in2
     integer                              :: ineigh
     integer                              :: ineigh_i
     integer                              :: itipo
     integer                              :: jtipo
     integer                              :: mbeta
     integer                              :: ikpts
     integer                              :: ispec
     integer                              :: issh
! Procedure
! ===========================================================================

     n_tipos = maxval(imass(1:natoms))
     allocate(dist_hop(0:n_tipos,0:n_tipos))
     dist_hop = 0.0d0

!!! write the atomo_i files for the STM
     do iatom = 1,natoms
       in1 = imass(iatom)

!       if (iatom .le. 9) write (ind_aux,'(i1)')iatom
!       if (iatom .gt. 9 .and. iatom .lt. 100) write (ind_aux,'(i2)') iatom
!       if (iatom .gt. 99 .and. iatom .lt. 1000) write (ind_aux,'(i3)') iatom
!       if (iatom .gt. 999 .and. iatom .lt. 10000) write (ind_aux,'(i4)') iatom
       if (iatom .gt. 9999) then
         write (*,*) ' too many atoms, must stop here!'
         stop
       endif
       write(ind_aux,'(i4.4)') iatom

       open (unit=10, file='Atomo_'//ind_aux, status='unknown')
       write (10,*) 'Niveles y hoppings intraatomicos del atomo', iatom
       write (10,*) num_orb(in1)
       write (10,*) efermi
       do iorb = 1,num_orb(in1)     
          write (10,'(15f16.8)') Hr_box(iorb,1:num_orb(in1),iatom,0)
       end do               

       write (10,*) 'Hoppings interatomicos'
       write (10,*) 'vecinos:', neighn_tot(iatom)-1
       ineigh_i = 0
       do ineigh = 1,neighn_tot(iatom)

         jatom = neighj_tot(ineigh,iatom)
         mbeta = neighb_tot(ineigh,iatom)
         if (jatom .eq. iatom .and. mbeta .eq. 0) then
         else
          ineigh_i = ineigh_i + 1
          in2 = imass(jatom)
          vector(1:3) = (xl(1:3,neighb_tot(ineigh,iatom))+ratom(1:3,jatom)   &
                      - ratom(1:3,iatom)) / lattice
          distance = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)

!          write (10,'(a8,2x,i3,8x,3f10.6)')'vecino:',ineigh_i, (vector(itipo),itipo=1,3)
          write (10,300) 'vecino:',jatom,mbeta,vector(1:3)
          write (10,*) num_orb(in2)            
          do iorb = 1,num_orb(in1)              
             write (10,'(15f16.8)') Hr_box(iorb,1:num_orb(in2),iatom,ineigh)
          end do                                
          if(distance .gt. dist_hop(in1,in2)) then
              dist_hop(in1,in2) = distance + 0.01d0
              dist_hop(in2,in1) = distance + 0.01d0
          end if
         end if
       end do                                  
       close(10)
     end do                          
    
!!! write the file hoppings (the cuttoff radius for each pair of types of atoms)

     open (4, file='hoppings.inp')
     write (4,*) '0.0'
     do itipo = 0, n_tipos
        write (4,*) '0.0'
        do jtipo = 0,itipo
           write (4,*) dist_hop(itipo,jtipo)
        enddo ! jtipo 
     enddo ! itipo
     close(4)

!!! write the file struc.inp (structural properties) for the STM calculation
     open (unit = 11, file = 'struc.inp', status = 'unknown')
! write number of atoms
     write (11,*) natoms
! write set of atoms to be considered for tunneling
     write (11,*) 1, natoms
! write max. number of neighbors
     write (11,*) num_neig_maxtot
! write atomic coordinates and atomic type
     if (ishiftO .eq. 1) then 
        do iatom = 1, natoms
!     write (11,'(3f13.6,X,2I5)') ratom(1:3,iatom)/lattice,imass(iatom)+2,iatom
           write (11,100) (ratom(1:3,iatom)-shifter(1:3))/lattice,imass(iatom)
        end do
     else 
        do iatom = 1, natoms
!     write (11,'(3f13.6,X,2I5)') ratom(1:3,iatom)/lattice,imass(iatom)+2,iatom
           write (11,100) ratom(1:3,iatom)/lattice,imass(iatom)
        end do
     endif
! write number of orbitals per each specie
     write (11,*) (num_orb(itipo),itipo=1,nspecies) 
     write (11,*) (nssh(issh), issh = 1, nspecies)
     do ispec = 1, nspecies
      write (11,*) (lssh(issh,ispec), issh = 1, nssh(ispec))
     enddo
     write (11,*) '0'
     write (11,*) '2 2'
! write lattice vector
     write (11,200) a1vec(1:3)/lattice 
     write (11,200) a2vec(1:3)/lattice 
     write (11,200) a3vec(1:3)/lattice 
     close (11)

! write out the kpts file in lattice parameter units
     open(11, file='KRCNST.dat')
     write(11,*) nkpoints
     do ikpts = 1,nkpoints
        write (11,400) special_k(1:3,ikpts)*lattice,weight_k(ikpts),ikpts
     end do

     close (11)

! Format Statements
! ===========================================================================
100     format (3f18.8, 2x, i5)
!100     format (3f18.8, 2x, 2i5)
200     format (3f18.8)
300     format (a8,2x,i4,2x,i4,6x,3f16.8)
400     format (4f16.6, i5)    
     return

   end subroutine writeout_atom
