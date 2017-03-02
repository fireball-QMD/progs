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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! hamilt_atom.f90
! Program Description
! ===========================================================================
! This subroutine writes down the hamiltonian in real space for each atom
! in format used by STM program.
! ===========================================================================
! Code written by:
! Cesar Gonzalez Pascual
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madird SPAIN    
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
! ==========================================================================

  subroutine hamilt_atom(weightk, k_temp)

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
  real,intent(in)           :: weightk      ! weight of the k
  real, intent (in), dimension (3)         :: k_temp ! k vector

! hr_box: hamiltonian in real space in boxes for each atom with its neighbours
! Local Variable Declaration and Description
! ==========================================================================
  integer            :: norbi
  integer            :: norbj
  integer            :: in1
  integer            :: in2
  integer            :: iatom
  integer            :: jatom
  integer            :: katom
  integer            :: ineigh
  integer            :: mbeta 

  norbi = 0
  do iatom = 1,natoms
    norbj = 0
    in1 = imass(iatom)
    do ineigh = 1,neighn_tot(iatom)

      jatom = neighj_tot(ineigh,iatom)
      mbeta = neighb_tot(ineigh,iatom)

      if (jatom .eq. iatom .and. mbeta .eq. 0) then
        Hr_box(1:num_orb(in1),1:num_orb(in1),iatom,0) =                      &
    &       Hr_box(1:num_orb(in1),1:num_orb(in1),iatom,0) +                  &
!    &      real(hamk(norbi+1:norbi+num_orb(in1),norbi+1:norbi+num_orb(in1))  &
!    &      *weightk)
    &       real((cexp((0.,1.)*(dot_product(k_temp(1:3),                   &
    &       (ratom(1:3,iatom) - ratom(1:3,jatom))))))*     &
    &       hamk(norbi+1:norbi+num_orb(in1),norbi+1:norbi+num_orb(in1))    &
    &       *weightk)
      else
        norbj=0
        do katom = 1, jatom-1
          in2 = imass(katom)
          norbj = norbj+num_orb(in2)
        end do
        in2 = imass(jatom)
        Hr_box(1:num_orb(in1),1:num_orb(in2),iatom,ineigh) =                 &
    &       Hr_box(1:num_orb(in1),1:num_orb(in2),iatom,ineigh)+              &
    &         real((cexp((0.,1.)*(dot_product(k_temp(1:3),                   &
    &         (ratom(1:3,iatom) - xl(1:3,mbeta) - ratom(1:3,jatom))))))*     &
    &         hamk(norbi+1:norbi+num_orb(in1),norbj+1:norbj+num_orb(in2))    &
    &         *weightk)
      end if
    end do !
    norbi = norbi + num_orb(in1)
  end do

  return

end subroutine hamilt_atom
