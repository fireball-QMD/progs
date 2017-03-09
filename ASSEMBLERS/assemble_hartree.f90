! copyright info:
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

! assemble_hartree.f90
! Program Description
! ===========================================================================
!       This routine assembles the Hartree-Fock potential to correctly
!	calculate the molecular gap
!
! ===========================================================================
! Code written by:
! Enrique Abad Gonzalez
! Dpto de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! Phone: +34914978648
! email: enrique.abad@uam.es
!
! ===========================================================================
        subroutine assemble_hartree ()
        use charges
        use density
        use dimensions
        use configuration
        use forces
        use interactions
        use neighbor_map
        use hartree_fock
!	use OO_variable
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Output
!        real, dimension (numorb_max, numorb_max, neigh_max, natoms), 		&
!     &		intent (out) :: hf_mat

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom
	integer ineigh
	integer mbeta
        integer in1, in2
        integer icount
	integer iorb, jorb
        real, dimension(3)  :: r1, r2, r21
        real    r21mod
        real, dimension(:,:), allocatable ::    alfa
        complex, dimension(:,:), allocatable :: nitimeshole
!	real, dimension(:,:,:,:), allocatable :: hf_aux

! Procedure
! ===========================================================================
! Initialize
!	allocate(hf_aux(numorb_max, numorb_max, neigh_max, natoms))
        hf_mat = 0.0
!	hf_aux = 0.0


! We are going to calculate now the fraction of hole in the molecule (alfa)

allocate(nitimeshole(numorb_max,natoms))
allocate(alfa(numorb_max,natoms))
nitimeshole = 0.0
alfa = 0.0
do iatom = natomhf_beg, natomhf_end
  in1 = imass(iatom)
  do iorb = 1, num_orb(in1)
    do jatom = natomhf_beg, natomhf_end
      in2 = imass(jatom)
      do jorb = 1, num_orb(in2)
        if ( (iatom.ne.jatom).or.(iorb.ne.jorb) ) then
        nitimeshole(iorb,iatom) = nitimeshole(iorb,iatom) + 0.25*(nij(iorb,iatom,jorb,jatom)*Conjg(nij(iorb,iatom,jorb,jatom)))
        end if
      end do
    end do
    if ( nij(iorb,iatom,iorb,iatom) .ne. 0.0  ) then
    alfa(iorb,iatom) = nitimeshole(iorb,iatom)/(0.25*nij(iorb,iatom,iorb,iatom)*(2.0-nij(iorb,iatom,iorb,iatom)))
    end if
    print *, 'iorb,iatom,ni(1-ni),sum_j nij,alfa',iorb,iatom,&
                   &0.25*nij(iorb,iatom,iorb,iatom)*(2.0-nij(iorb,iatom,iorb,iatom)),&
                   &nitimeshole(iorb,iatom),alfa(iorb,iatom)
  end do
end do

! alfa calculated

	do iatom = natomhf_beg, natomhf_end ! Loop over the atoms in the central cell.
          r1(:) = ratom(:,iatom)
          in1 = imass(iatom)

          do ineigh = 1, neighn(iatom)       ! <==== loop over i's neighbors
            jatom = neigh_j(ineigh,iatom)
	    in2 = imass(jatom)
            mbeta = neighb_tot(ineigh,iatom)

            icount = 0
            if(iatom.le.natomhf_end.and.iatom.ge.natomhf_beg) icount = icount+1
            if(jatom.le.natomhf_end.and.jatom.ge.natomhf_beg) icount = icount+1
            if (icount.eq.2) then
              if((mbeta.eq.0).and.(iatom.eq.jatom)) then
	        do iorb = 1, num_orb(in1)
	          do jorb = 1, num_orb(in2)
! First of all we calculate the HF matrix in the Lowdin basis
	            if (iorb.ne.jorb) then
                    hf_mat(iorb,jorb,ineigh,iatom) =         &
      &          -  Uisigma(iorb,jorb,in1)*0.5*nij(iorb,iatom,jorb,jatom)  ! Yannick
!	&               0.0
	            else
		    hf_mat(iorb,jorb,ineigh,iatom) =         &
      &     +  betha*alfa(iorb,iatom)*Jialpha(iorb,iatom)*(0.5-(0.5*nij(iorb,iatom,iorb,iatom)))  ! This 0.5 nij is due to the spin
                    end if
	          end do
	        end do
              else
                r2(:) = ratom(:,jatom) + xl(:,mbeta)
                r21(:) = r2(:) - r1(:)
                r21mod = sqrt(r21(1)**2+r21(2)**2+r21(3)**2)
	        do iorb = 1, num_orb(in1)
	          do jorb = 1, num_orb(in2)
! First of all we calculate the HF matrix in the Lowdin basis
                    hf_mat(iorb,jorb,ineigh,iatom) =         &
!      &          -  Jijsigma(iorb,in1,jorb,in2)*0.5*nij(iorb,iatom,jorb,jatom)  ! This 0.5 nij is due to the spin
      &       -  betha*(Jijsigma(iorb,in1,jorb,in2)/r21mod)*0.5*nij(iorb,iatom,jorb,jatom)  ! This 0.5 nij is due to the spin
	          end do
	        end do
              end if
! Now we change to the atomic basis
! Not neccesary by now
            end if
          end do
        end do ! End loop over iatom.

!        hf_mat = hf_aux






! Deallocate arrays
! ===========================================================================

!  deallocate(hf_aux)

! Format Statements
! ===========================================================================
        return
        end subroutine assemble_hartree
