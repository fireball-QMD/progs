! copyright info:
!
!                             @Copyright 2009
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

! build_nij.f90
! Program Description
! ===========================================================================
!       This routine calculates the nij charges we need for the HF stuff
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
!
! Program Declaration
! ===========================================================================

     subroutine buildnij(ioccupy_k,foccupy,bmix,Kscf)

     use hartree_fock
     use interactions
     use kpoints
     use constants_fireball
     use density
     use configuration
     implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

     integer, intent(in), dimension (norbitals, nkpoints) :: ioccupy_k
     real, intent(in), dimension (norbitals, nkpoints) :: foccupy
     real, intent(in) :: bmix
     integer, intent(in) :: Kscf





! Local Parameters and Data Declaration
! ===========================================================================

     real, dimension(numorb_max, natoms, numorb_max, natoms) :: nij_old


! Local Variable Declaration and Description
! ===========================================================================

     complex   :: auxi,aux3
     real      :: aux1,aux2
     integer   :: iatom, jatom
     integer   :: in1,in2
     integer   :: ikpoint
     integer   :: iorbital
     integer   :: imu,jmu
     integer   :: issh,jssh
     integer   :: mmu,lmu
     integer   :: mqn,lqn
     integer   :: iorb,jorb

     auxi=(0.0d0,1.0d0)


! We save nij of the previous scf step for mixing between new and old nij

  if ( Kscf .eq. 1 ) then
    nij_old = 0.0
  else
    nij_old = nij
  end if

! First of all we initialize nij to zero
! Otherwise, nij will be nij_new + nij_old and that's not what we want

nij = 0.0


     do iatom = natomhf_beg, natomhf_end
     in1 = imass(iatom)

! Loop over the special k points.
      do ikpoint = 1, nkpoints
       aux1 = weight_k(ikpoint)*spin
       do iorbital = 1, norbitals
        if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
         aux2 = aux1*foccupy(iorbital,ikpoint)

! Finally the imu loop.
        imu = 0
        do issh = 1, nssh(in1)
         do mqn = 1, 2*lssh(issh,in1) + 1
           imu = imu + 1
           mmu = imu + degelec(iatom)

           do jatom = natomhf_beg, natomhf_end
             in2 = imass(jatom)
             jmu = 0
             do jssh = 1, nssh(in2)
              do lqn = 1, 2*lssh(jssh,in2) + 1
               jmu = jmu + 1
               lmu = jmu + degelec(jatom)
               aux3 = aux2*(blowre(mmu,iorbital,ikpoint)+ auxi*blowim(mmu,iorbital,ikpoint))* &
    &                      (blowre(lmu,iorbital,ikpoint)- auxi*blowim(lmu,iorbital,ikpoint))
!               aux3 = aux2*(bbnkre(mmu,iorbital,ikpoint)+ auxi*bbnkim(mmu,iorbital,ikpoint))* &
!    &                      (bbnkre(lmu,iorbital,ikpoint)- auxi*bbnkim(lmu,iorbital,ikpoint))
               nij(imu,iatom,jmu,jatom) = nij(imu,iatom,jmu,jatom) + aux3
              end do
             end do
!                Qout(issh,iatom) = Qout(issh,iatom) + aux3
!                QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux3
           end do
         end do
        end do
        end if

! End loop over orbitals and kpoints
       end do
      end do

! End loop over atoms
     end do

! Now we will write nij for restarting purposes


     open(271, file = './nijmatrix', status='unknown')

     if ( Kscf .eq. 1 ) then
     do iatom = natomhf_beg, natomhf_end
       do jatom = natomhf_beg, natomhf_end
         do iorb = 1, numorb_max
	   do jorb = 1, numorb_max
!	     nij(iorb,iatom,jorb,jatom) = bmix*nij(iorb,iatom,jorb,jatom) + (1-bmix)*nij_old(iorb,iatom,jorb,jatom)
	     write(271,*) nij(iorb,iatom,jorb,jatom)  ! we want not to mix the input and output nij (because nij_old is zero)
	   end do
	 end do
       end do
     end do
     else
     do iatom = natomhf_beg, natomhf_end
       do jatom = natomhf_beg, natomhf_end
         do iorb = 1, numorb_max
	   do jorb = 1, numorb_max
	     nij(iorb,iatom,jorb,jatom) = bmix*nij(iorb,iatom,jorb,jatom) + (1-bmix)*nij_old(iorb,iatom,jorb,jatom)
	     write(271,*) nij(iorb,iatom,jorb,jatom)  ! we want to mix the input and output nij (as in CHARGES)
	   end do
	 end do
       end do
     end do
     end if

     close(271)

     return

     end subroutine
