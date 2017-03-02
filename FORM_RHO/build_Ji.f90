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

! assemble_hartree.f90
! Program Description
! ===========================================================================
!       This routine calculates the ji parameters we need for the HF stuff
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
     subroutine build_ji

     use hartree_fock
     use interactions
     use configuration
     implicit none

     real   :: denominador
     real   :: numerador
     real,dimension(3)  :: r21
     real   :: r21mod
     real   :: nij2
     integer  :: iatom, jatom
     integer  :: iorb, jorb
     integer  :: in1, in2


     do iatom = natomhf_beg, natomhf_end
       in1 = imass(iatom)

! Finally the imu loop.
!       imu = 0
       do iorb = 1, num_orb(in1)
!       do mqn = 1, 2*lssh(issh,in1) + 1
!        imu = imu + 1
!        mmu = imu + degelec(iatom)

         numerador = 0.0d0
         denominador = 0.0d0
!         do jatom = 1, natoms
         do jatom = natomhf_beg, natomhf_end
          in2 = imass(jatom)
!         jmu = 0
          do jorb = 1, num_orb(in2)
!           do lqn = 1, 2*lssh(jssh,in2) + 1
!           jmu = jmu + 1
!           lmu = jmu + degelec(jatom)
            nij2 = real(nij(iorb,iatom,jorb,jatom)*Conjg(nij(iorb,iatom,jorb,jatom)))
            denominador = denominador + nij2
            if(iatom.eq.jatom) then
              numerador = numerador +         &
      &            Uisigma(iorb,jorb,in1)*nij2
            else
              r21(:) = ratom(:,jatom) - ratom(:,iatom)
              r21mod = sqrt(r21(1)**2+r21(2)**2+r21(3)**2)
              numerador = numerador + nij2*            &
      &             Jijsigma(iorb,in1,jorb,in2)/r21mod
!      &             Jijsigma(iorb,in1,jorb,in2)
            end if
           end do
          end do
          Jialpha(iorb,iatom) = numerador/denominador ! check
        end do

! End loop over atoms
     end do

     return

     end subroutine
