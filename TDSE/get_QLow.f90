! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University of Utah - James P. Lewis, Chair
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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! get_QLow.f90
! Program Description
! ===========================================================================
!       This routine calculates Lowdin charges from time-dependent orthogonal
! basis set
!
! ===========================================================================
! Code rewritten by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine get_QLow (ifixcharge)

        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        use tdse

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ifixcharge


! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ielec
        integer ikpoint
        integer imu
        integer inu
        integer ineigh
        integer in1
        integer iorbital
        integer issh
        integer jatom
        integer jneigh
        integer mqn
        integer mmu

        real aux1
        real aux2


! Procedure
! ===========================================================================
! Initialize some things

        aux1 = 0.0d0
        aux2 = 0.0d0

        write (*,*)
        write (*,*) '   Compute Lowdin charges '
        write (*,*) ' ****************************************************** '

! ****************************************************************************
!
!  C O M P U T E    L O W D I N    C H A R G E S
! ****************************************************************************
! Initialize

        if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           QLowdin_TOT(iatom) = 0.0d0
           do issh = 1, nssh(in1)
            QLowdin_TOT(iatom) = QMulliken_TOT(iatom) + Qin(issh,iatom)
           end do
          end do

        else

         Qin = 0.0d0
         QLowdin_TOT = 0.0d0
! Loop over electrons
         do ielec = 1, nelec
! Loop over atoms
          do iatom = 1, natoms
           in1 = imass(iatom)
! Loop over the special k points.
           do ikpoint = 1, nkpoints
            aux1 = weight_k(ikpoint)
! Finally the imu loop
             imu = 0
             do issh = 1, nssh(in1)
              do mqn = 1, 2*lssh(issh,in1) + 1
               imu = imu + 1
               mmu = imu + degelec(iatom)
               aux2 = aux1*(real(psi(mmu,ielec,ikpoint))**2                &
     &                        + aimag(psi(mmu,ielec,ikpoint))**2)
               Qin(issh,iatom) = Qin(issh,iatom) + aux2
               QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux2
              end do ! mqn
             end do ! issh
! End loop over orbitals and kpoints
           end do ! ikpoint
! End loop over atoms
          end do ! iatom
         end do ! ielec
        end if       ! endif ifixcharges

! aux write out of Lowdin charges
        write (*,*) ' QLow =',(QLowdin_TOT(iatom),iatom=1,natoms)


! Format Statements
! ===========================================================================

        return
      end subroutine get_QLow

