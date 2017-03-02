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

! get_QMul.f90
! Program Description
! ===========================================================================
!       This routine calculates Mulliken charges from time-dependent
! non-orthogonal basis set.
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
        subroutine get_QMul (ifixcharge)

        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use interactions
        use neighbor_map


        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

        integer, intent (in) :: ifixcharge

! Output

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: halfspin = 0.5d0

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer inu
        integer ineigh
        integer in1
        integer in2
        integer issh
        integer jatom
        integer jneigh
        integer mqn

        real, dimension (numorb_max, natoms) :: QMulliken

! Procedure
! ===========================================================================
! Initialize some things

        write (*,*) '   Compute Mulliken charges '
        write (*,*) ' ****************************************************** '

! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N    C H A R G E S
! ****************************************************************************
! Compute Mulliken charges.
         Qin = 0.0d0
         QMulliken = 0.0d0
         QMulliken_TOT = 0.0d0

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           QMulliken_TOT(iatom) = 0.0d0
           do issh = 1, nssh(in1)
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qin(issh,iatom)
           end do
          end do

         else

          do iatom = 1, natoms
           in1 = imass(iatom)

! Loop over neighbors
           do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            jneigh = neigh_back(iatom,ineigh)
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in2)
              QMulliken(imu,iatom) = QMulliken(imu,iatom)                    &
     &        + 0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) &
     &        + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))*halfspin
             end do
            end do

! End loop over neighbors
           end do

! Finally the imu loop.
           imu = 0
           do issh = 1, nssh(in1)
            do mqn = 1, 2*lssh(issh,in1) + 1
               imu = imu + 1
               Qin(issh,iatom) = Qout(issh,iatom) + QMulliken(imu,iatom)
            end do
               QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qin(issh,iatom)
           end do
! End loop over atoms
          end do
         end if     ! endif of ifixcharges

!         write (*,*) 'QMul =',(QMulliken_TOT(iatom),iatom=1,natoms)

! Format Statements
! ===========================================================================

        return
      end subroutine get_QMul

