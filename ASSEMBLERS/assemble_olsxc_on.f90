! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio State University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! assemble_snxc_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles the two-center exchange-correlation
! (on-site - atom case) for the average density approximation. 
!
! This subroutine could be easily incorporated in assemble_2c.f90 (JOM)
! The double-counting xc (uxcdcc) is also calculated here.
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
        subroutine assemble_olsxc_on (natoms, nprocs, my_proc, iordern,      &
     &                                itheory, uxcdcc)
        use charges
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iordern
        integer, intent (in) :: itheory
        integer, intent (in) :: my_proc
        integer, intent (in) :: natoms
        integer, intent (in) :: nprocs

! Output
        real, intent (out) :: uxcdcc
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer imu
        integer in1, in3
        integer inu
        integer matom
        integer natomsp

        real, dimension (numorb_max, numorb_max) :: bcxcx
        real xc

! Procedure
! ===========================================================================
! Initialize
        vxc = 0.0d0
        if (itheory .eq. 1) vxc_ca = 0.0d0
        uxcdcc = 0.0d0
  
! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
     &                + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if

! Loop over the atoms in the central cell.
        bcxcx  = 0.0d0
!!$omp parallel do private (matom, in1, bcxcx, xc, in3, inu, imu)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         in1 = imass(iatom)

! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
         if (itheory .eq. 1) then
! Dogs
          call build_ca_olsxc_on (in1, iatom, bcxcx, xc)
! double-counting xc correction
!!$omp atomic
          uxcdcc = uxcdcc + xc
          in3 = in1
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            vxc_ca(imu,inu,matom,iatom) =                                       &
     &       vxc_ca(imu,inu,matom,iatom) + bcxcx(imu,inu)
           end do
          end do
          
         else
! Harris + ext Hubbard
          call build_olsxc_on (in1, iatom, bcxcx, xc)
! double-counting xc correction
          uxcdcc = uxcdcc + xc
          in3 = in1
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            vxc(imu,inu,matom,iatom) =                                    &
     &       vxc(imu,inu,matom,iatom) + bcxcx(imu,inu)
           end do
          end do
         end if
        end do ! End loop over iatom.

! Deallocate arrays
! ===========================================================================

! Format Statements
! ===========================================================================
 
        return
        end subroutine assemble_olsxc_on
