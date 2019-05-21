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


! Program Description
! ===========================================================================
!       This routine assembles the two-center exchange-correlation
! (on-site - atom case) for the average density approximation. 
!
! This subroutine could be easily incorporated in assemble_2c.f90 (JOM)
! The double-counting xc (uxcdcc_dc) is also calculated here.
!
! ===========================================================================
! Code written by:
!
! ===========================================================================
        subroutine assemble_zw_on_na (natoms, nprocs, my_proc, iordern,      &
     &                                itheory)
        use charges
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use energy, only : uxcdcc_zw
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
        uxcdcc_zw = 0.0d0   !this quantity is initialized here
  
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
          call build_zw_on_na (in1, iatom, bcxcx, xc)
! double-counting xc correction
          uxcdcc_zw = uxcdcc_zw + xc
          in3 = in1
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            vxc(imu,inu,matom,iatom) =                                    &
     &       vxc(imu,inu,matom,iatom) + bcxcx(imu,inu)
           end do
          end do
        end do ! End loop over iatom.

! Deallocate arrays
! ===========================================================================
! Format Statements
! ===========================================================================
 
        return
        end subroutine assemble_zw_on_na
