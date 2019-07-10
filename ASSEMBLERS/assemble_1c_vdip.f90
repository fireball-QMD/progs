! copyright info:
!
!                             @Copyright 2002
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

 
! assemble_1c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the one-center exchange-correlation
! interactions. The results are stored in vxc_1c and etotxc_1c.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_1c_vdip (natoms, itheory, iforce)
        use charges
        use dimensions
        use interactions
        use neighbor_map
        use energy
        use density, only: rho
        use constants_fireball, only: eq2
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: iforce
        integer, intent(in) :: itheory
        integer, intent(in) :: natoms
 
! Output
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu, inu, ialpha, ibeta
        integer in1
        integer matom
        integer ixc
        integer indx
 
        real dccexc_1c
        real exc_1c
        real muexc_1c
        real I
        real Integral
        real, dimension (numorb_max,numorb_max) :: mu1xc
 
! Procedure
! ===========================================================================

!..........................NEW TRY

 
       Vdip_1c = 0.0d0
       dc_v_intra_dip_1c = 0.0d0
         do iatom = 1,natoms

          in1 = imass(iatom)
          matom = neigh_self(iatom)

          do indx = 1,Nlines_vdip1c(in1)
            !Here imu,inu,ialpha,ibeta are assumed to be pairwise different
            imu    = muR(indx,in1)
            inu    = nuR(indx,in1)
            ialpha = alphaR(indx,in1)
            ibeta  = betaR(indx,in1)
            I      = IR(indx,in1)
            if (imu .ne. inu) then
            !eq2?
            Vdip_1c(imu,inu,iatom) = (rho(ialpha,ibeta,matom,iatom)+rho(ibeta,ialpha,matom,iatom))*I*eq2
            Vdip_1c(inu,imu,iatom) = (rho(ialpha,ibeta,matom,iatom)+rho(ibeta,ialpha,matom,iatom))*I*eq2
             
            dc_v_intra_dip_1c = dc_v_intra_dip_1c- &
            & 0.5*(rho(imu,inu,matom,iatom)*Vdip_1c(imu,inu,iatom)+ &
            &  rho(inu,imu,matom,iatom)*Vdip_1c(inu,imu,iatom))
            end if
            !There's another DC term remaining..!!  --> -xi^2 xi xi
            !arising from rho_0*rho_dip   

          end do !end do indx


             end do !end do iatom

! Format Statements
! ===========================================================================
!100      format (2x, ' <|vxc_1c|> = ', 9f8.3) 
        return
      end subroutine assemble_1c_vdip
 
