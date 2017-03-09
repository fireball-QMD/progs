! copyright info:
!
!                             @Copyright 2002
!                            Fireball Committee
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


! ceperley_alder.f90
! Program Description
! ===========================================================================
!       This routine compute the ceperley-alder form of the LDA as
! parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981).  The units
! of this program are in atomic units, so the density but be changed to atomic
! units after input and the final answer converted to eV-Angstrom units.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ceperley_alder (rho_in, epsxc, potxc, dpotxc, drvexc)
        use constants_fireball
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: rho_in

! Output
        real, intent (out) :: dpotxc
        real, intent (out) :: drvexc
!        real, intent (out) :: epsx
        real, intent (out) :: epsxc
!        real, intent (out) :: potx
        real, intent (out) :: potxc

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: tolerance = 1.0d-10

! Local Variable Declaration and Description
! ===========================================================================
        real epsx
        real potx
!        real drvexc

        real density
        real densityp
        real densitypp
        real dpotc
        real dpotx
        real depsc
        real ddepsc
        real potc
        real rho
        real rs
        real rsl
        real sqrs

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Convert density to Angstrom units.
        rho = rho_in*abohr**3
        if (rho .le. tolerance) then
         epsx = 0.d0
         potx = 0.d0
         epsxc = 0.d0
         potxc = 0.d0
         dpotxc = 0.0d0
         return
        end if

! Initialize some constants_fireball related to density.
        rs = 0.62035049d0/rho**(1.0d0/3.0d0)
        if (rho .lt. 0.23873241d0) then
         sqrs = sqrt(rs)
         density = 1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs
         epsxc = -0.4581652d0/rs - 0.1423d0/density
         potxc = epsxc - rs*(0.15273333d0/rs**2                               &
     &                + (0.02497128d0/sqrs + 0.01581427d0)/density**2)

         densityp =  1.0529d0/(2.0d0*sqrs) + 0.3334d0
         densitypp = -0.5d0*1.0529d0/(2.0d0*rs*sqrs)
         depsc = 0.1423d0*densityp/(density*density)
         ddepsc = - 2.0d0*0.1423d0*densityp*densityp/density**3              &
     &          + 0.1423d0*densitypp/(density*density)
        else
         rsl = log(rs)
         epsxc = - 0.4581652d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs    &
     &           + 0.002d0*rs*rsl
         potxc = epsxc - rs*(0.15273333d0/rs**2 + 0.01036667d0/rs            &
     &                - 0.003866667d0 + 0.00066667d0*(1.0d0 + rsl))
         depsc = 0.0311d0/rs - 0.0116d0 + 0.0020d0*(rsl + 1.0d0)
         ddepsc = -0.0311d0/(rs*rs) + 0.0020d0/rs
        end if

! Exchange-only energy and potential
        epsx = - 0.7385587664d0*rho**(1.0d0/3.0d0)
        potx = 4.0d0/3.0d0*epsx

! Extended hubbard additions.
        drvexc = (potxc - epsxc)/rho

! Now dpotxc; we compute dpot/dn. We use dpot/dn = 2*dexc/dn + n*d2(exc)/dn2.
! Here dexc/dn = drvexc, and
! dpot/dn = (-2/(9*n*n))*ex + 4*rs/(9*n*n)*dec + rs*rs/(9*n*n)*ddec
! Let dpotc = dpot/dn = 4*rs/(9*n*n)*dec + rs*rs/(9*n*n)*ddec
        dpotc = (4.0d0*rs/(9.0d0*rho*rho))*depsc                             &
     &         + (rs*rs/(9.0d0*rho*rho))*ddepsc
        dpotx = - (2.0d0/(9.0d0*rho*rho))*epsx
        dpotxc = 2.0d0*drvexc + rho*(dpotx + dpotc)

! Convert output to eV-Angstrom units
        epsx = epsx*Hartree
        epsxc = epsxc*Hartree
        potx = potx*Hartree
        potxc = potxc*Hartree
        drvexc = drvexc*Hartree*abohr**3
        dpotxc = dpotxc*Hartree*abohr**3

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end
