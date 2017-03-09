! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
! University of Utah - James P. Lewis, Chair
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
! Ohio State University - Dave Drabold

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

 
! cepal.f90
! Program Description
! ===========================================================================
!       This routine calculates the lda exchange-correlation energy,
! potential and the derivative of the potential.  The form of the
! functional is that of Ceperley-Alder as parameterized by Perdew-Zunger.
! Phys. Rev. B23, 5048 (1981). Units are Hartree a.u. (but see below)
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine cepal (rh, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: rh

! Output
        real, intent (out) :: dexc
        real, intent (out) :: d2exc
        real, intent (out) :: dmuxc
        real, intent (out) :: exc
        real, intent (out) :: muxc
        real, intent (out) :: d2muxc
 
! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: eps = 1.0d-3
        real, parameter :: delta_rh = 1.0d-6
 
! Local Variable Declaration and Description
! ===========================================================================
        real d2nec
        real d2nex
        real d3nec
        real d3nex
        real dec
        real ddec
        real d2dec
        real den
        real dden
        real d2den
        real d3den
        real ex
        real rho_third
        real rho
        real rhx
        real rs
        real rsl
        real sqrs
        real hartree1
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize things to zero (thero in Spanish).
        exc = 0.0d0
        muxc = 0.0d0
        dexc = 0.0d0
        dmuxc = 0.0d0
        d2muxc = 0.0d0
        d2exc = 0.0d0

        rhx = sqrt(rh*rh + delta_rh)

! Convert to a.u.
        rho = rhx*(abohr**3)

! Find rho^(1/3) 
        rho_third=rho**(1.0e0/3.0e0)

! Effective radius
        rs = 0.62035049d0/rho_third

! Find the energy, potential, and the deerivative of the potential.
        if (rho .lt. 0.23873241d0) then
         sqrs = sqrt(rs)

         den = 1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs       ! Effective density
         exc = -0.4581652d0/rs - 0.1423d0/den
         muxc = exc - rs*(0.15273333d0/rs**2                                 & 
     &              + (0.02497128d0/sqrs + 0.01581427d0)/den**2)

! Stuff for dmuxc
         dden = 1.0529d0/(2.0d0*sqrs) + 0.3334d0
         d2den = (-0.5d0)*1.0529d0/(2.0d0*rs*sqrs)
         d3den = (0.75d0)*1.0529d0/(2.0d0*rs*rs*sqrs)
         dec = 0.1423d0*dden/(den*den)
         ddec = -2.0d0*0.1423d0*dden*dden/(den**3) + 0.1423d0*d2den/(den*den)
         d2dec = 6.0d0*0.1423d0*(dden*3)/(den**4) - 6.0d0*0.1423d0*          &
     &       dden*d2den/(den**3) + 0.1423d0*d3den/(den*den) 

        else

         rsl = log(rs)
         exc = -0.4581652d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs       &
     &        + 0.002d0*rs*rsl
         muxc = exc - rs*(0.15273333d0/rs**2 + 0.01036667d0/rs               & 
     &                    - 0.003866667d0 + 0.00066667d0*(1.0d0 + rsl))

         dec = 0.0311d0/rs - 0.0116d0 + 0.0020d0*(rsl + 1.0d0)
         ddec = -0.0311d0/(rs*rs) + 0.0020d0/rs
         d2dec = 2.0d0*0.0311d0/(rs*rs*rs) - 0.0020d0/(rs*rs)
        end if

! Exchange-only energy and potential
        ex = -0.7385587664d0*rho_third 

! Now find dmuxc 
        dexc = (muxc - exc)/rho
 
! Compute dmu/dn. Use d(mu)/dn = 2*dexc/dn + n*d2(exc)/dn2.
        d2nec = (4.0d0*rs/(9.0d0*rho*rho))*dec + (rs*rs/(9.0d0*rho*rho))*ddec
        d2nex = -(2.0d0/(9.0d0*rho*rho))*ex
        dmuxc = 2.0d0*dexc + rho*(d2nex + d2nec)
        
! Compute d2mu/dn2, using d2(mu)/dn2 = 3*d2(exc)/dn2 + n*d3(exc)/dn3 
        d3nec = (-28.0d0*rs/(27.0d0*rho*rho*rho))*dec + (-4.0d0*rs*rs/       &
     &   (9.0d0*rho*rho*rho))*ddec + (rs*rs*rs/(-27.0d0*rho*rho*rho))*d2dec
        d3nex = (10.0d0/(27.0d0*rho*rho*rho))*ex
        d2muxc = 3.0*(d2nex + d2nec) + rho*(d3nex + d3nec)
        d2exc = d2nex + d2nec

! Convert output to eV (exc and muxc) and to eV*(Angstrom**3) (dexc and dmuxc)
        hartree1 = eq2/abohr
        exc = exc*hartree1
        muxc = muxc*hartree1
        dexc = dexc*hartree1*(abohr)**3
        d2exc = d2exc*hartree1*(abohr)**6
        dmuxc = dmuxc*hartree1*(abohr)**3
        d2muxc = d2muxc*hartree1*(abohr)**6

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================

        return
        end subroutine cepal
