! Copyright info:
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


! assembleG_dcc.f90
! Program Description
! ===========================================================================
! it calcs double counting corrections of the XC and the (part) Hartree terms
! 1. XC term:  int (exc[n]-muxc[n])n(r) dr
! 2. Hartree term: -1/2 int dVH*n(r) dr
! 3. Hartree term: -1/2 int Vna*dn(r) dr
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine assemble_KS_dcc (uxcdcc, uhdcc)

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
   
!Output
   real, intent (inout)         :: uxcdcc
   real, intent (inout)         :: uhdcc

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
  integer i
  integer iatom
  integer jatom
  integer matom
  integer imu
  integer inu
  integer in1
  integer in2
  integer ineigh
  real dens 
  real exc
  real dexc 
  real muxc
  real dmuxc
  real vnadn
  real dvhn0
  real dvhdn
  real vnan
  real dvhn

! Procedure
! ===========================================================================

! reset variables
   uxcdcc = 0.0d0
   uhdcc = 0.0d0
   vnadn = 0.0d0
   vnan  = 0.0d0
   dvhn0 = 0.0d0
   dvhdn = 0.0d0
   dvhn  = 0.0d0
   
! Loop over atoms 
   do i = 0,nrm-1
    dens = drhoG(i) + rhoG0(i) 
! calc XC potential 
    call ceperley_alder (dens, exc, muxc, dmuxc, dexc)
    uxcdcc = uxcdcc + (exc-muxc)*dens
! Hartree part
    uhdcc = uhdcc + dens*vcaG(i) 
    vnadn = vnadn + drhoG(i)*vnaG(i)
    vnan  = vnan  + dens*vnaG(i)
    dvhn0 = dvhn0 + vcaG(i)*rhoG0(i)
    dvhdn = dvhdn + vcaG(i)*drhoG(i)  
    dvhn  = dvhn  + vcaG(i)*dens 
   enddo  

! multiply by integraation factor
   uxcdcc = uxcdcc*dvol 
   dvhdn = -0.5d0*dvhdn*dvol  
   dvhn0 = -1.0d0*dvhn0*dvol
   uhdcc = dvhdn + dvhn0
   
   write (*,100)
   write (*,200) uxcdcc
   write (*,310) dvhdn
   write (*,315) dvhn0

!   uhdcc = uhdcc + comp_VNA
   write (*,300) uhdcc
   write (*,*)  'VNA*dn =',-0.5d0*vnadn*dvol
   write (*,*)  'VNA*n  =',-0.5d0*vnan*dvol
   write (*,*)  'dVH*n  =',-0.5d0*dvhn*dvol

! Format Statements
! ===========================================================================
100 format ( ' Double counting correction terms from the Grid:') 
200 format ( ' uXCdcc =',f14.6) 
300 format ( ' uHdcc  =',f14.6) 
310 format ( ' dVH*dn  =',f14.6) 
315 format ( ' dVH*n0  =',f14.6) 
 
   return
 end subroutine assemble_KS_dcc

