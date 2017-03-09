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


! den2mesh.f90
! Program Description
! ===========================================================================
!       Project density on the mesh.
!
! ===========================================================================
! Code written by:
! ===========================================================================
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! Program Declaration
! ===========================================================================
 subroutine assemble_KS_den0 ()

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use outputs
   use constants_fireball
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input

!Output

 
! Local Parameters and Data Declaration
! ===========================================================================
   interface 
    subroutine writeout_xsf (xsfname, message, aa)
     real, dimension (:), pointer, intent (in) :: aa
     character (len=40) xsfname
     character (len=30) message
    end
  end interface
  
! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer imu
   integer in1
   integer index
   integer index0
   integer ind
   integer i, j, k
   integer imesh
   integer hit
   integer issh

   real qtot
   real renorm
   real dens
   real psi2


   real, dimension (:), pointer   :: pmat
   character (len=40) filename
   character (len=30) mssg

!   real, dimension (3,natoms)     :: ratom2g
   real, dimension (3,3)          :: aamat
 
!                 + X0 (iatom)
!                / \     u1X = g1 - X0
! uX0 = X0- g0  /   \
!              /     + g1 (nearest grid point to iatom)
!             /
!            /
!           +
!          g0 (origin of grid coords) 
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    X0g0
!    u1X = g1 - X0
!    
! Procedure
! ===========================================================================

! reset variables
   rhoG0 = 0.0d0
   
! integration checking 
   renorm = 0.0d0
   hit = 0
   
! Loop over atoms 
   do iatom = 1, natoms

    in1 = imass(iatom)

! Loop over points in the atomic mesh gP
    do imesh = 1, nam 

! get index within the mesh
     ind = ipsi2m (imesh,iatom)
          
! assemble density
     dens = 0.0d0
     do imu = 1, num_orb(in1)
      psi2 = psi2m(imu,imesh,iatom)
      dens = dens + rhoA(imu,iatom)*psi2
      renorm = renorm + dvol
      hit = hit + 1
     enddo ! do imu
! store variation of density at given point
     rhoG0(ind) = dens + rhoG0(ind)
     write (302,300) ind, rhoG0(ind)
    end do ! do imesh
   end do ! do iatom


! test for hydrogen s-orbital the quality of the integration; 
! the ratio it should go to one with higher Ecut
   write (*,*) 'Rc_max =',Rc_max
   write (*,*) 'vol    =',dvol
   write (*,*) '# hit  =',hit
   dens = 4.0d0*3.141592653589793238462643*Rc_max**3/3.0d0
   write (*,'(a,5f14.7)') 'SPHERE =',dvol,renorm, dens, renorm-dens,renorm/dens

! calc the total denstity
   dens = 0.0d0
   do i = 0,nrm-1
    dens = dens + rhoG0(i)*dvol
   enddo
   write (*,*) ' -- Total atomic density before renormalization =',dens

! sum total charge of the system
   qtot = 0.0d0
   do iatom = 1, natoms
     in1 = imass(iatom)
     do issh = 1, nssh(in1)
       qtot = qtot + Qneutral(issh,in1)
     end do
   end do

! the renormalization factor
   renorm = qtot/dens

! check total density after renormalization
   dens = 0.0d0
   do i = 0,nrm-1
    rhoG0(i) = rhoG0(i)*renorm
    dens = dens + rhoG0(i)*dvol
   enddo
   write (*,*) ' -- Total atomic density after renormalization =',dens
   
! write out files
   if (iwrtxsf .eq. 1) then 
! write out rho into xsf file
    pmat => rhoG0
    filename = 'density_atm.xsf'
    mssg = 'density_3D'
    call writeout_xsf (filename, mssg, pmat) 
   endif

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6) 
200 format (e14.6) 
300 format (i8,f16.8)

   return
 end subroutine assemble_KS_den0

