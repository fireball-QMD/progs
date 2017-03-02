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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! assembleG_mat.f90
! Program Description
! ===========================================================================
!  Assemble matrix elements of orbitals.
!
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
 subroutine assemble_KS_mat ()

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

 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer imu, inu
   integer in1, in2
   integer jatom
   integer mbeta
   integer ineigh
   integer index
   integer ind
   integer i, j, k
   integer imesh


   real psi12
   real dens


!   real, dimension (numorb_max,numorb_max)   :: integ1
!   real, dimension (numorb_max,numorb_max)   :: integ2
   real ex, mux, exc, muxc, dexc, d2exc, dmuxc, d2muxc

!                 + X0 (iatom)
!                /|\     u1X = g1 - X0
! uX0 = X0- g0  / | \
!              /  |  + g1 (nearest grid point to iatom)
!             /   | /
!            /    |/
!           /     + Y0
!          +      
!          g0 (origin of grid coords) 
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    uX0 = X0 - g0
!    u1X = g1 - X0
!    r21 = Y0 - X0
!    u1Y = g1 - Y0 = g1 - Y0 - X0 + X0 = u1X - r21
!    
! Procedure
! ===========================================================================

! reset variables
   vca = 0.0d0
   vxc = 0.0d0
   vxc_ca = 0.0d0
   

! Loop over atoms 
   do iatom = 1, natoms
! make a copy of the elem grid lattice vector
    in1 = imass(iatom)

! Loop over the neighbors
    do ineigh = 1, neighn(iatom)

      jatom = neigh_j(ineigh,iatom)
      mbeta = neigh_b(ineigh,iatom)
      in2 = imass(jatom)

! case self-atom      
      if (iatom .eq. jatom .and. mbeta .eq. 0) then
 ! Loop over mesh points shared by iatom and jatom     
        do imesh = 1, nam

! get index within the mesh
         ind = ipsi2m (imesh,iatom)

! get the total density at the point 
         dens = drhoG(ind) + rhoG0(ind)
! calc XC potential 
         call ceperley_alder (dens, exc, muxc, dmuxc, dexc)
         vxcG(ind) = muxc 
         
! assemble matrices
         do inu = 1, num_orb(in1)
          psi12 = psi2m(inu,imesh,iatom)
! vca  dV_H
          vca(inu,inu,ineigh,iatom) =  vca(inu,inu,ineigh,iatom)        &
    &           + vcaG(ind)*psi12*dvol
! vxc  Vxc
          vxc(inu,inu,ineigh,iatom) =  vxc(inu,inu,ineigh,iatom)        &
    &           + vxcG(ind)*psi12*dvol
! rho
          vxc_ca(inu,inu,ineigh,iatom) =  vxc_ca(inu,inu,ineigh,iatom)  &
    &            + psi12*dvol
    
         enddo ! do inu                 
        enddo ! do imesh
! jatom .ne. iatom      
      else
      
! Loop over mesh points shared by iatom and jatom
       do imesh = 1, npsi22m(ineigh,iatom)
! map the point from the extended mesh into the normal mesh
	    ind = ipsi22m(imesh,ineigh,iatom)

! assemble matrices
        do inu = 1, num_orb(in1)
         do imu = 1, num_orb(in2)
           psi12 = psi22m(inu,imu,imesh,ineigh,iatom)
! vca  dV_H
           vca(inu,imu,ineigh,iatom) =  vca(inu,imu,ineigh,iatom)        &
    &           + vcaG(ind)*psi12*dvol
! vxc  Vxc
           vxc(inu,imu,ineigh,iatom) =  vxc(inu,imu,ineigh,iatom)        &
    &           + vxcG(ind)*psi12*dvol
! rho
           vxc_ca(inu,imu,ineigh,iatom) =  vxc_ca(inu,imu,ineigh,iatom)  &
    &            + psi12*dvol
         enddo ! do imu
        enddo ! do inu
       end do ! do imesh
      endif ! if (iatom .eq. jatom) 
    end do ! do ineigh 
   end do ! do iatom
      
!vca=0.0d0
!vxc=0.0d0

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6) 
101 format (<numorb_max>f14.6) 
300 format (i8,f16.8)
   return
 end subroutine assemble_KS_mat

