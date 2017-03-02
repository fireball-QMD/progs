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

! assembleG_den.f90
! Program Description
! ===========================================================================
!       Project density on the mesh + double counting correction terms
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
 subroutine assemble_KS_den ()

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use outputs
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input

!Output

 
! Local Parameters and Data Declaration
! ===========================================================================
   real, parameter ::  Debye = 0.208194
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
   integer imu, inu
   integer in1, in2
   integer jatom
   integer mbeta
   integer ineigh
   integer index
   integer ind
   integer i, j, k
   integer imesh
   real x,y,z
   real dip_x
   real dip_y
   real dip_z
   real dip_tot
   real dqi
   real distX
   real distY
   real renorm
   real dens
   real adens
   real psi2

   real, target, dimension (:), allocatable  :: rhotmp

   real, dimension (:), pointer   :: pmat
   character (len=40) filename
   character (len=30) mssg
   
   real dmax
   real, target, allocatable, dimension (:) :: resf
!
!
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


   allocate (resf(0:nrm-1))   
! save previous step
   do index = 0,nrm-1
    resf(index) = drhoG(index)
   enddo
! reset variables
   drhoG = 0.0d0
   

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
          
! assemble density
         dens = 0.0d0
         do imu = 1, num_orb(in1)
          psi2 = psi2m(imu,imesh,iatom)
          dens = dens + rho(imu,imu,ineigh,iatom)*psi2
         enddo ! do imu
! store variation of density at given point
         drhoG(ind) = dens + drhoG(ind)
        enddo ! do imesh
        
! jatom .ne. iatom      
      else
      
! Loop over mesh points shared by iatom and jatom
       do imesh = 1, npsi22m(ineigh,iatom)
! map the point from the extended mesh into the normal mesh
	    ind = ipsi22m(imesh,ineigh,iatom)

! assemble density
        dens = 0.0d0
        do inu = 1, num_orb(in1)
         do imu = 1, num_orb(in2)
          dens = dens + rho(inu,imu,ineigh,iatom)                &
      &        *psi22m(inu,imu,imesh,ineigh,iatom)
         enddo ! do inu
        enddo ! do imu

! store variation of density at given point
        drhoG(ind) = dens + drhoG(ind)
       end do ! do imesh
      endif ! if (Rc_max) 
    end do ! do ineigh 
   end do ! do iatom

! calc the total denstity
   dens = 0.0d0
   do i = 0,nrm-1
    dens = dens + drhoG(i)*dvol
   enddo
   write (*,*) ' -- Total density before renormalization =',dens

! the renormalization factor
   renorm = ztot/dens

! check total density after renormalization
   dens = 0.0d0
   do i = 0,nrm-1
    drhoG(i) = drhoG(i)*renorm
    dens = dens + drhoG(i)*dvol
   enddo
   write (*,*) ' -- Total density after renormalization =',dens

! get drho (rest atomic charge)
   drhoG = drhoG - rhoG0

! evaluate residua of drhoG
   dmax = 0.0d0
   do index = 0, nrm-1
     dmax = max(dmax, abs(drhoG(index)-resf(index)))
   enddo
   write (*,*) ' residua drhoG = ',dmax, dmax*dvol   

! also check if delta density goes to zero
   dens = 0.0d0
   adens = 0.0d0
   do i =0,nrm-1
    adens = adens + abs(drhoG(i)*dvol)
    dens = dens + drhoG(i)*dvol
   enddo
   write (*,*) ' -- Check sum drho should go to zero =',dens
   write (*,*) ' -- |drho| =', adens

! calc dipole with the unit cell
   index = 0
   dip_x = 0.0d0
   dip_y = 0.0d0
   dip_z = 0.0d0
   do k = 0, rm3-1
    do j = 0, rm2-1
     do i = 0, rm1-1
!     dqi = drhoG(index) - rhoG0(index)
      dqi = drhoG(index) 
      x = i*elvec(1,1) + j*elvec(2,1) + k*elvec(3,1)
      y = i*elvec(1,2) + j*elvec(2,2) + k*elvec(3,2)
      z = i*elvec(1,3) + j*elvec(2,3) + k*elvec(3,3)
      dip_x = dip_x + x*dqi
      dip_y = dip_y + y*dqi
      dip_z = dip_z + z*dqi
!     dip_x = dip_x + (x-g0(1))*dqi
!     dip_y = dip_y + (y-g0(2))*dqi
!     dip_z = dip_z + (z-g0(3))*dqi
      index = index + 1
     enddo ! i
    enddo ! j
   enddo ! k
   dip_x = dip_x * dvol
   dip_y = dip_y * dvol
   dip_z = dip_z * dvol
   dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
   write (*,301) dip_x/Debye
   write (*,302) dip_y/Debye
   write (*,303) dip_z/Debye
   write (*,304) dip_tot/Debye

! -----------
! dipole units (charge x distance)
! 1 Debye = 0.208194 eAng = 0.393430 eBohr
! 

! write out files
   if (iwrtxsf .eq. 1) then 

! total density   
    allocate (rhotmp (0:nrm-1))
! get total density
    rhotmp = drhoG + rhoG0
! write out rho into xsf file
    pmat => rhotmp
    filename = 'density.xsf'
    mssg = 'density_3D'
    call writeout_xsf (filename, mssg, pmat) 
    deallocate (rhotmp)  
  
! write out drho into xsf file
    pmat => drhoG
    filename = 'ddensity.xsf'
    mssg = 'ddensity_3D'
    call writeout_xsf (filename, mssg, pmat) 
   endif
   
! deallocate 
    deallocate (resf)    

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6) 
200 format (e14.6) 
301 format (2x,'Dipole_x =',e14.6,'  [D] ') 
302 format (2x,'Dipole_y =',e14.6,'  [D] ') 
303 format (2x,'Dipole_z =',e14.6,'  [D] ') 
304 format (2x,'Dipole_tot =',e14.6,'  [D] ') 

   return
 end subroutine assemble_KS_den

