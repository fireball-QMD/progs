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
 subroutine assemble_KS_den (icluster)

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
   integer, intent (in) :: icluster

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
   integer index0
   integer ind
   integer i, j, k
   integer i0, j0, k0
   integer lmu
   integer issh
   integer l
   integer imesh
   integer, dimension (3) :: ipiv
   integer, dimension (3) :: nr

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


   real, dimension (3) :: r1
   real, dimension (3) :: r2
   real, dimension (3) :: r21
   real, dimension (3) :: u
   real, dimension (3) :: dXr
   real, dimension (3) :: dYr
   real, dimension (3) :: X0
   real, dimension (3) :: Y0
   real, dimension (3) :: g1
   real, dimension (3) :: u1X
   real, dimension (3) :: uX0

!   real, dimension (nspec_max):: rcutoff_max
   real, dimension (numorb_max)   :: psi1
   real, dimension (3,numorb_max) :: dpsi1
   real, dimension (numorb_max)   :: psi2
   real, dimension (3,numorb_max) :: dpsi2
   real :: psiR
   real :: dpsiR
   real, dimension (5)            :: psiL
   real, dimension (3,5)          :: dpsiL

   real, target, dimension (:), allocatable  :: rhotmp
   real, dimension (3,3)          :: amat
   real, dimension (3,3)          :: inva

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

! set nr(:)
   nr(1) = rm1
   nr(2) = rm2
   nr(3) = rm3

! make a copy of the elem grid lattice vector
! we need to solve this linear eq.
!
!  | a1x  a2x  a3x |   |n1|   |x|
!  | a1y  a2y  a3y | x |n2| = |y|
!  ! a1z  a2z  a3z |   |n3|   |z|
!
! copy and invert original elvec to get form written above
   amat = transpose(elvec)
! inverse A: solving A*n=x -> n=A-1*x
   call inv3x3 (amat,inva)

   renorm = 0.0d0
! Loop over atoms
   do iatom = 1, natoms

    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)

! vector between the atom and initial point of the grid
    do i = 1,3
     u(i) = ratom(i,iatom) - g0(i)
    enddo ! i
! get n-vector
    call mult3x1 (inva,u)

! round coefficients to get the position of the nearest grid point g1 to the iatom X1
! i,j,k can be positive or negative exceeding rmX (it means not centered in the unit cell)
    i0 = nint( u(1) )
    j0 = nint( u(2) )
    k0 = nint( u(3) )

! find the vector u1 between the iatom X1 and the nearest point g1
    u1X(1) = u(1) - real(i0)
    u1X(2) = u(2) - real(j0)
    u1X(3) = u(3) - real(k0)

! check if the nearest grid point is located within the unit cell of the grid coords
! if not, let's map it within
!i0
    if (u(1) .lt. 0.0d0) then
     i0 = i0 + rm1*(int(abs(i0/rm1)) + 1)
    else
     i0 = i0 - rm1*int(i0/rm1)
    endif
!j0
    if (u(2) .lt. 0.0d0) then
     j0 = j0 + rm2*(int(abs(j0/rm2)) + 1)
    else
     j0 = j0 - rm2*int(j0/rm2)
    endif
!k0
    if (u(3) .lt. 0.0d0) then
     k0 = k0 + rm3*(int(abs(k0/rm3)) + 1)
    else
     k0 = k0 - rm3*int(k0/rm3)
    endif

! find the coordinates of the nearest point g1 witihin the grid coords
    g1(1) = i0*elvec(1,1) + j0*elvec(2,1) + k0*elvec(3,1)
    g1(2) = i0*elvec(1,2) + j0*elvec(2,2) + k0*elvec(3,2)
    g1(3) = i0*elvec(1,3) + j0*elvec(2,3) + k0*elvec(3,3)

! evaluate coordinates of the iatom in the grid coords
    X0(1) = g1(1) + u1x(1)*elvec(1,1) + u1x(2)*elvec(2,1) + u1x(3)*elvec(3,1)
    X0(2) = g1(2) + u1x(1)*elvec(1,2) + u1x(2)*elvec(2,2) + u1x(3)*elvec(3,2)
    X0(3) = g1(3) + u1x(1)*elvec(1,3) + u1x(2)*elvec(2,3) + u1x(3)*elvec(3,3)

! vector pointing from g1 to X0
    u1X(1) = g1(1) - X0(1)
    u1X(2) = g1(2) - X0(2)
    u1X(3) = g1(3) - X0(3)

! save iatom coord within the grid unit cell
    ratom2g(:,iatom) = X0(:)

! find index of the gX point within the extended mesh
    index0 = 1 + (i0+emx1) + em1*(j0+emx2) + em1*em2*(k0+emx3)

! Loop over the neighbors
    do ineigh = 1, neighn(iatom)

     jatom = neigh_j(ineigh,iatom)

      mbeta = neigh_b(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)

      do i = 1,3
       r21(i) = r2(i) - r1(i)
      enddo

! Loop over points in the atomic mesh gP
      do imesh = 1, nam
! restore index of the given mesh point gP within the extended mesh
       index = index0 + am2rc(imesh)

! evaluate the vector between the iatom and the mesh point gP
       do i = 1,3
        dXr(i) = ram2rc(i,imesh) + u1X(i)
       enddo
! evaluate the vector between the jatom and the mesh point gP
       do i = 1,3
        dYr(i) = dXr(i) - r21(i)
       enddo
! distance between the mesh point and jatom
       distY = sqrt(dYr(1)**2 + dYr(2)**2 + dYr(3)**2)

! check if jatom overlap with the gP mesh point
       if (distY .lt. Rc_max) then

! distance between the mesh point and iatom
        distX = sqrt(dXr(1)**2 + dXr(2)**2 + dXr(3)**2)

        psi1 = 0.0d0
! restore wavefunctions of iatom
        imu = 1
        do issh = 1,nssh(in1)
! get radial part of wf.
         call getpsi(in1,issh,distX,psiR,dpsiR)
! angular momentum
         l = lssh(issh,in1)
         call getYlm(l,dXr,psiL,dpsiL)
         do lmu = 1, (2*l+1)
! get spherical harmonics part of wf.
          psi1(imu) = psi1(imu) + psiL(lmu)*psiR
          imu = imu + 1
         enddo ! do lmu
        enddo ! do issh


        psi2 = 0.0d0
! restore wavefunctions of jatom
        imu = 1
        do issh = 1,nssh(in2)
! get radial part of wf.
         call getpsi(in2,issh,distY,psiR,dpsiR)
! angular momentum
         l = lssh(issh,in2)
         call getYlm(l,dYr,psiL,dpsiL)
         do lmu = 1, (2*l+1)
! get spherical harmonics part of wf.
          psi2(imu) = psi2(imu) + psiL(lmu)*psiR
          imu = imu + 1
         enddo ! do lmu
        enddo ! do issh

! map the point from the extended mesh into the normal mesh
	    ind = e2r(index) - 1

! assemble density
        dens = 0.0d0
        do inu = 1, num_orb(in1)
         do imu = 1, num_orb(in2)
          dens = dens + rho(inu,imu,ineigh,iatom)                &
      &        *psi1(inu)*psi2(imu)
         enddo ! do inu
        enddo ! do imu

! store variation of density at given point
        drhoG(ind) = dens + drhoG(ind)
       endif ! if (Rc_max)
      end do ! do imesh
    end do ! do ineigh
   end do ! do iatom

! test for hydrogen s-orbital the quality of the integration;
! the ratio it should go to one with higher Ecut
   dens = 4.0d0*3.141592653589793238462643*Rc_max**3/3.0d0

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
   do i =0,nrm-1
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

