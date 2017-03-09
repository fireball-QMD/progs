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


! psi22mesh.f90
! Program Description
! ===========================================================================
!      build up list of overlaping mesh points between two atoms
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
 subroutine psi22mesh ()

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
   integer npmax
   integer npatm
   integer ip


   real distX
   real distY


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


! set nr(:)
   nr(1) = rm1
   nr(2) = rm2
   nr(3) = rm3
! ======================================================================
! Find maximal number of overlaping points (exclude self-interaction)
! ======================================================================

   npmax = 0
! Loop over atoms
   do iatom = 1, natoms

! make a copy of the elem grid lattice vector
    do i = 1,3
     amat(i,:) = elvec(i,:)
    enddo
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)

! vector between the atom and initial point of the grid
    do i = 1,3
     u(i) = ratom(i,iatom) - g0(i)
    enddo ! i

! solve linear set of equations (3x3)
    call dgesl(amat,3,ipiv,u)
! round coefficients to get the position of the nearest point of the grid g1
    i0 = nint( u(1) )
    j0 = nint( u(2) )
    k0 = nint( u(3) )

! find the vector u between the iatom X0 and the nearest point g1
    u1X(1) = u(1) - real(i0)
    u1X(2) = u(2) - real(j0)
    u1X(3) = u(3) - real(k0)

! check if the nearest grid point is located within the unit cell of the grid
! if not, let map it there
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

! find the coordinates of gX
    g1(1) = i0*elvec(1,1) + j0*elvec(2,1) + k0*elvec(3,1)
    g1(2) = i0*elvec(1,2) + j0*elvec(2,2) + k0*elvec(3,2)
    g1(3) = i0*elvec(1,3) + j0*elvec(2,3) + k0*elvec(3,3)

! evaluate coordinates of iatom in the grid coord system
    X0(1) = g1(1) + u1X(1)*elvec(1,1) + u1X(2)*elvec(2,1) + u1X(3)*elvec(3,1)
    X0(2) = g1(2) + u1X(1)*elvec(1,2) + u1X(2)*elvec(2,2) + u1X(3)*elvec(3,2)
    X0(3) = g1(3) + u1X(1)*elvec(1,3) + u1X(2)*elvec(2,3) + u1X(3)*elvec(3,3)

! vector pointing from gX to X0
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
! skip iatom itself
     if (iatom .eq. jatom .and. mbeta .eq. 0) then
! nothing to do
     else
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)

      do i = 1,3
       r21(i) = r2(i) - r1(i)
      enddo

! reser counter of overlaping mesh points
      npatm = 0
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
       if (distY .lt. Rc_max) npatm = npatm + 1
      end do ! do imesh
! update max. number of overlaping points
      npmax = max(npatm, npmax)
     endif ! if (iatom)
    end do ! do ineigh
   end do ! do iatom

! the array reallocation regarding a change of neigh_max
! needs to be fixed for MD-loop properly
! probably movement to other block of the code.
   if (.not. allocated(psi22m)) then
    allocate (psi22m (numorb_max, numorb_max, npmax, neigh_max, natoms))
    allocate (ipsi22m (npmax, neigh_max, natoms))
    allocate (npsi22m (neigh_max, natoms))
   endif

! Loop over atoms
   do iatom = 1, natoms
! make a copy of the elem grid lattice vector
    do i = 1,3
     amat(i,:) = elvec(i,:)
    enddo
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)

! vector between the atom and initial point of the grid
    do i = 1,3
     u(i) = ratom(i,iatom) - g0(i)
    enddo ! i

! solve linear set of equations (3x3)
    call dgesl(amat,3,ipiv,u)
! round coefficients to get the position of the nearest point of the grid g1
    i0 = nint( u(1) )
    j0 = nint( u(2) )
    k0 = nint( u(3) )

! find the vector u between the iatom X0 and the nearest point g1
    u1X(1) = u(1) - real(i0)
    u1X(2) = u(2) - real(j0)
    u1X(3) = u(3) - real(k0)

! check if the nearest grid point is located within the unit cell of the grid
! if not, let map it there
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

! find the coordinates of gX
    g1(1) = i0*elvec(1,1) + j0*elvec(2,1) + k0*elvec(3,1)
    g1(2) = i0*elvec(1,2) + j0*elvec(2,2) + k0*elvec(3,2)
    g1(3) = i0*elvec(1,3) + j0*elvec(2,3) + k0*elvec(3,3)

! evaluate coordinates of iatom in the grid coord system
    X0(1) = g1(1) + u1X(1)*elvec(1,1) + u1X(2)*elvec(2,1) + u1X(3)*elvec(3,1)
    X0(2) = g1(2) + u1X(1)*elvec(1,2) + u1X(2)*elvec(2,2) + u1X(3)*elvec(3,2)
    X0(3) = g1(3) + u1X(1)*elvec(1,3) + u1X(2)*elvec(2,3) + u1X(3)*elvec(3,3)

! vector pointing from gX to X0
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
      ! skip iatom itself
     if (iatom .eq. jatom .and. mbeta .eq. 0) then
! nothing to do
     else

      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      do i = 1,3
       r21(i) = r2(i) - r1(i)
      enddo

! reset local index
      ip = 0
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

! increase local index
        ip = ip + 1

! reset variables
        psi1 = 0.0d0
        psi2 = 0.0d0

! evaluate the vector between the iatom and the mesh point gP
        distX = sqrt(dXr(1)**2 + dXr(2)**2 + dXr(3)**2)

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

! renormalize wavefunctions
        psi1(:) = psi1(:)*dloc
        psi2(:) = psi2(:)*dloc

! map the point from the extended mesh into the normal mesh
	    ind = e2r(index) - 1
	    ipsi22m (ip,ineigh,iatom) = ind


! store the wf values
        do inu = 1, num_orb(in1)
         do imu = 1, num_orb(in2)
          psi22m(inu,imu,ip,ineigh,iatom) = psi1(inu)*psi2(imu)
         enddo ! do inu
        enddo ! do imu


       endif ! if (Rc_max)
      end do ! do imesh
      npsi22m (ineigh,iatom) = ip
     endif ! if (iatom)
    end do ! do ineigh
   end do ! do iatom


! Format Statements
! ===========================================================================

   return
 end subroutine psi22mesh

