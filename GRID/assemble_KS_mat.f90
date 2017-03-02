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
 subroutine assemble_KS_mat (icluster)

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
   integer, intent (in) :: icluster


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
   integer ineighij,jneighij
   integer index
   integer index0
   integer ind
   integer indFD
   integer i, j, k
   integer i0, j0, k0
   integer lmu
   integer issh
   integer l
   integer imesh
   integer, dimension (3) :: ipiv
   integer, dimension (3) :: nr

   integer ii, jj, kk
   integer iD
   real maxD


   real distX
   real distY
   real qtot
   real renorm
   real dens
   real ddens
   real psi12
   real fx
   real fy
   real fz


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
   real, dimension (3) :: df

     

!   real, dimension (nspec_max):: rcutoff_max
   real, dimension (numorb_max)   :: psi1
   real, dimension (3,numorb_max) :: dpsi1
   real, dimension (numorb_max)   :: psi2
   real, dimension (3,numorb_max) :: dpsi2
   real :: psiR
   real :: dpsiR
   real, dimension (5)            :: psiL
   real, dimension (3,5)          :: dpsiL


!   real, dimension (3,natoms)     :: ratom2g
   real, dimension (3,3)          :: amat
   real, dimension (3,3)          :: inva

!   real, dimension (numorb_max,numorb_max)   :: integ1
!   real, dimension (numorb_max,numorb_max)   :: integ2
   real ex, mux, exc, muxc, dexc, d2exc, dmuxc, d2muxc
   real, dimension (numorb_max)   :: psiS
   real, dimension (5)            :: psiLS
   real, dimension (3,5)          :: dpsiLS

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

!   integ1 = 0.0d0
!   integ2 = 0.0d0


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
! modified by honza        
        ind = e2r(index) - 1
        indFD = n2e(ind)
        dens = drhoG(ind) + rhoG0(ind)

! first setup directional derivatives in df
! those are actually coordinate derivatives in terms of coordinates induced by the grid geometry        
        ddens = 0.0d0
        df = 0.0d0
        do ineighij = 1,3
           do jneighij = 1,2
              i = e2n(indFD + neighij(ineighij,jneighij))
              df(ineighij) = df(ineighij) + (drhoG(i) + rhoG0(i))  * d1f(ineighij,jneighij)
           end do
        end do

        fx = 0.0d0
        fy = 0.0d0
        fz = 0.0d0
! Now get the partial derivatives (i.e. directional derivatives in direction of canonical base vectors).
! Take note that the summation goes "over columns".        
        do i=1,3
           fx = fx + ervec(i,1)*df(i)
           fy = fy + ervec(i,2)*df(i)
           fz = fz + ervec(i,3)*df(i)
        end do
      
! Now sum all the partial derivatives together.
        ddens = fx + fy + fz

! modified by honza - end                   
        
        
! calc XC potential
        call ceperley_alder (dens, exc, muxc, dmuxc, dexc)
        vxcG(ind) = muxc

! assemble matrices
        do inu = 1, num_orb(in2)
         do imu = 1, num_orb(in1)
! vca  dV_H
           vca(imu,inu,ineigh,iatom) =  vca(imu,inu,ineigh,iatom)        &
!    &           + psi1(imu)*vcaG(ind)*psi2(inu)*dvol
    &           + psi1(imu)*vcaG(ind)*psi2(inu)*dvol
! vxc  Vxc
           vxc(imu,inu,ineigh,iatom) =  vxc(imu,inu,ineigh,iatom)        &
!    &           + psi1(imu)*vxcG(ind)*psi2(inu)*dvol
    &           + psi1(imu)*vxcG(ind)*psi2(inu)*dvol
          psi12 = psi1(imu)*psi2(imu)
! delta_rho
           vxc_ca(imu,inu,ineigh,iatom) =  vxc_ca(imu,inu,ineigh,iatom)  &
    &            + psi1(imu)*psi2(inu)*dvol
!    &            + psi1(imu)*vnaG(ind)*psi2(inu)*dvol
!          integ1(imu,inu) =  integ1(imu,inu) +                       &
!      &       psi1(imu)*(dmuxc*psiS(1)**2)*psi2(inu)*dvol
!      &       psi1(imu)*(dexc*psiS(1)**2)*psi2(inu)*dvol
!      &       psi1(imu)*muxc*psi2(inu)*dvol
!          integ2(imu,inu) =  integ2(imu,inu) +                       &
!      &       psi1(imu)*(dmuxc*psiS(5)**2)*psi2(inu)*dvol
!      &       psi1(imu)*(dexc*psiS(5)**2)*psi2(inu)*dvol
!      &       psi1(imu)*exc*psi2(inu)*dvol

         enddo ! do inu
        enddo ! do imu

       endif ! if (Rc_max)
      end do ! do imesh
    end do ! do ineigh
   end do ! do iatom
!vca=0.0d0
!vxc=0.0d0

   
! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
!101 format (<numorb_max>f14.6)
300 format (i8,f16.8)
301 format (i8,2f16.8)
   return
 end subroutine assemble_KS_mat

