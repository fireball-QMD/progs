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
   integer im, jm, km
   integer lmu
   integer issh
   integer l
   integer imesh
   integer info
   integer hit
   integer, dimension (3) :: ipiv
   integer, dimension (3) :: nr

   real distX
   real qtot
   real renorm
   real dens
   real vna0
   real psi2
   real, dimension (3) :: dvna0

   real, dimension (3) :: r1
   real, dimension (3) :: r2
   real, dimension (3) :: u
   real, dimension (3) :: dXr
   real, dimension (3) :: g1
   real, dimension (3) :: X0
   real, dimension (3) :: u1X
   real, dimension (3) :: uX0


   real, dimension (numorb_max)   :: psi1
   real, dimension (3,numorb_max) :: dpsi1
   real :: psiR
   real :: dpsiR
   real, dimension (5)            :: psiL
   real, dimension (3,5)          :: dpsiL
   real, dimension (:), pointer   :: pmat
   character (len=40) filename
   character (len=30) mssg

!   real, dimension (3,natoms)     :: ratom2g
   real, dimension (3,3)          :: lmat
   real, dimension (3,3)          :: invl

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
   vnaG = 0.0d0
   psi2 = 0.0d0

! set nr(:)
   nr(1) = rm1
   nr(2) = rm2
   nr(3) = rm3

! make a copy of the elem grid lattice vector
! elvec
!  | a1x a1y a1z |
!  | a2x a2y a2z |
!  | a3x a3y a3z |
! position of atom is described as:
!  x = n1*a1x + n2*a2x + n3*a3x
! we need to solve this linear eq.
!
!  | a1x  a2x  a3x |   |n1|   |x|
!  | a1y  a2y  a3y | x |n2| = |y|
!  ! a1z  a2z  a3z |   |n3|   |z|
!
! copy and invert original elvec to get form written above
   lmat = transpose(elvec)
! inverse A: solving A*n=x -> n=A-1*x
   call inv3x3 (lmat,invl)

    write (*,*) 'elvec'
    write (*,400) elvec(1,:)
    write (*,400) elvec(2,:)
    write (*,400) elvec(3,:)
! integration checking
   renorm = 0.0d0
   hit = 0

! Loop over atoms
   do iatom = 1, natoms

    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)

! vector between the iatom (not centered in the unit cell yet) and
! the initial grid point
    do i = 1,3
     u(i) = ratom(i,iatom) - g0(i)
    enddo ! i

! get n-vector
    call mult3x1 (invl,u)
! check it
    u1X = u
!    write (*,*) ' check in'
!    write (*,*) 'lmat'
!    write (*,400) lmat(1,:)
!    write (*,400) lmat(2,:)
!    write (*,400) lmat(3,:)
!   write (*,*) 'vec'
!    write (*,*) ratom(:,iatom) - g0(:)
!    write (*,*) 'ni'
!    write (*,400) u1X(:)
    call mult3x1 (lmat,u1X)
    do i = 1,3
      if (abs(u1X(i)-(ratom(i,iatom) - g0(i))) .gt. 0.000001d0 ) then
         write (*,*) ' check out'
         write (*,*) 'vec'
         write (*,400) u1X(:)
         write (*,*) 'diff = ',abs(u1X(i)-(ratom(i,iatom) - g0(i)))
         write (*,*) 'stop'
         stop
       endif
    enddo

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

!    do i = 1,3
!      if (abs(X0(i)-(ratom(i,iatom) - g0(i))) .gt. 0.000001d0 ) then
!         write (*,*) ' inconsistent'
!         write (*,400) X0(:)
!         write (*,400) (ratom(:,iatom) - g0(:))
!         write (*,*) 'diff = ',abs(X0(i)-(ratom(i,iatom) - g0(i)))
!         write (*,*) 'stop'
!         stop
!       endif
!    enddo

! vector pointing from g1 to X0
    u1X(1) = g1(1) - X0(1)
    u1X(2) = g1(2) - X0(2)
    u1X(3) = g1(3) - X0(3)

! save iatom coord within the grid unit cell
    ratom2g(:,iatom) = X0(:)

! find index of the gX point within the extende mesh
    index0 = 1 + (i0+emx1) + em1*(j0+emx2) + em1*em2*(k0+emx3)

! skip certain atoms to avoid the overcounting; it means all jatoms having
! the identification number less then iatom because those pairs
! have been counted previously


! Loop over points in the atomic mesh gP
    do imesh = 1, nam

! restore index of the given mesh point gP within the extended mesh
     index = index0 + am2rc(imesh)

! evaluate the vector between the iatom X0 and the mesh point gP
     do i = 1,3
      dXr(i) = ram2rc(i,imesh) + u1X(i)
     enddo

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

! map the point from the extended mesh into the normal mesh
     ind = e2r(index) - 1
! assemble density
     dens = 0.0d0
     do imu = 1, num_orb(in1)
      dens = dens + rhoA(imu,iatom)*psi1(imu)*psi1(imu)
      renorm = renorm + dvol
      hit = hit + 1
      psi2 = psi2 + psi1(imu)*psi1(imu)
     enddo ! do imu
! get vna potential
     call getvna (in1, distX, vna0, dvna0)
     vnaG(ind) = vnaG(ind) + vna0
! store variation of density at given point
     rhoG0(ind) = dens + rhoG0(ind)

    end do ! do imesh
   end do ! do iatom

   vnaG = vnaG*eq2

! test for hydrogen s-orbital the quality of the integration;
! the ratio it should go to one with higher Ecut
   write (*,*) 'Rc_max =',Rc_max
   write (*,*) ' Psi^2 =',psi2*dvol
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
   do i =0,nrm-1
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
400 format (3f12.6)
410 format (3i5)

   return
 end subroutine assemble_KS_den0



 subroutine inv3x3 (amatrix, ainverse)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent (in), dimension (3, 3) :: amatrix

! Output:
        real, intent (out), dimension (3, 3) :: ainverse

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        real determinant

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
        determinant =  &
     &    amatrix(1,1)*(amatrix(2,2)*amatrix(3,3) - amatrix(2,3)*amatrix(3,2))&
     &  + amatrix(1,2)*(amatrix(3,1)*amatrix(2,3) - amatrix(2,1)*amatrix(3,3))&
     &  + amatrix(1,3)*(amatrix(2,1)*amatrix(3,2) - amatrix(3,1)*amatrix(2,2))

! Now calculate inverse of inertia tensor
! ainv(i,j) = (-1**(i+j)) * cofactor(j,i) / det(a)

        if (abs(determinant) .gt. 1.0d-5) then
         ainverse(1,1) =    &
     &    (amatrix(2,2)*amatrix(3,3) - amatrix(3,2)*amatrix(2,3))/determinant
         ainverse(2,1) =    &
     &    - (amatrix(2,1)*amatrix(3,3) - amatrix(3,1)*amatrix(2,3))/determinant
         ainverse(3,1) =    &
     &    (amatrix(2,1)*amatrix(3,2) - amatrix(3,1)*amatrix(2,2))/determinant
         ainverse(1,2) =    &
     &    - (amatrix(1,2)*amatrix(3,3) - amatrix(3,2)*amatrix(1,3))/determinant
         ainverse(2,2) =    &
     &    (amatrix(1,1)*amatrix(3,3) - amatrix(3,1)*amatrix(1,3))/determinant
         ainverse(3,2) =    &
     &    - (amatrix(1,1)*amatrix(3,2) - amatrix(3,1)*amatrix(1,2))/determinant
         ainverse(1,3) =    &
     &    (amatrix(1,2)*amatrix(2,3) - amatrix(1,3)*amatrix(2,2))/determinant
         ainverse(2,3) =    &
     &    - (amatrix(1,1)*amatrix(2,3) - amatrix(1,3)*amatrix(2,1))/determinant
         ainverse(3,3) =    &
     &    (amatrix(2,2)*amatrix(1,1) - amatrix(1,2)*amatrix(2,1))/determinant

       else
         ainverse = 0.0d0
         write (*,*) ' ********* WARNING ********* '
         write (*,*) ' The determinant of the amatrix in invert3x3 is '
         write (*,*) ' equal to zero. Be careful and make sure that '
         write (*,*) ' you really want to continue. '
        end if

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
 end subroutine

 subroutine mult3x1 (amat, bvec)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent (in), dimension (3, 3) :: amat

! Output:
        real, intent (inout), dimension (3) :: bvec

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        real, dimension (3) :: cvec

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================

! multiply

! 1. column
       cvec(1) = amat (1,1)*bvec(1) + amat(1,2)*bvec(2) + amat(1,3)*bvec(3)
       cvec(2) = amat (2,1)*bvec(1) + amat(2,2)*bvec(2) + amat(2,3)*bvec(3)
       cvec(3) = amat (3,1)*bvec(1) + amat(3,2)*bvec(2) + amat(3,3)*bvec(3)



! copy the result into bmat matrix
       bvec = cvec

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
 end subroutine

 subroutine mult3x3 (amat, bmat)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent (in), dimension (3, 3) :: amat

! Output:
        real, intent (inout), dimension (3,3) :: bmat

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        real, dimension (3,3) :: cmat

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================

! multiply

! 1. column
       cmat(1,1) = amat (1,1)*bmat(1,1) + amat(1,2)*bmat(2,1) + amat(1,3)*bmat(3,1)
       cmat(1,2) = amat (1,1)*bmat(1,2) + amat(1,2)*bmat(2,2) + amat(1,3)*bmat(3,2)
       cmat(1,3) = amat (1,1)*bmat(1,3) + amat(1,2)*bmat(2,3) + amat(1,3)*bmat(3,3)
! 2. column
       cmat(2,1) = amat (2,1)*bmat(1,1) + amat(2,2)*bmat(2,1) + amat(2,3)*bmat(3,1)
       cmat(2,2) = amat (2,1)*bmat(1,2) + amat(2,2)*bmat(2,2) + amat(2,3)*bmat(3,2)
       cmat(2,3) = amat (2,1)*bmat(1,3) + amat(2,2)*bmat(2,3) + amat(2,3)*bmat(3,3)
! 3. column
       cmat(3,1) = amat (3,1)*bmat(1,1) + amat(3,2)*bmat(2,1) + amat(3,3)*bmat(3,1)
       cmat(3,2) = amat (3,1)*bmat(1,2) + amat(3,2)*bmat(2,2) + amat(3,3)*bmat(3,2)
       cmat(3,3) = amat (3,1)*bmat(1,3) + amat(3,2)*bmat(2,3) + amat(3,3)*bmat(3,3)


! copy the result into bmat matrix
       bmat = cmat

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
 end subroutine

