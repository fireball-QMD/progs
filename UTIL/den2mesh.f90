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
 subroutine den2mesh (icluster)

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
   real qtot
   real renorm
   real dens

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
!   real, dimension (3,natoms)     :: ratom2g
   real, dimension (3,3)          :: lmat
   real, dimension (3,3)          :: invl

   real, dimension (:), pointer   :: pmat
   character (len=40) filename
   character (len=30) mssg


! Procedure
! ===========================================================================

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
   lmat = transpose(elvec)
! inverse A: solving A*n=x -> n=A-1*x
   call inv3x3 (lmat,invl)

   renorm = 0.0d0
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

! find index of the gX point within the extende mesh
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

! evaluate the vector between the iatom and the mesh point gP
       distX = sqrt(dXr(1)**2 + dXr(2)**2 + dXr(3)**2)
       psi1 = 0.0d0
! restore wavefunctions on iatom
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
! restore wavefunctions on jatom
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

! sum total charge of the system
  qtot = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      qtot = qtot + Qin(issh,iatom)
    end do
  end do

! the renormalization factor
  renorm = qtot/dens

! check total density after renormalization
  dens = 0.0d0
  do i =0,nrm-1
   drhoG(i) = drhoG(i)*renorm
   dens = dens + drhoG(i)*dvol
  enddo
  write (*,*) ' -- Total density after renormalization =',dens

! write out xsf files
! write out total density rho into xsf file
  pmat => drhoG
  filename = 'density.xsf'
  mssg = 'density_3D'
  call writeout_xsf (filename, mssg, pmat)

! get drho (rest atomic charge)
  drhoG = drhoG - rhoG0

! calc dipole with the unit cell
  index = 0
  dip_x =0.0d0
  dip_y =0.0d0
  dip_z =0.0d0
  do k = 0, rm3-1
   do j = 0, rm2-1
    do i = 0, rm1-1
     dqi = drhoG(index)
     x = i*elvec(1,1) + j*elvec(2,1) + k*elvec(3,1)
     y = i*elvec(1,2) + j*elvec(2,2) + k*elvec(3,2)
     z = i*elvec(1,3) + j*elvec(2,3) + k*elvec(3,3)
     dip_x = dip_x + x*dqi
     dip_y = dip_y + y*dqi
     dip_z = dip_z + z*dqi
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

! reminder
! drhoG ... variation of the SCF density (rho - rho0)
! rhoG0 ... the neutral atom density

! write out xsf file
! write out drho into xsf file
  pmat => drhoG
  filename = 'ddensity.xsf'
  mssg = 'ddensity_3D'
  call writeout_xsf (filename, mssg, pmat)

! write out vna into xsf file
  pmat => vnaG
  filename = 'vna.xsf'
  mssg = 'vna_3D'
  call writeout_xsf (filename, mssg, pmat)



! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)
301 format (2x,'Dipole_x =',e14.6,'  [D] ')
302 format (2x,'Dipole_y =',e14.6,'  [D] ')
303 format (2x,'Dipole_z =',e14.6,'  [D] ')
304 format (2x,'Dipole_tot =',e14.6,'  [D] ')

   return
 end subroutine den2mesh

! Updated 10/24/2001.
!
!cccccccccccccccccccccccc     Program 4.3     cccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
! Please Note:                                                         c
!                                                                      c
! (1) This computer program is part of the book, "An Introduction to   c
!     Computational Physics," written by Tao Pang and published and    c
!     copyrighted by Cambridge University Press in 1997.               c
!                                                                      c
! (2) No warranties, express or implied, are made for this program.    c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
     SUBROUTINE  DGESL (A,N,INDX,B)
!
! An example of solving linear equation set A(N,N)*X(N) = B(N)
! with the partial-pivoting Gaussian elimination scheme.  The
! numerical values are for the Wheatstone bridge example discussed
! in Section 4.1 in the book with all resistances being 100 ohms
! and the voltage 200 volts.
!
!!      PARAMETER (N=3)
      DIMENSION X(N),B(N),A(N,N),INDX(N)
!!      DATA B/200.0,0.0,0.0/,
!!     *     ((A(I,J), J=1,N),I=1,N) /100.0,100.0,100.0,-100.0,
!!     *                   300.0,-100.0,-100.0,-100.0, 300.0/
!!
      CALL LEGS (A,N,B,X,INDX)
!
!!      WRITE (6, 999) (X(I), I=1,N)
!!      STOP
        B(:) = X(:)

        RETURN
  999 FORMAT (F16.8)
      END
!
      SUBROUTINE LEGS(A,N,B,X,INDX)
!
! Subroutine to solve the equation A(N,N)*X(N) = B(N) with the
! partial-pivoting Gaussian elimination scheme.
!
      DIMENSION A(N,N),B(N),X(N),INDX(N)
!
      CALL ELGS(A,N,INDX)
!
      DO     100 I = 1, N-1
        DO    90 J = I+1, N
            B(INDX(J)) = B(INDX(J)) - A(INDX(J),I)*B(INDX(I))
   90   CONTINUE
  100 CONTINUE
!
      X(N) = B(INDX(N))/A(INDX(N),N)
      DO     200 I = N-1, 1, -1
        X(I) = B(INDX(I))
        DO   190 J = I+1, N
          X(I) = X(I)-A(INDX(I),J)*X(J)
  190   CONTINUE
          X(I) =  X(I)/A(INDX(I),I)
  200 CONTINUE
!
      RETURN
      END
!
      SUBROUTINE ELGS(A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed
! matrix plus the pivoting element ratios below the diagonal in
! the output.  INDX(N) records the pivoting order.

!
      DIMENSION A(N,N),INDX(N),C(N)
!
! Initialize the index
!
      DO     50    I = 1, N
        INDX(I) = I
   50 CONTINUE
!
! Find the rescaling factors, one from each row
!
        DO     100   I = 1, N
          C1= 0.0
          DO    90   J = 1, N
            C1 = AMAX1(C1,ABS(A(I,J)))
   90     CONTINUE
          C(I) = C1
  100   CONTINUE
!
! Search the pivoting (largest) element from each column
!
      DO     200   J = 1, N-1
        PI1 = 0.0
        DO   150   I = J, N
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
  150   CONTINUE
!
! Interchange the rows via INDX(N) to record pivoting order
!
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
        DO   170   I = J+1, N
          PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
          A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
          DO 160   K = J+1, N
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
  160     CONTINUE
  170   CONTINUE
  200 CONTINUE
!
      RETURN
      END

