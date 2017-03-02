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

! assemble_KS_vna.f90
! Program Description
! ===========================================================================
!       Project the neutral atomic potential on the mesh.
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
 subroutine assemble_KS_vna ()

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
   integer iloc

   real distX
   real vna0
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
   vnaG = 0.0d0

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

! map the point from the extended mesh into the normal mesh
     ind = e2r(index) - 1

! get vna potential
     call getvna (in1, distX, vna0, dvna0)
     vnaG(ind) = vnaG(ind) + vna0
    end do ! do imesh
   end do ! do iatom

   vnaG = vnaG*eq2


! write out files
   if (iwrtxsf .eq. 1) then
! write out rho into xsf file
    pmat => vnaG
    filename = 'vna.xsf'
    mssg = 'VNA_3D'
    call writeout_xsf (filename, mssg, pmat)
   endif

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)

   return
 end subroutine assemble_KS_vna

