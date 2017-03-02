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
!       Project bands on the mesh.
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
 subroutine ew2mesh (icluster)

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use kpoints
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in) :: icluster

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
   integer ikpoint
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
   integer job
   integer, dimension (3) :: ipiv
   integer, dimension (3) :: nr
   integer ii
   integer nnu
   integer mmu
   integer file
   integer iband
   integer nmax

   real distX
   real distY
   real dens
!   real densmax
   real dot
   real gutr

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

!   real, dimension (3,natoms)     :: ratom2g
   real, dimension (3,3)          :: lmat
   real, dimension (3,3)          :: invl

   real,    target, dimension (:), allocatable :: ewfaux
   integer nbandsin
   integer, dimension         (:), allocatable :: pbandsin
   integer, dimension         (:), allocatable :: pkpointsin

   complex phase
   complex step1
   complex step2
   complex ai

   character(40)   :: namewf
   character(4)    :: name
   real, dimension (:), pointer   :: pmat
   character (len=30) mssg

! Procedure
! ===========================================================================


   write (*,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
   write (*,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "

   write (*,*) " HEY ew2mesh here, iewform = ", iewform

   if (  iewform .eq. 5)                            call  ew2mesh_fourier (icluster)
   if (  iewform .eq. 6)                            call  ew2mesh_kscan (icluster)
   if (  iewform .eq. 7)                            call  ew2mesh_ARPES (icluster)
   if (( iewform .eq. 3) .OR. ((iewform .eq. 4)))   call  ew2mesh_gamma (icluster)

! reset variables
   drhoG = 0.0d0
   job = 0

! allocate aux arrays
   allocate ( ewfaux(0:nrm-1))

   ai = cmplx(0.0d0, 1.0d0)

! set nr(:)
   nr(1) = rm1
   nr(2) = rm2
   nr(3) = rm3

   if (iewform .eq. 2) then
    ! count bands within energy window
    nbandsin = 0
    do iband = 1, norbitals_new
      do ikpoint = 1, nkpoints
         if ((ewfewin_min .le. eigen_k(iband,ikpoint)) .and. (ewfewin_max .ge. eigen_k(iband,ikpoint))) then
         	nbandsin = nbandsin + 1
         endif
      enddo  ! do ikpoint
    enddo ! do iband
    if (nbandsin .eq. 0) then
     write (*,*) '   +++++++++++++++++++++++++++++++++++++++++++++++'
     write (*,*) '            NO bands within the interval '
     write (*,*) '              skip the band projection '
     write (*,*) '   +++++++++++++++++++++++++++++++++++++++++++++++'
     return
    endif
    allocate ( pbandsin   (nbandsin+1) )
    allocate ( pkpointsin (nbandsin+1) )
    pbandsin   (nbandsin+1) = -1
! find bands within energy window ( 2nd pass )
    ii = 1
    do iband = 1, norbitals_new
      do ikpoint = 1, nkpoints
         if ((ewfewin_min .le. eigen_k(iband,ikpoint)) .and. (ewfewin_max .ge. eigen_k(iband,ikpoint))) then
      		pbandsin  (ii) = iband
      		pkpointsin(ii) = ikpoint
            ii = ii + 1
         	write (*,'(A,i6,A,i6,A,f10.5)') 'Selected band ',iband,'  kpoint ',ikpoint,' E = ',eigen_k(iband,ikpoint) 
         endif
      enddo  ! do ikpoint
    enddo ! do iband

   else ! (iewform .eq. 2)
     nbandsin = npbands*nkpoints
     allocate ( pbandsin   (nbandsin+1) )
     allocate ( pkpointsin (nbandsin+1) )
     pbandsin   (nbandsin+1) = -1
     ii = 1
     do iband = 1,npbands
     	do ikpoint = 1,nkpoints
      	pbandsin  (ii) = pbands(iband)
      	pkpointsin(ii) = ikpoint
        ii = ii + 1
		end do  
     end do
   endif ! (iewform .eq. 2)

! ===============================
! =====     Projection     ======
! ===============================
 
write (*,*) " DEBUG : 1 "

   ewfaux = 0.0d0
   do ii = 1, nbandsin
!write (*,*) " DEBUG : 2 ii ", ii,iband,ikpoint
    iband   = pbandsin   (ii)
    ikpoint = pkpointsin (ii)

    write (*,'(A,i6,A,i6,A,f10.5)') 'Projecting band ',iband,'  kpoint ',ikpoint,' E = ',eigen_k(iband,ikpoint) 
!    densmax = 0.0d0

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

! Loop over the neighbors
      do ineigh = 1, neighn(iatom)

       jatom = neigh_j(ineigh,iatom)
       mbeta = neigh_b(ineigh,iatom)


       r2(:) = ratom(:,jatom) + xl(:,mbeta)
       in2 = imass(jatom)
       do i = 1,3
        r21(i) = r2(i) - r1(i)
       enddo

!====================================================
! ===== Loop over points in the atomic mesh gP ======
!====================================================

       do imesh = 1, nam
! restore index of the given mesh point gP within the extended mesh
        index = index0 + am2rc(imesh)
        do i = 1,3
         dXr(i) = ram2rc(i,imesh) + u1X(i) ! iatom to the mesh point gP
        enddo
        do i = 1,3
         dYr(i) = dXr(i) - r21(i) ! jatom to the mesh point gP
        enddo
        distY = sqrt(dYr(1)**2 + dYr(2)**2 + dYr(3)**2) ! dist. of mesh point to jatom

        if (distY .lt. Rc_max) then ! in cutoff
         distX = sqrt(dXr(1)**2 + dXr(2)**2 + dXr(3)**2) ! dist. of mesh point to jatom

! get wavefunctions of iatom
         imu = 1
         do issh = 1,nssh(in1)
          call getpsi(in1,issh,distX,psiR,dpsiR)  ! rarial basis func 
          l = lssh(issh,in1)
          call getYlm(l,dXr,psiL,dpsiL)  ! angular basis func
          do lmu = 1, (2*l+1)
           psi1(imu) = psiL(lmu)*psiR
           imu = imu + 1
          enddo ! do lmu
         enddo ! do issh

! get wavefunctions of jatom
         imu = 1
         do issh = 1,nssh(in2)
          call getpsi(in2,issh,distY,psiR,dpsiR)
          l = lssh(issh,in2)
          call getYlm(l,dYr,psiL,dpsiL)
          do lmu = 1, (2*l+1)
           psi2(imu) = psiL(lmu)*psiR
           imu = imu + 1
          enddo ! do lmu
         enddo ! do issh

! map the point from the extended mesh into the normal mesh
         ind = e2r(index) - 1

        dens = 0.0d0
        if (icluster .ne. 1) then
        ! Periodic System
          dot = special_k(1,ikpoint)*r21(1) + special_k(2,ikpoint)*r21(2) + special_k(3,ikpoint)*r21(3)
          phase = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)
          do imu = 1, num_orb(in1)
          	mmu = imu + degelec(iatom)
          	step1 = phase*(bbnkre(mmu,iband,ikpoint) - ai*bbnkim(mmu,iband,ikpoint))
         	do inu = 1, num_orb(in2)
          		nnu = inu + degelec(jatom)
         		step2 = step1*(bbnkre(nnu,iband,ikpoint) + ai*bbnkim(nnu,iband,ikpoint))
         		gutr = real(step2)
         		dens = dens + gutr*psi1(imu)*psi2(inu)
          end do ! imu
         end do ! inu
        else ! icluster
        ! Molecule
         do imu = 1, num_orb(in1)
         	mmu = imu + degelec(iatom)
         	step1 = bbnkre(mmu,iband,ikpoint)
         	do inu = 1, num_orb(in2)
         		nnu = inu + degelec(jatom)
         		step2 = step1*bbnkre(nnu,iband,ikpoint)
         		gutr = real(step2)
           		dens = dens + gutr*psi1(imu)*psi2(inu)
         		end do ! do inu
         	end do ! do imu
        end if ! if (icluster)
        ewfaux(ind) = ewfaux(ind) + dens
!       densmax = max(dens,densmax)

       endif ! if (Rc_max)

      end do ! do imesh
     end do ! do ineigh
    end do ! do iatom


! ===========================
! write out eigenfunctions into bandplotXXX.xsf file (format of xcrysden visual code)
!  for details see www.xcrysden.org

    if ((iewform .eq. 1) .and. ( iband .ne. pbandsin(ii+1) )) then
     file = 100 + iband
     write (name,'(i4.4)') iband
     namewf = 'bandplot_'//name//'.xsf'
     write (*,*) '  writting down band no.',iband,' into the file ',namewf
     pmat => ewfaux
     mssg = 'density_3D'
     call writeout_xsf (namewf, mssg, pmat)
     ewfaux = 0.0d0
    end if ! iewform =1

   enddo !do ii (iband)

   if (iewform .eq. 2) then
    namewf = 'bandplot_EW.xsf'
    pmat => ewfaux
    mssg = 'density_3D'
    call writeout_xsf (namewf, mssg, pmat)
   endif ! iewform = 2

   deallocate (ewfaux)
   deallocate (pbandsin)
   

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)

   return
 end subroutine ew2mesh

