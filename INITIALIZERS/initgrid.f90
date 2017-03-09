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


! initgrid.f90
! Program Description
! ===========================================================================
!       This routine initializes the real-space grid defined in terms of
! Ecut = (Pi /dx)^2 / 2.0 [a.u.]
! ===========================================================================
! Code written by:
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine initgrid (icluster)

   use grid
   use constants_fireball
   use interactions
   use configuration
   use dimensions
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in)                     :: icluster

! Output

! Local Parameters and Data Declaration
! ===========================================================================

   real, parameter :: droff = 0.2d0
!   real, parameter :: droff = 1.5d0

! Local Variable Declaration and Description
! ===========================================================================

   integer n1
   integer n2
   integer n3
   integer ne1
   integer ne2
   integer i
   integer j
   integer k
   integer ix
   integer index
   integer index0
   integer ispec
   integer ipoint
   integer iatom
   integer ii
   integer jj
   integer kk

   real r1
   real r2
   real r3
   real vol
   real dist
   real xmin
   real xmax
   real denom
   real drmax
   real coxm

   real, dimension (3)       :: cmass
   real, dimension (3)       :: cvec
   real, dimension (3)       :: xvec
   real, dimension (3,3)     :: avec
   real, dimension (3,3)     :: rvec
   integer, dimension (3)    :: ipiv
   integer, dimension (3)    :: np
   real, dimension (3)     :: er1vec
   real, dimension (3)     :: er2vec
   real, dimension (3)     :: er3vec

   
   real, dimension (3)       :: acrossb
   real, dimension (3)       :: bcrossc
   real, dimension (3)       :: ccrossa


! Allocate arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Compute grid size from cutoff energy.

! find Rc_max
   Rc_max = 0.0d0
   do ispec = 1, nspecies
    do i = 1, nssh(ispec)
     if (rcutoff(ispec,i) .gt. Rc_max) Rc_max = rcutoff(ispec,i)
    end do
   end do

   
! calculate volume of elementar unit cell
   cvec(1)= a1vec(2)*a2vec(3) - a1vec(3)*a2vec(2)
   cvec(2)= a1vec(3)*a2vec(1) - a1vec(1)*a2vec(3)
   cvec(3)= a1vec(1)*a2vec(2) - a1vec(2)*a2vec(1)
   vol = abs (a3vec(1)*cvec(1) + a3vec(2)*cvec(2) + a3vec(3)*cvec(3))


! non-periodic boundary conditions
   if (icluster .eq. 1) then

    avec = 0.0d0
    do ix = 1,3
! find max and min position of atom in ix-axis direction
     xmax = -8000.0d0
     xmin = 8000.0d0
     do iatom = 1,natoms
      if (xmax .lt. ratom(ix,iatom)) xmax = ratom(ix,iatom)
      if (xmin .gt. ratom(ix,iatom)) xmin = ratom(ix,iatom)
     enddo ! do iatom

! define lattice vector
     avec(ix,ix) = xmax - xmin + 2*Rc_max + 2*droff
! define intial point
     if (ifixg0 .eq. 0) g0(ix) = xmin - Rc_max - droff

    enddo ! do ix

! setup lattice vector
    a1vec(:) = 0.0d0
    a1vec(1) = avec(1,1)
    a2vec(:) = 0.0d0
    a2vec(2) = avec(2,2)
    a3vec(:) = 0.0d0
    a3vec(3) = avec(3,3)

! find center of mass of atoms in the unit cell
    do i = 1,natoms
     do j = 1,3
      cmass(j) = cmass(j) + ratom(j,i)
     enddo ! do j
    enddo ! do i
    cmass(:) = cmass(:)/float(natoms)

   else

    avec(1,:) = a1vec(:)
    avec(2,:) = a2vec(:)
    avec(3,:) = a3vec(:)
! find center of mass of atoms in the unit cell
    do i = 1,natoms
     do j = 1,3
      cmass(j) = cmass(j) + ratom(j,i)
     enddo ! do j
    enddo ! do i
    cmass(:) = cmass(:)/float(natoms)

    if (ifixg0 .eq. 0) then
! find  initial point of the mesh. defined through the center of mass
     do i = 1,3
      g0(i) = cmass(i) - 0.5d0*(a1vec(i) + a2vec(i) + a3vec(i))
     enddo
    endif ! if (ifixg0 .eq. 1)

   endif ! if (icluster .eq. 1)

! Gauss diagonalization with pivoting
   write (*,*) ' ---------------------------------'
   write (*,*) '          Lattice vector          '
   write (*,*) ' ---------------------------------'
   do i = 1,3
    write (*,*) (avec(i,j),j=1,3)
   enddo
   write (*,*)  'atomic units'
   do i = 1,3
    write (*,*) (avec(i,j)/abohr,j=1,3)
   enddo
! Find reciprocal cell vectors (multiplied by 2*pi)
   call crossx(a2vec,a3vec,bcrossc)
   call crossx(a3vec,a1vec,ccrossa)
   call crossx(a1vec,a2vec,acrossb)
   denom=a1vec(1)*bcrossc(1)+a1vec(2)*bcrossc(2)+a1vec(3)*bcrossc(3)
   do i=1,3
    rlvec(1,i)=2.0d0*pi*bcrossc(i)/denom
    rlvec(2,i)=2.0d0*pi*ccrossa(i)/denom
    rlvec(3,i)=2.0d0*pi*acrossb(i)/denom
   enddo

   
   write (*,*)  'Reciporocal Lattice vector [Angstrom]:'
   do i = 1, 3
   write (*,*) (rlvec(i,j),j=1,3)
   enddo
   write (*,*) 'Volume [Ang^3]= ',denom

   write (*,*)  'Reciporocal Lattice vector [abohr]:'
   do i = 1, 3
   write (*,*) (rlvec(i,j)*abohr,j=1,3)
   enddo



! calc number of mesh points along axis
   do i = 1,3
! calc the size of the unit cell along the i-axis  
    cvec(i) = sqrt (rlvec(i,1)**2 + rlvec(i,2)**2 + rlvec(i,3)**2)      
! estimate the number of points along the i-axis
    np(i) = int (2 * sqrt(Ecut) / (cvec(i)*abohr) + 1)     

! iterate until the number of points is multiplier of 2
    do
! number of points must be multiply of 2!!
      if ( mod( np(i), 2 ) .eq. 0 ) exit
      np(i) = np(i) + 1
    enddo

! in the reciprocal space one has to multiply by np
    ervec(i,:) = (rlvec(i,:) * np(i)) / (2.0d0*pi)

! the elementary distance at given axis
! Should we multiply by np  here as well?
! (In order to compute drmax a few lines later)
    cvec(i) = cvec(i) / np(i)
 enddo ! do i

 
 er1vec(:) = ervec(1,:)
 er2vec(:) = ervec(2,:)
 er3vec(:) = ervec(3,:)

! Find direct elementary vectors
   call crossx(er2vec,er3vec,bcrossc)
   call crossx(er3vec,er1vec,ccrossa)
   call crossx(er1vec,er2vec,acrossb)
   denom=er1vec(1)*bcrossc(1)+er1vec(2)*bcrossc(2)+er1vec(3)*bcrossc(3)
   do i=1,3
    elvec(1,i)=bcrossc(i)/denom
    elvec(2,i)=ccrossa(i)/denom
    elvec(3,i)=acrossb(i)/denom
 enddo

    write (*,*)  'Reciporocal Elementary Lattice vector [Angstrom]:'
   do i = 1, 3
   write (*,*) (ervec(i,j),j=1,3)
   enddo
   write (*,*) 'Volume [Ang^3]= ',denom

   write (*,*)  'Reciporocal Elementary Lattice vector [abohr]:'
   do i = 1, 3
   write (*,*) (ervec(i,j),j=1,3)
   enddo
     

! save max dr
   drmax = 0.0d0
   do i = 1,3
    if (drmax .gt. cvec(i)) drmax = cvec(i)
   enddo

! store the division of the regular mesh along the axis
   rm1 = np(1)
   rm2 = np(2)
   rm3 = np(3)

! total number of points on the regular mesh
   nrm = rm1 * rm2 * rm3

! elementary volume
   dvol = abs(elvec(1,1)*(elvec(2,2)*elvec(3,3)-elvec(2,3)*elvec(3,2))    &
 &           +elvec(1,2)*(elvec(2,3)*elvec(3,1)-elvec(2,1)*elvec(3,3))    &
 &           +elvec(1,3)*(elvec(2,1)*elvec(3,2)-elvec(2,2)*elvec(3,1)))

   
! define the extended mesh at each direction
   do i = 1,3
! get the minimal distance
      
    dist = sqrt ( elvec(i,1)**2 + elvec(i,2)**2 + elvec(i,3)**2 )
! calc the trashold
    ipiv(i) = Rc_max / dist + 1
! total points in extended mesh
    np(i) = np(i) + 2*ipiv(i)
   enddo ! do i

! save offset of the extended mesh
   emx1 = ipiv(1)
   emx2 = ipiv(2)
   emx3 = ipiv(3)

! save division of the extended mesh
   em1 = np(1)
   em2 = np(2)
   em3 = np(3)
   nem = em1*em2*em3

   write (*,*) '  ---------  Begin mesh Info --------- '
   write (*,*) ''
   write (*,'(a,f8.6)') ' Rc_max = ',Rc_max
   write (*,300) Ecut
   write (*,380) (g0(i),i=1,3)
   write (*,'(a,4i9)') ' Regular mesh: ',rm1, rm2, rm3, nrm
   write (*,'(a,f16.8)') 'dVol :',dvol
   write (*,'(a,4i9)') ' Extended mesh: ',em1, em2, em3, nem
   write (*,*) ' Elementary grid lattice vector :'
   write (*,2000)  (elvec(1,i),i=1,3)
   write (*,2000)  (elvec(2,i),i=1,3)
   write (*,2000)  (elvec(3,i),i=1,3)
   write (*,*) '  ---------   End mesh Info --------- '
! ==============================================================
! =======   Map the Extended mesh to the Regular grid   ========
! ==============================================================

   
! allocate arrays
  allocate (e2r(nem))

! Loops over each axis
  do k = 0, em3-1
   do j = 0, em2-1
    do i = 0, em1-1

! rest the offset
     ii = i - emx1
     jj = j - emx2
     kk = k - emx3

! find the point in the regular mesh
     ii = mod( ii + 100*rm1, rm1)
     jj = mod( jj + 100*rm2, rm2)
     kk = mod( kk + 100*rm3, rm3)
    
! calc index of the extended mesh point
     index = 1 + i + em1*j + em1*em2*k
! calc index of the regular mesh point
     index0 = 1 + ii + rm1*jj + rm1*rm2*kk
     
! save the value of the point in the regular mesh context
     e2r(index) = index0

    enddo ! do i
   enddo ! do j
  enddo ! do k

! ==============================================================
! =======             Setup Atomic Mesh                 ========
! ==============================================================
! generate atomic mesh (sphere) within Rc_max
   write (*,*) 'atomic mesh'
   nam = 0
   do k = -emx3,emx3
    do j = -emx2,emx2
     do i = -emx1,emx1
        
! distance to the mesh cell
      do ix = 1,3
       cvec(ix) = elvec(1,ix)*i + elvec(2,ix)*j + elvec(3,ix)*k
      enddo
      dist = sqrt (cvec(1)**2 + cvec(2)**2 + cvec(3)**2)
! distance from point to mesh cell
      if ((Rc_max + drmax) .gt. dist) then
       nam = nam + 1
      endif
     enddo ! do i
    enddo ! do j
   enddo ! do k

! allocate arrays related to atomic mesh
   allocate (am2rc(nam))
   allocate (ram2rc(3,nam))
   allocate (ratom2g(3,natoms))

   index = 0
   do k = -emx3,emx3
    do j = -emx2,emx2
     do i = -emx1,emx1

! find the total coordinate of the point
       do ix = 1,3
        cvec(ix) = elvec(1,ix)*i + elvec(2,ix)*j + elvec(3,ix)*k
       enddo
! distance to the mesh cell
       dist = sqrt (cvec(1)**2 + cvec(2)**2 + cvec(3)**2)

! ?? is the point within the Rc_max ??
       if ((Rc_max + drmax) .gt. dist) then
        index = index + 1         
        am2rc(index) = i + em1*j + em1*em2*k       
        ram2rc(:,index) = cvec(:)
       endif
     enddo ! do i
    enddo ! do j
   enddo ! do k

   write (*,'(a,4i9)') ' Atomic mesh: ',2*emx1 + 1,2*emx2 + 1, 2*emx3&
        & + 1, nam

! ==============================================================
! =======              SET UP FDM GRID                  ========
! ==============================================================
! set extended mesh (used for finite difference method)
   noff = 1
   mfd1 = rm1 + 2*noff
   mfd2 = rm2 + 2*noff
   mfd3 = rm3 + 2*noff
   nmfd = mfd1 * mfd2 * mfd3

! write out information about grids
   write (*,*) '  ----    the FDM grid informations   ----- '
   write (*,300) Ecut
   write (*,*) '  Normal FDM mesh: ',nrm
   write (*,400) rm1,rm2,rm3
   write (*,*) '  Extended FDM mesh: ', nmfd
   write (*,400) mfd1,mfd2,mfd3

! modified by honza   
! set up auxilliary array keeping track of neighbors
! list of neighbors - 7,8 are currently redundant
   nneighij = 8
! relative coordinates
   neighij(1,1) = 1
   neighij(1,2) = -1
   neighij(2,1) = mfd1
   neighij(2,2) = -1*mfd1
   neighij(3,1) = mfd1*mfd2
   neighij(3,2) = -1*mfd1*mfd2
   neighij(4,1) = 0
   neighij(4,2) = 0

  
   
! recalculate the real space step of the mesh & convert it into a.u.
! remeber: only rectangular mesh can be used !!



! experimenting with (possibly) non-rectangular grid

! Vectors elvec(i,:) form the covariant base of the geometry.
! Vectors ervec(i,:) form the contravariant base of the geometry.
!   
! dvol is actually kind of volume form.
! Explicitly: dvol = sqrt(det(metric)) (where metric would be the gramm-matrix (<elvec(i,:),elvec(j,:)>;i,j=1,3)
!

! 1st derivative coefficients   
! No need to introduce any normalization factor - it's already present in the base vectors.
   d1f(1,1) = 1.0d0 / 2.0d0
   d1f(1,2) = -1.0d0 / 2.0d0
   d1f(2,1) = 1.0d0 / 2.0d0
   d1f(2,2) = -1.0d0 / 2.0d0
   d1f(3,1) = 1.0d0 / 2.0d0
   d1f(3,2) = -1.0d0 / 2.0d0
   d1f(4,1) = 0.0d0
   d1f(4,2) = 0.0d0


! end-modified by honza
   
! allocate arrays
   allocate (e2n (0:nmfd-1))
   allocate (n2e (0:nrm-1))

! allocate arrays
   allocate (vnaG (0:nrm-1))
   allocate (drhoG (0:nrm-1))
   allocate (rhoG0 (0:nrm-1))
   allocate (vcaG (0:nrm-1))
   allocate (vxcG (0:nrm-1))

   vnaG = 0.0d0
   drhoG = 0.0d0
   rhoG0 = 0.0d0
   vcaG = 0.0d0
   vxcG = 0.0d0

! --------------------------------------------------------------
!         set up mapping NORMAL(REGULAR) mesh INTO EXTENDED-FD mesh
! --------------------------------------------------------------
   
   index = 0
   index0 = mfd1*(mfd2+1)
   do k = 0, rm3-1
      do j = 0, rm2-1
         do i = 0, rm1-1
            index0 = index0 + 1
            n2e(index) = index0
!            write (*,*) ' norm : ',index,'ext : ',n2e(index)
            index = index + 1
         enddo ! do i
         index0 = index0 + 2*noff
      enddo ! do j
      index0 = index0 + 2*mfd1
   enddo ! do k

!*************************************************************
! --------------------------------------------------------------
!         Mapping EXTENDED-FD mesh TO NORMAL(REGULAR) mesh
! --------------------------------------------------------------   

  do k = 0, mfd3-1
   do j = 0, mfd2-1
    do i = 0, mfd1-1

! rest the offset
     ii = i - noff
     jj = j - noff
     kk = k - noff

! find the point in the basic mesh
     ii = mod( ii + 1000*rm1, rm1)
     jj = mod( jj + 1000*rm2, rm2)
     kk = mod( kk + 1000*rm3, rm3)

! calc index of the fdm mesh point
     index =  i + mfd1*j + mfd1*mfd2*k
! calc index of the regular mesh point
     index0 = ii + rm1*jj + rm1*rm2*kk

! save the value of the point in the basic mesh context
     e2n(index) = index0

    enddo ! do i
   enddo ! do j
  enddo ! do k




! Format Statements
! ===========================================================================
100     format (2x, i5, 2x, a40, 2x, i2)
200     format (2x, 70('='))
225     format (2x,  3f15.7)
250     format (2x, ' Vol  = ', f16.6, ' [Ang^3] ')
300     format (2x, ' Ecut = ', f16.6, ' [Ry] ')
350     format (2x, ' dr = ', 3f16.8, ' [Ang]')
370     format (2x, ' dV = ', f16.8, ' [Ang^3]')
380     format (2x, ' Virtual initial grid point :', 3f16.8)
400     format (2x, ' Number of points = '3i9)
500     format (2x, ' Virtual center of mass : '3f16.8)
2000    format (2x, '  '3f16.8)

   return
 end subroutine initgrid

      subroutine crossx(a,b,c)
       implicit double precision(a-h,o-z)
       dimension a(3),b(3),c(3)
! computes c = a X b.
       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
      end


!      SUBROUTINE RECLAT (A,B,IOPT)
!
!!  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
!!  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1
!
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      DOUBLE PRECISION A(3,3),B(3,3)
!      PI=ACOS(-1.D0)
!      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
!      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
!      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
!      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
!      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
!      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
!      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
!      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
!      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
!      C=1.D0
!      IF (IOPT.EQ.1) C=2.D0*PI
!      DO 20 I=1,3
!         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
!         B(1,I)=B(1,I)*CI
!         B(2,I)=B(2,I)*CI
!         B(3,I)=B(3,I)*CI
!  20  CONTINUE
!      END


      SUBROUTINE NFFT( N )

! CHANGES N INTO THE NEXT INTEGER ALLOWED BY THE FFT ROUTINE

      PARAMETER (NP = 3, NMAX = 1000000)
      INTEGER IPRIME(NP)
      DATA IPRIME / 2, 3, 5 /

      NMIN = N
      DO N = NMIN, NMAX
        NREM = N
        DO IP = 1,NP
   10     CONTINUE
          IF ( MODULO( NREM, IPRIME(IP) ) .EQ. 0 ) THEN
            NREM = NREM / IPRIME(IP)
            GOTO 10
          ENDIF
        ENDDO
        IF (NREM .EQ. 1) RETURN
      ENDDO
      WRITE(6,*) 'NFFT: NO SUITABLE INTEGER FOUND FOR N =', NMIN
      STOP
      END

