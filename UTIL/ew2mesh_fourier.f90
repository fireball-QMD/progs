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
 subroutine ew2mesh_fourier (icluster)

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use kpoints

   implicit none
   INCLUDE 'fftw3.f'  


! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in) :: icluster

!Output


! Local Parameters and Data Declaration
! ===========================================================================

   INTEGER, PARAMETER :: DPC = 8

   interface
    subroutine writeout_xsf (xsfname, message, aa)
     real, dimension (:), pointer, intent (in) :: aa
     character (len=40) xsfname
     character (len=30) message
    end
  function inter(pvar,x,y,z)
   !real, dimension (:,:,:), pointer, intent (in) :: pvar
   COMPLEX(8), ALLOCATABLE, DIMENSION(:,:,:),intent (in)  :: pvar
   real, intent (in) :: x
   real, intent (in) :: y
   real, intent (in) :: z
   real inter
   end
  end interface

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer in1
   integer ikpoint
   integer imu
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
   integer mmu
   integer file
   integer iband
   integer nmax

   real distX
   real dens
!   real densmax

   real, dimension (3) :: r1
   real, dimension (3) :: u
   real, dimension (3) :: dXr
   real, dimension (3) :: X0
   real, dimension (3) :: Y0
   real, dimension (3) :: g1
   real, dimension (3) :: u1X
   real, dimension (3) :: uX0

!   real, dimension (nspec_max):: rcutoff_max
   real, dimension (numorb_max)   :: psi1
   real, dimension (3,numorb_max) :: dpsi1
   real :: psiR
   real :: dpsiR
   real, dimension (5)            :: psiL
   real, dimension (3,5)          :: dpsiL

!   real, dimension (3,natoms)     :: ratom2g
   real, dimension (3,3)          :: lmat
   real, dimension (3,3)          :: invl

   real, target, dimension (:), allocatable :: ewfaux
   real, target, dimension (:), allocatable :: ewfRe
   real, target, dimension (:), allocatable :: ewfIm
  ! real, target, dimension (:), allocatable :: ewfacum

!  export xsf 
   character(40)   :: namewf
   character(4)    :: name
   character(15)    :: name2
   real, dimension (:), pointer   :: pmat
   character (len=30) mssg

!  print image
   real maxv
   integer cR,cB 
   real, dimension (:,:), allocatable :: positive
   real, dimension (:,:), allocatable :: negative

! ========= Fourier ==========

 integer hrm

 !variable used by FFTW
   INTEGER(DPC) :: plan
   
   COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: in3D
  ! COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: four3D
   COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: four3DX
   ! COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: acum3D
   ! real, target, allocatable, dimension (:) :: resf


! ========= get max variables

  real ValMax
  integer kmax, jmax, imax
  real rkmax,rjmax,rimax
  real recstep

! ========= control varibles

 integer, ALLOCATABLE, DIMENSION (:) :: myBands

 real Emax, Emin
 real alat
 integer orbnum 
 integer PlotXSFs
 integer PlotPPNs
 integer PlotImag
 integer GetMaxPos
 integer byEnerg
 integer ValIn
 
! ksamples
 real Val
 integer  nksamples,iksample
 real, dimension (:,:), allocatable  :: ksamples 
 real x,y,z



! Procedure
! ===========================================================================
 

  Write (*,*)  " ================= ew2mesh_fourier ================= "
  write (*,*) "pbands ="
  do i = 1, npbands
	   
       write(*,'(i10)',advance='no') pbands (i)
  end do


  open (unit = 27, file = 'kmap.optional', status = 'unknown')
  read (27,*) byEnerg
  read (27,*) Emin, Emax
  read (27,*) PlotXSFs
  read (27,*) PlotPPNs
  read (27,*) PlotImag
  read (27,*) GetMaxPos
  read (27,*) alat
  read (27,*) ValIn
  close (27)

  if (ValIn .eq. 1) then
    open (unit = 28, file = 'ksamples.optional', status = 'unknown')
     read (28,*) nksamples
     allocate(ksamples(3,nksamples))
     do iksample = 1,nksamples
       read (28,*) ksamples(:,iksample)    
     end do   
    close (28)
    open (unit = 1028, file = 'ksamples.dat', status = 'unknown')
  end if




  if (PlotImag .eq. 1) then
     allocate ( ewfIm(0:nrm-1))
     allocate ( ewfRe(0:nrm-1))
  end if 

   Write (*,*) "DEBUG:"
   write (*,'(A,3f16.8)') "a1vec", a1vec
   write (*,'(A,3f16.8)') "a2vec", a2vec
   write (*,'(A,3f16.8)') "a3vec", a3vec

 if (GetMaxPos .gt. 0) then
    open (unit = 33, file = 'kmaxs.dat', status = 'unknown')
 end if

  Write (*,*) "HERE 3"

! reset variables
   job = 0

! allocate aux arrays
   allocate ( ewfaux(0:nrm-1))

   allocate ( positive(rm3,rm2))
   allocate ( negative(rm3,rm2))

! fourier allocation
   ALLOCATE (in3D(rm1,rm2,rm3))
  ! ALLOCATE (four3D(rm1,rm2,rm3))
   ALLOCATE (four3DX(rm1,rm2,rm3))
  ! ALLOCATE (axum3D(rm1,rm2,rm3))

  ! ai = cmplx(0.0d0, 1.0d0)

! set nr(:)
   nr(1) = rm1
   nr(2) = rm2
   nr(3) = rm3



allocate (myBands(norbitals))
myBands=0
if ( byenerg .eq. 1 ) then
    do iband = 1, norbitals
       if( (eigen_k(iband,1) .lt. Emax) .AND. ( eigen_k(iband,1) .gt. Emin) ) myBands(iband)=1
    end do
else
    do i = 1, npbands
       myBands(pbands (i))=1
    end do
end if 

!
Write (*,*) "DEBUG >>>"
    do iband = 1, norbitals
       if( myBands(iband) .eq. 1 ) write (*,*) iband 
    end do
Write (*,*) "<<<< DEBUG"
!


do iband = 1, norbitals
  if ( myBands(iband) .eq. 1) then
       Write (*,*) iband, eigen_k(iband,1) 
       ewfaux = 0.0d0

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

! evaluate the vector between the iatom and the mesh point gP
        do i = 1,3
         dXr(i) = ram2rc(i,imesh) + u1X(i)
        enddo
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
           psi1(imu) = psiL(lmu)*psiR
           imu = imu + 1
          enddo ! do lmu
         enddo ! do issh

! map the point from the extended mesh into the normal mesh
         ind = e2r(index) - 1

         dens = 0.0d0

! Loop over the special k points.
         do ikpoint = 1, nkpoints
! new, with phase
          do imu = 1, num_orb(in1)
           mmu = imu + degelec(iatom)
           dens = dens + bbnkre(mmu,iband,ikpoint)*psi1(imu)
          end do ! do inu

! Finished loop over kpoints
        end do ! do ikpoint

! sum up density at mesh point
        ewfaux(ind) = ewfaux(ind) + dens
!        densmax = max(dens,densmax)

      end do ! do imesh
    end do ! do iatom




 ! ==================== START FOURIER TRANSOFRM ==================

  !  write (*,*) 'DEBUG: Starting FFT'

   ! in3D   = 0.0
   ! four3D = 0.0



   DO k = 1, rm3
      DO j = 1, rm2
         DO i = 1, rm1
            index = i + rm1*(j-1) + rm1*rm2*(k-1) - 1
            in3D(i,j,k) = ewfaux(index)   
         END DO
      END DO
   END DO

 !  CALL dfftw_plan_dft_c2c_3d(plan,rm1,rm2,rm3,in3D,four3D,FFTW_ESTIMATE)  !FFTW_MEASURE)
 
   CALL dfftw_plan_dft_3d(plan,rm3,rm2,rm1,in3D,in3D,FFTW_FORWARD,FFTW_MEASURE)  !FFTW_MEASURE)
   CALL dfftw_execute(plan,in3D,in3D)
   CALL dfftw_destroy_plan(plan)


   DO k = 1, rm3
      DO j = 1, rm2
         DO i = 1, rm1
                four3DX(i,j,k) = in3D( MOD((i+rm1/2),rm1)+1, MOD((j+rm2/2),rm2)+1, MOD((k+rm3/2),rm3)+1 )   
         END DO
      END DO
   END DO

   DO k = 1, rm3
      DO j = 1, rm2
         DO i = 1, rm1
            index = i + rm1*(j-1) + rm1*rm2*(k-1) - 1
            ewfaux(index) = abs(four3DX(i,j,k))   
            ! ewfacum(index) = ewfacum(index) + ewfaux(index)
         END DO
      END DO
   END DO

if (PlotImag .eq. 1) then
   DO k = 1, rm3
      DO j = 1, rm2
         DO i = 1, rm1
            index = i + rm1*(j-1) + rm1*rm2*(k-1) - 1
            ewfIm(index) = Real (four3DX(i,j,k))   
            ewfRe(index) = aImag(four3DX(i,j,k))  
         END DO
      END DO
   END DO
end if


! ======================== END FOURIER ==========================

 ValMax = 0.0
 if (GetMaxPos .gt. 0) then
     DO k = 1, rm3
      DO j = 1, rm2
       DO i = 1, rm1
          if (abs(four3DX(i,j,k)) .gt.  ValMax ) then
            ValMax = abs(four3DX(i,j,k))
            kmax = k
            jmax = j
            imax = i
          end if  
       END DO
      END DO
     END DO
    kmax=kmax-(rm3/2)
    jmax=jmax-(rm2/2)
    imax=imax-(rm1/2)
    rimax = imax * ( alat / a1vec(1) )
    rjmax = jmax * ( alat / a2vec(2) )
    rkmax = kmax * ( alat / a3vec(3) )
   write (33,'(i10,f16.8,3i10,4f16.8)') iband, eigen_k(iband,1), imax,jmax,kmax, ValMax, rimax, rjmax, rkmax
 end if

if (ValIn .eq. 1) then
  write (1028,'(f10.5)', advance='no')   eigen_k(iband,1)
  do iksample = 1,nksamples
          x = rm1*0.5  + ksamples(1,iksample)*a1vec(1)/ alat 
          y = rm2*0.5  + ksamples(2,iksample)*a2vec(2)/ alat
          z = rm3*0.5  + ksamples(3,iksample)*a3vec(3)/ alat
       Val = inter(four3DX,     x,y,z )
       write (1028,'(f16.8)',advance='no')    val
  end do  
  write (1028,*)
end if



if (PlotPPNs .gt. 0) then

! print image
   pmat => ewfaux
   ind = 0
   maxv = 0.0
   positive = 0.0
   negative = 0.0
   write (name,'(i4.4)') iband
   namewf = 'orbfft_'//name//'.ppm'
   Write (*,*) " Write out image", namewf
   open ( unit = 302, file = namewf, status = 'unknown' )
   write (302,'(A2)') "P3"
   write (302,'(1i6,1i6)') ,rm2,rm3
   write (302,*) "255"
   do k = 1, rm3
    do j = 1, rm2
     do i = 1, rm1
       if (pmat(ind) .gt.   0.00) positive(k,j)=max(positive(k,j),pmat(ind))
       if (pmat(ind) .lt. 0.00) negative(k,j)=max(negative(k,j),-pmat(ind))
      ind = ind + 1
     enddo ! i
    enddo ! j
   enddo ! k
   do k = 1, rm3
    do j = 1, rm2
      maxv = max(maxv, positive(k,j))
      maxv = max(maxv, negative(k,j))
    enddo ! j
   enddo ! k
!   write (302,*) maxv
  do k = 1, rm3
    do j = 1, rm2
      cR=int(255*negative(k,j)/maxv)
      cB=int(255*positive(k,j)/maxv)
      ! write (302,*) cR," 0 ",cB
      ! write (302,*) 255-cB,int(max(0.0,float(255-cR-cB))),255-cR
      ! write (302,*) 255-cB,255-int(abs(float(cR-cB))),255-cR
       write (302,*) 255-cB,int(max(0.0,float(255-cR-cB))),255-cR
    enddo ! j
   enddo ! k
   close(302)
endif

if (PlotXSFs .gt. 0) then
Write (*,*) " Write OUT GRID"
     ! file = 100 + iband
     write (name,'(i4.4)') iband
     namewf = 'fourier_'//name//'.xsf'
     Write (*,*) " Write out grid ", namewf
     pmat => ewfaux
     mssg = 'wavefunc_3D'
     call writeout_xsf (namewf, mssg, pmat)

if (PlotImag .eq. 1) then
     write (name,'(i4.4)') iband
     namewf = 'fft_Im_'//name//'.xsf'
     pmat => ewfIm
     mssg = 'wavefunc_3D'
     call writeout_xsf (namewf, mssg, pmat)

     write (name,'(i4.4)') iband
     namewf = 'fft_Re_'//name//'.xsf'
     pmat => ewfRe
     mssg = 'wavefunc_3D'
     call writeout_xsf (namewf, mssg, pmat)

end if
endif

write (*,*) "DEBUG: 3"
       end if ! Emin < eigen_k(iband,1) < Emax
    enddo ! iband
! =============== END OF BIG CYCLE ============

 if (GetMaxPos .gt. 0) then
    close(33)
 end if

if (ValIn .eq. 1) then
  close(1028)
end if


deallocate (myBands)

   deallocate ( in3D)
   ! deallocate ( four3D)
   deallocate ( four3DX)
   deallocate (ewfaux)
  if (PlotImag .eq. 1) then
     deallocate ( ewfIm)
     deallocate ( ewfRe)
  end if 



   deallocate ( positive)
   deallocate ( negative)

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)

   return
 end subroutine ew2mesh_fourier




!===================
!===    UTILS    ===
!===================

real function inter(pvar,x,y,z)
   !real, dimension (:,:,:), pointer,     intent (in) :: pvar
     COMPLEX(8), ALLOCATABLE, DIMENSION(:,:,:),    intent (in) :: pvar
   real, intent (in) :: x
   real, intent (in) :: y
   real, intent (in) :: z
   real val

   COMPLEX(8) CVAL
   real dx,dy,dz,mdx,mdy,mdz
   integer ix,iy,iz
   ix = int(x)
   iy = int(y)
   iz = int(z)
   dx = x - ix
   dy = y - iy
   dz = z - iz
   mdx=1.0-dx
   mdy=1.0-dy
   mdz=1.0-dz
 !  write (*,*),"========="
 !  write (*,*) x,ix,dx
 !  write (*,*) y,iy,dy
 !  write (*,*) z,iz,dz
 CVAL =   pvar(ix+1,iy+1,iz+1)*( dx)*( dy)*( dz)      & 
 &    + pvar(ix+1,iy+1,iz  )*( dx)*( dy)*(mdz)      &
 &    + pvar(ix+1,iy  ,iz+1)*( dx)*(mdy)*( dz)      &
 &    + pvar(ix+1,iy  ,iz  )*( dx)*(mdy)*(mdz)      &
 &    + pvar(ix  ,iy+1,iz+1)*(mdx)*( dy)*( dz)      &
 &    + pvar(ix  ,iy+1,iz  )*(mdx)*( dy)*(mdz)      & 
 &    + pvar(ix  ,iy  ,iz+1)*(mdx)*(mdy)*( dz)      &
 &    + pvar(ix  ,iy  ,iz  )*(mdx)*(mdy)*(mdz)   
  inter = abs(CVAL)   
   return
end

