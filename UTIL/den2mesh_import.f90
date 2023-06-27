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
! ===========================================================================
! Code written by:
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine den2mesh_import (icluster, etrans, iphase)

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use outputs
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in) :: icluster, iphase
   real, intent (in) :: etrans


!Output
! Local Parameters and Data Declaration
! ===========================================================================
   real, parameter ::  Debye = 0.208194
   real, parameter ::  pi = 3.14159265358
   interface
    subroutine writeout_xsf (xsfname, message, aa)
     real, dimension (:), pointer, intent (in) :: aa
     character (len=40) xsfname
     character (len=30) message
    end
   end interface
! Local Variable Declaration and Description
! ===========================================================================
!  For importing rho
   integer iplace, jplace, i1, j1, imax, jmax
   complex*16, dimension (norbitals, norbitals) :: flatrho
!   real, dimension (norbitals, norbitals) :: rhomag
!   real dimension (norbitals, norbitals,2) :: rhophase
   real phasemunu, phase
   real deV
   real maxrho, mingrid, maxgrid
   complex*8 cmaxrho
   integer energyindex
   character (len = 25) estring
   character (len = 25) phstring
! end rho importing variables

   integer iatom
   integer imu, inu,  nnu, mmu
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

   integer iD
   real maxD

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


   real, dimension (3) :: rvec
   real, dimension (3) :: r1
   real, dimension (3) :: r2
   real, dimension (3) :: r21
   real, dimension (3) :: u
   real, dimension (3) :: dXr
   real, dimension (3) :: dYr
   real, dimension (3) :: gX
   real, dimension (3) :: gP
   real, dimension (3) :: X0
   real, dimension (3) :: Y0


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
   real, dimension (3,3)          :: amat

   real, dimension (:), pointer   :: pmat
   character (len=40) filename
   character (len=30) mssg


! Procedure
! ===========================================================================
! Import rho
  write(*,*) '================ Welcome to den2mesh_import ================'

   write(estring,'(f4.2)') etrans
   estring = trim(estring)
!   write(*,*)'estring', estring
   phase = (iphase-1)*22.5
   write(*,'(a6,f6.1)') 'Phase', phase
!write(*,*) 'Phase', phase
!   write(phstring,'(i4.4)') phase
!   phstring = trim(phstring)
   if (iphase .eq. 1) then !read only for first phase
	   write(*,*) 'Reading rho from file rhoin for ', norbitals, 'orbitals'
	    open (unit = 18, file = 'rhoin', status = 'old')
	   read(18,*) deV !separation between energy states recorded
	   energyindex =idnint(etrans/deV + 1)
	   write(*,*) 'Reading rho matrix from rhoin at position', energyindex
	   write(*,'(a29,f5.3)') 'Nearest transition energy at ', (energyindex-1)*deV

	   do j1 = 1,energyindex !last one read is one used
	    do iplace = 1,norbitals
	        read(18, *) (flatrho(iplace,jplace), jplace = iplace,norbitals)
	    end do
	   end do
	   close(18)
	 ! Fill out lower diagonal.  Rho for density use is symmetric, not hermitian
	 ! The phase is not a quantum phase here, but phase delay in oscillation
      maxrho = 0.0d0
	   do iplace = 1,norbitals
	    do jplace = iplace,norbitals
	      flatrho(jplace,iplace) = flatrho(iplace,jplace)
	      if ( abs(flatrho(iplace,jplace)) .gt. maxrho) then
	        maxrho = abs(flatrho(iplace,jplace))
	        cmaxrho = flatrho(iplace,jplace)
	        imax = iplace
	        jmax = jplace
	      endif
	    end do
	   end do
!	   write(*,'(f6.3)') 'Maximum rho magnitude', maxrho
!	   write(*,*) 'Maximum rho magnitude', imax,jmax, maxrho, cmaxrho
   	   write(*,*) 'Maximum rho magnitude', maxrho, cmaxrho

   end if !phase = 000
! Convert to 4-index rho
   rho = 0.0d0
   maxrho = 0.0d0
   do iatom = 1, natoms
     in1 = imass(iatom)
     do ineigh = 1, neighn(iatom)
       mbeta = neigh_b(ineigh,iatom)
       jatom = neigh_j(ineigh,iatom)
       in2 = imass(jatom)
       do imu = 1, num_orb(in1)
          mmu = imu + degelec(iatom)
          do inu = 1, num_orb(in2)
!          write(*,*) 'imu,inu,ineigh,iatom'
!           write(*,*)imu,inu,ineigh,iatom
            nnu = inu + degelec(jatom)
!          write(*,*) 'mmu,nnu', mmu,nnu
!           write(*,*)'flatrho', abs(flatrho(mmu,nnu))
            phasemunu = phase*pi/180.0 +  	&
               & atan2(aimag(flatrho(mmu,nnu)),real(flatrho(mmu,nnu)))
!               write(*,*) 'phase', phasemunu
            rho(imu,inu,ineigh,iatom) = abs(flatrho(mmu,nnu))*cos(phasemunu)
!	        write(*,*) 'rho',rho(imu,inu,ineigh,iatom)
	      if ( abs(rho(imu,inu,ineigh,iatom)) .gt. maxrho) maxrho = rho(imu,inu,ineigh,iatom)
          end do
      end do
   end do
  end do
  write(*,'(a29,f6.1)') 'Maximum phased rho magnitude', maxrho

! Grid calculation begins
! reset variables
   drhoG = 0.0d0
   job = 0

! set nr(:)
   nr(1) = rm1
   nr(2) = rm2
   nr(3) = rm3

   renorm = 0.0d0
! Loop over atoms
   do iatom = 1, natoms
! make a copy of the elem grid lattice vector
    do i = 1,3
     amat(:,i) = elvec(i,:)
    enddo
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)

! vector between the atom and initial point of the grid
    do i = 1,3
     rvec(i) = ratom(i,iatom) - g0(i)
! save vector for next
     u(i) = rvec(i)
    enddo ! i
! solve linear set of equations (3x3)
    call dgesl(amat,3,ipiv,u,job)
! round coefficients to get the position of the nearest point of the grid gX
    i0 = nint( u(1) )
    j0 = nint( u(2) )
    k0 = nint( u(3) )

! find the vector u between the iatom X0 and the nearest point gX
    u(1) = u(1) - real(i0)
    u(2) = u(2) - real(j0)
    u(3) = u(3) - real(k0)

! check if the nearest grid point is located within the unit cell of the grid
! if not, let map it there
!i0
    if (i0 .lt. 0) then
     i0 = i0 + rm1*(int(abs(i0/rm1)) + 1)
    else
     i0 = i0 - rm1*int(i0/rm1)
    endif
!j0
    if (j0 .lt. 0) then
     j0 = j0 + rm2*(int(abs(j0/rm2)) + 1)
    else
     j0 = j0 - rm2*int(j0/rm2)
    endif
!k0
    if (k0 .lt. 0) then
     k0 = k0 + rm3*(int(abs(k0/rm3)) + 1)
    else
     k0 = k0 - rm3*int(k0/rm3)
    endif


! find the coordinates of gX
    gX(1) = i0*elvec(1,1) + j0*elvec(2,1) + k0*elvec(3,1)
    gX(2) = i0*elvec(1,2) + j0*elvec(2,2) + k0*elvec(3,2)
    gX(3) = i0*elvec(1,3) + j0*elvec(2,3) + k0*elvec(3,3)


! evaluate coordinates of iatom in the grid coord system
    X0(1) = gX(1) + u(1)*elvec(1,1) + u(2)*elvec(2,1) + u(3)*elvec(3,1)
    X0(2) = gX(2) + u(1)*elvec(1,2) + u(2)*elvec(2,2) + u(3)*elvec(3,2)
    X0(3) = gX(3) + u(1)*elvec(1,3) + u(2)*elvec(2,3) + u(3)*elvec(3,3)

! vector pointing from gX to X0
    u(1) = X0(1) - gX(1)
    u(2) = X0(2) - gX(2)
    u(3) = X0(3) - gX(3)

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
        dXr(i) = ram2rc(i,imesh) - u(i)
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
          psi2(imu) = psiL(lmu)*psiR
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
  write (*,*) ' -- Total charge on grid =',dens

! calc dipole with the unit cell
! bch also find maximum and minimum on grid
  index = 0
  maxgrid = 0.0
  mingrid = 0.0
  dip_x =0.0d0
  dip_y =0.0d0
  dip_z =0.0d0
  do k = 0, rm3-1
   do j = 0, rm2-1
    do i = 0, rm1-1
     dqi = drhoG(index)
     if (dqi .gt. maxgrid) maxgrid = dqi
     if (dqi .lt. mingrid) mingrid = dqi
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
  write (*,*) 'Maximum density', maxgrid
  write (*,*) 'Minimum density', mingrid


! -----------
! dipole units (charge x distance)
! 1 Debye = 0.208194 eAng = 0.393430 eBohr
!

! reminder
! drhoG ... variation of the SCF density (rho - rho0)
! rhoG0 ... the neutral atom density

! write out xsf files
! total density
  allocate (rhotmp (0:nrm-1))
! get imported density
  rhotmp = drhoG !imported density gets written to 'density.xsf'

! write out rho into xsf file
  pmat => rhotmp
!  write(filename,'(a7,a,a4,a,a3,a4)') 'density','_',estring,'_',phstring,'.xsf'
  write(filename,'(a7,a,a4,a,i2.2,a4)') 'density','_',estring,'_',iphase,'.xsf'
  mssg = 'density_3D'
  call writeout_xsf (filename, mssg, pmat)
  deallocate (rhotmp)





! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)
301 format (2x,'Dipole_x =',e14.6,'  [D] ')
302 format (2x,'Dipole_y =',e14.6,'  [D] ')
303 format (2x,'Dipole_z =',e14.6,'  [D] ')
304 format (2x,'Dipole_tot =',e14.6,'  [D] ')
600 format (10000F11.6)
   return
 end subroutine den2mesh_import

