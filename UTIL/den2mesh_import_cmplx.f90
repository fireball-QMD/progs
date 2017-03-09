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


! den2mesh_import_cmplx.f90
! Program Description
! ===========================================================================
!       Project density on the mesh. Density read from file, and sum over
! energies from etrans1 to etrans2.
!
! ===========================================================================
! Code written by:
! ===========================================================================
!
! ===========================================================================
 subroutine den2mesh_import_cmplx (icluster, etrans1, etrans2, projectyn)

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
   real, intent (in) :: etrans1, etrans2
   character (len = 1), intent(in) :: projectyn


!Output
! Local Parameters and Data Declaration
! ===========================================================================
   real, parameter ::  Debye = 0.208194
   real, parameter ::  pi = 3.14159265358
   real, parameter ::  gridtol = 1.0d-2
   complex, dimension(:), allocatable :: drhoCplx
   complex, dimension(:,:,:,:), allocatable :: rhoc
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
   integer iplace, jplace, i1, j1, imax, jmax, idummy, iw, n1, indexmax
   real sum1, phasemax
   complex ai
   complex, dimension (:,:) , allocatable :: flatrho
   real, dimension (:) , allocatable :: temp
   complex, dimension (:) , allocatable :: rhow
   real, dimension (:,:), allocatable  :: rhotemp


   real deV
   real maxrho, mingrid, maxgrid, temp1,temp2
   complex cmaxrho
   integer energyindex1, energyindex2, ienergy, nenergy_max, inonzero, norb
   character (len = 25) estring1, estring2, estring
   !Vasp
   integer, dimension (3) :: ngridvasp
   real, dimension (3,3) :: lvecvasp
   integer natomvasp
   real, dimension(:), allocatable :: rhovasp
   integer irow, nrows

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
   complex dens


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
   character (len=50) filename
   character (len=30) mssg


! Procedure
! ===========================================================================
! Import rho
 allocate( flatrho(norbitals,norbitals) )
 write(*,*) '================ Welcome to den2mesh_import ================'
 if (etrans1 .lt. 10.0) then
      write(estring1,'(f4.2)') etrans1
      else
      write(estring1,'(f5.2)') etrans1
 end if
 if (etrans2 .lt. 10.0) then
      write(estring2,'(f4.2)') etrans2
      else
      write(estring2,'(f5.2)') etrans2
 end if
 estring = trim(estring1)//'-'//trim(estring2)
 estring = trim(estring)
! write(*,*) 'estring', estring
 if (etrans1 .eq. 0.0) then !import flatrho directly from NxN file
   allocate( rhotemp(norbitals,norbitals) )
   write(*,*) 'File name?'
   read(*,*) filename
   open (unit = 18, file = filename, status = 'old')
   do iplace = 1,norbitals
        read(18,*) (rhotemp(iplace, jplace), jplace = 1, norbitals)! will fill only real part
   end do
   close(18)
	flatrho = cmplx(rhotemp, 0.0d0)
   deallocate( rhotemp )
  else if (etrans1 .gt. 0.0) then
   flatrho = cmplx(0.0d0, 0.0d0)
   write(*,*) 'Reading rho from file rhoin for ', norbitals, 'orbitals'
   open (unit = 18, file = 'rhoin', status = 'old')
   read(18,*) !'norbitals         energy step        skip_factor        nenergy'
   read(18,*) norb, deV, idummy, nenergy_max   !dev: separation between energy states recorded in FT
!   write(*,*) norb, deV, idummy, nenergy_max   !dev: separation between energy states recorded in FT

   if (norb .ne. norbitals) then
   	write(*,*) 'Number of orbitals does not match file in den2mesh_import.xxx ', norbitals, norb
   	return
   end if
   allocate( rhow(nenergy_max) )
      allocate( temp(2*nenergy_max) )
   energyindex1 =idnint(etrans1/deV + 1)
   energyindex2 =idnint(etrans2/deV + 1)
   write(*,*) 'Reading rho matrix from rhoin, positions', energyindex1, 'to', energyindex2
   write(*,'(a33,f5.2,a4,f5.2)')'Summing over transition energies ', (energyindex1-1)*deV , ' to ', (energyindex2-1)*deV

   if (energyindex2 .gt. nenergy_max) then
   	write(*,*) 'Transition position number too large for rhoin '
   	return
   end if
   do iplace = 1,norbitals
   	  if (mod(iplace,100) .eq. 0 ) then
	    write(*,*), 'Reading orbital', iplace,'of', norbitals
	  end if

      do jplace =1, iplace !lower diagonal is stored
        read(18,*) !iplace, jplace
        read(18,*) (temp(iw), iw = 1, 2*energyindex2)
        do iw =  1, energyindex2
           rhow(iw) = cmplx(temp(2*iw-1),temp(2*iw))
        end do
!        write(*,*) iplace, jplace
!   		write(*,*) (rhow(iw), iw = 1, energyindex2)
        flatrho(iplace, jplace) = flatrho(iplace, jplace)  &
         & + sum(rhow(energyindex1:energyindex2))
!   	  write(*,*) 'sum',  sum(rhow(energyindex1:energyindex2))
   	  end do
   end do
   close(18)

	 ! Fill out upper diagonal.  Rho for density use is symmetric, not hermitian
	 ! The phase is not a quantum phase here, but phase delay in oscillation
   maxrho = 0.0d0
   do iplace = 1,norbitals
	    do jplace = iplace+1,norbitals
	      flatrho(iplace,jplace) = flatrho(jplace,iplace)
	      if ( abs(flatrho(iplace,jplace)) .gt. maxrho) then
	        maxrho = abs(flatrho(iplace,jplace))
	        cmaxrho = flatrho(iplace,jplace)
	        imax = iplace
	        jmax = jplace
	      endif
	    end do
   end do

   deallocate (rhow)
   deallocate (temp)
 else !etrans1 then must be less than 0.0; read in grid from VASP file
! 	   allocate( rhotemp(norbitals,norbitals) )
   write(*,*) 'Reading from CHG Vasp file'
   open (unit = 18, file = 'CHG', status = 'old')
   read(18,*) !Comment
   read(18,*) !Scale factor for unit cell
   do iplace = 1,3 !read in lattice vectors
    read(18,*) (lvecvasp(iplace,jplace), jplace = 1, 3)
   end do
   read(18,*) natomvasp !works only for clusters of one atom type (Vasp has several entries for several atoms)
   read(18,*) !direct
   !read atomic position
   do iplace = 1,natomvasp
         read(18,*) !read past atomic positions
   end do
   read(18,*) !comment
   read(18,*) (ngridvasp(iplace), iplace = 1, 3)
   rm1 = ngridvasp(1)
   rm2 = ngridvasp(2)
   rm3 = ngridvasp(3)
   nrm = rm1 * rm2 * rm3
   write(*,*)'Grid points', nrm
   allocate (drhoCplx (0:nrm-1))
   allocate (rhovasp (0:nrm-1))
   nrows = int(nrm/10)
   write(*,*)'Number of rows', nrows
   do irow = 1,nrows+1
!     write(*,*), 'irow', irow
   	 if (irow .lt. nrows+1) then
      if (mod(iplace,10000) .eq. 0 ) then
	    write(*,*), 'Reading grid point', irow * 10
	  end if
!	  read(18,'(E12.5$)')( rhovasp((irow-1)*10+iplace-1),iplace = 1,10)
	  read(18,*)( rhovasp((irow-1)*10+iplace-1),iplace = 1,10)
	 else ! if (irow .eq. nrows+1
	  read(18,*)( rhovasp((irow-1)*10+iplace-1),iplace = 1,mod(nrm,10))
	 end if !irow
   end do
   close(18)
   drhoCplx = cmplx(rhovasp, 0.0d0)
   deallocate(rhovasp)
   !grid lattice:
   do iplace = 1,3
   elvec(:,iplace) = lvecvasp(:,iplace) / ngridvasp(iplace)
!   write(*,*) 'elvec', iplace
!   write(*,*)  elvec(:,iplace)
!  write(*,*) 'ngridvasp',  ngridvasp(iplace)
   end do
   go to 1000

 end if !etrans1 .eq. 0.0
	   ! Note that off diagonal rho elements alone have little meaning in nonorthogonal basis...must
	   ! be multiplied by smat to judge true density contribution.
!	   write(*,*) 'Maximum rho magnitude read in', imax,jmax, maxrho, cmaxrho
!   	   write(*,*) 'Maximum rho magnitude', maxrho, cmaxrho
!'test

!	     open (unit = 18, file = 'flatrho.out', status = 'unknown')
!	     do iplace = 1,norbitals
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace = 1, norbitals)!
!	     end do
!	     close(18)

!'test

!		open (unit = 18, file = 'flatrho.out', status = 'unknown')
!	     do iplace = 1,norbitals
!	     	do jplace = 1,iplace
!	     	 write (18,*), iplace, jplace
!	     	 write (18,600), iplace, flatrho(iplace, jplace)
!	        end do
!	     end do
!	     close(18)
   ! do we project out only a portion of rho?
!	open (unit = 18, file = 'junk1.out', status = 'unknown')
!		 n1 = 2
!		 write(18,*)  n1,n1+2
!		 sum1 = 0.0
!	     do iplace = n1,n1+2
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace =  n1,n1+2)!
!	        sum1 = sum1 + real(flatrho(iplace, iplace))
!	     end do
!	     write(18,*) sum1
!		 write(18,*)
!
!		 n1 = 6
!		 write(18,*)  n1,n1+2
!		 sum1 = 0.0
!	     do iplace = n1,n1+2
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace =  n1,n1+2)!
!	        sum1 = sum1 + real(flatrho(iplace, iplace))
!	     end do
!	     write(18,*) sum1
!		 write(18,*)
!
!		 n1 = 102
!		 write(18,*)  n1,n1+2
!		 sum1 = 0.0
!	     do iplace = n1,n1+2
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace =  n1,n1+2)!
!	        sum1 = sum1 + real(flatrho(iplace, iplace))
!	     end do
!	     write(18,*) sum1
!		 write(18,*)
!
!		 n1 = 118
!		 write(18,*)  n1,n1+2
!		 sum1 = 0.0
!	     do iplace = n1,n1+2
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace =  n1,n1+2)!
!	        sum1 = sum1 + real(flatrho(iplace, iplace))
!	     end do
!	     write(18,*) sum1
!		 write(18,*)
!	     close(18)
!
!   if (projectyn .eq. 'y') then
!       call  projectrho (flatrho, norbitals)
!'test
!	     open (unit = 18, file = 'flatrho2.out', status = 'unknown')
!	     do iplace = 1,norbitals
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace = 1, norbitals)!
!	     end do
!	     close(18)
!
!		 open (unit = 18, file = 'junk2.out', status = 'unknown')
!		 n1 = 2
!		 write(18,*)  n1,n1+2
!		 sum1 = 0.0
!	     do iplace = n1,n1+2
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace =  n1,n1+2)!
!	        sum1 = sum1 + real(flatrho(iplace, iplace))
!	     end do
!	     write(18,*) sum1
!		 write(18,*)
!
!		 n1 = 6
!		 write(18,*)  n1,n1+2
!		 sum1 = 0.0
!	     do iplace = n1,n1+2
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace =  n1,n1+2)!
!	        sum1 = sum1 + real(flatrho(iplace, iplace))
!	     end do
!	     write(18,*) sum1
!		 write(18,*)
!
!		 n1 = 102
!		 write(18,*)  n1,n1+2
!		 sum1 = 0.0
!	     do iplace = n1,n1+2
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace =  n1,n1+2)!
!	        sum1 = sum1 + real(flatrho(iplace, iplace))
!	     end do
!	     write(18,*) sum1
!		 write(18,*)
!
!		 n1 = 118
!		 write(18,*)  n1,n1+2
!		 sum1 = 0.0
!	     do iplace = n1,n1+2
!	        write(18,600) (real(flatrho(iplace, jplace)), jplace =  n1,n1+2)!
!	        sum1 = sum1 + real(flatrho(iplace, iplace))
!	     end do
!	     write(18,*) sum1
!		 write(18,*)
!	     close(18)
!
!	end if

! Convert to 4-index complex rho
   allocate (drhoCplx (0:nrm-1))
   allocate (rhoc (numorb_max, numorb_max, neigh_max, natoms))
   rhoc = 0.0d0
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
            nnu = inu + degelec(jatom)
             rhoc(imu,inu,ineigh,iatom) = flatrho(mmu,nnu)
	      if ( abs(rhoc(imu,inu,ineigh,iatom)) .gt. maxrho) maxrho = abs(rhoc(imu,inu,ineigh,iatom))

          end do
      end do
    end do
   end do
!   write(*,'(a33,f6.1)') ' Maximum magnitude in rhoc', maxrho

! Grid calculation begins
! reset variables
   drhoCplx = 0.0d0
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
          dens = dens + rhoc(inu,imu,ineigh,iatom)                &
      &        *psi1(inu)*psi2(imu)
         enddo ! do inu
        enddo ! do imu

! store variation of density at given point
        drhoCplx(ind) = dens + drhoCplx(ind)
       endif ! if (Rc_max)
      end do ! do imesh
    end do ! do ineigh
   end do ! do iatom

! test for hydrogen s-orbital the quality of the integration;
! the ratio it should go to one with higher Ecut
   dens = 4.0d0*3.141592653589793238462643*Rc_max**3/3.0d0

1000 continue
! calc the total denstity
  dens = 0.0d0
  do i = 0,nrm-1
   dens = dens + drhoCplx(i)*dvol
  enddo
  write (*,*) ' -- Total charge on grid =',dens
  write (*,*) ' -- Number of grid points ', (rm3*rm2*rm1)
! Write grid to file
! bch also find maximum and minimum on grid
  inonzero = 0
  index = 0
  maxgrid = 0.0
  mingrid = 0.0
  do k = 0, rm3-1
   do j = 0, rm2-1
    do i = 0, rm1-1
     dqi = abs(drhoCplx(index))
     if (dqi .gt. maxgrid) then
       maxgrid = dqi
       indexmax = index
     end if !dqi
     if (dqi .lt. mingrid) mingrid = dqi

!     if (dqi .gt. gridtol * maxgrid) inonzero = inonzero + 1
     index = index + 1
    enddo ! i
   enddo ! j
  enddo ! k

!  write (*,*) 'Maximum density magnitude', maxgrid
!  write (*,*) 'Minimum density magnitude', mingrid
!  write (*,*) 'Number of nonzero grid points', inonzero


! write out rho into xsf file
  ! create real density at phase of maximum amplitude
  allocate (rhotmp (0:nrm-1))
  phasemax = atan(aimag(drhoCplx(indexmax))/real(drhoCplx(indexmax)))
  ai = cmplx(0.0d0,1.0d0)
  index = 0
  do k = 0, rm3-1
   do j = 0, rm2-1
    do i = 0, rm1-1
!    real(exp(i*phasemax)*denout)
     rhotmp(index) = real(exp(ai*phasemax) * drhoCplx(index))
     index = index + 1
    enddo ! i
   enddo ! j
  enddo ! k
  !
  pmat => rhotmp
!  write(filename,'(a7,a,a4,a,a3,a4)') 'density','_',estring,'_',phstring,'.xsf'
  write(filename,'(a7,a,a4,a4)') 'density','_',estring,'.xsf'
  mssg = 'density_3D'
  call writeout_xsf (filename, mssg, pmat)
  deallocate (rhotmp)


! simple grid and geodesic
  filename ='density'//'_'//trim(estring)//'.gridc'
  call writeout_drhoc (filename,drhoCplx, maxgrid, gridtol)
  filename ='density'//'_'//trim(estring)//'.geo'
!  !if (norbitals .ge. 240) call writeout_drhoc_geo (filename,drhoCplx, maxgrid, gridtol)
   call writeout_drhoc_geo (filename,drhoCplx, maxgrid, gridtol)

  if (etrans1 .ge. 0.0) then
   deallocate (drhoCplx)
   deallocate (rhoc)
   deallocate (flatrho)
  endif



	write(*,*) 'Exiting den2mesh_import_cmplx'


! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)
301 format (2x,'Dipole_x =',e14.6,'  [D] ')
302 format (2x,'Dipole_y =',e14.6,'  [D] ')
303 format (2x,'Dipole_z =',e14.6,'  [D] ')
304 format (2x,'Dipole_tot =',e14.6,'  [D] ')
400 format (10000F8.3)
600  	format (10000F22.14)

   return
 end subroutine den2mesh_import_cmplx

