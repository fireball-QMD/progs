! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek
!
! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! writeout_drhoc_geo.f90
! Program Description
! ===========================================================================
!       The subroutine projects complex grid densities onto geodesic verticies
! to represent a spherical density (for C60).
!
! ===========================================================================
! Code written by:
! B. Hess
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine writeout_drhoc_geo (filename, dgridc, maxgrid, gridtol)

   use configuration
   use grid
   use charges
   use interactions
   use options
   use constants_fireball
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

   complex,   dimension (0:nrm-1), intent (in) :: dgridc
!   character (len=40), intent(in) :: filename
   character (len=50) :: filename
   real :: maxgrid, gridtol
!, intent(in)

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer i
   integer j
   integer k, kv
   integer i0
   integer j0
   integer k0
   integer index, inonzero
   real xx, yy, zz
   real xavg, yavg, zavg

	integer Nvert, Ntri, header

	real, dimension(:,:), allocatable :: rvert
	complex, dimension(:), allocatable :: denout, denin
	integer, dimension(:,:), allocatable :: itri

	real angwidth, cutoff, radius, rnorm, dalpha
	complex dengridout, dengridin, denvertout, denvertin

	complex, dimension(:), allocatable :: denintemp, denouttemp
	real, dimension(3,3) :: eps2
	real theta, dalphamin
	real, dimension(:,:), allocatable :: rrotvert
	integer kvmin, jv

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================


   write (*,*) '  Project complex imported density onto a sphere',filename
   write (*,100)


! open geodesic file
!     close all
   open (unit = 18, file = 'geo.out', status = 'old')
   read(18,*) !header line
   read(18,*) Nvert, Ntri !Number geodesic vertices
   allocate(rvert(Nvert,3))
   allocate(denin(Nvert), denout(Nvert))
   allocate(itri(Ntri,4))
!  read vertices
   do i = 1,Nvert !
	   read(18, *) (rvert(i,j), j = 1,3)
   end do

!  read triangles.  Note there is first a 3 then the indices
   do i = 1,Ntri !
	   read(18, *) (itri(i,j), j = 1,4)
   end do
   close(18)
!Find center of grid to subtract off.
   xavg = 0.0d0
   yavg = 0.0d0
   zavg = 0.0d0
   index = 0
   do k = 0, rm3-1
    do j = 0, rm2-1
     do i = 0, rm1-1
      index = index + 1
      xx = i*elvec(1,1) + j*elvec(2,1) + k*elvec(3,1)
      yy = i*elvec(1,2) + j*elvec(2,2) + k*elvec(3,2)
      zz = i*elvec(1,3) + j*elvec(2,3) + k*elvec(3,3)
   	  xavg = xavg + xx
   	  yavg = yavg + yy
   	  zavg = zavg + zz
     enddo ! i
    enddo ! j
  enddo ! k
  xavg = xavg/index
  yavg = yavg/index
  zavg = zavg/index

! Project onto sphere: geodesic verticies.  Use a gaussian spread, based on how far
! apart the unit vectors are.
   write(*,*) 'Projecting onto sphere, i.e. geodesic with', Nvert, 'vertices'
!   angwidth =0.02* 3.1459 ! width for gaussian point spread
   angwidth =0.04* 3.1459 ! width for gaussian point spread
   cutoff = angwidth*3 !
   radius = 6.83/2 !radius of C60
   index = 0
   denout = 0.0d0
   denin = 0.0d0
   dengridout = 0.0d0
   dengridin = 0.0d0
   do k = 0, rm3-1
    do j = 0, rm2-1
     do i = 0, rm1-1
     if (abs(dgridc(index)) .gt. gridtol * maxgrid) then
!      xx = i*elvec(1,1) + j*elvec(2,1) + k*elvec(3,1) - xavg
!      yy = i*elvec(1,2) + j*elvec(2,2) + k*elvec(3,2) - yavg
!      zz = i*elvec(1,3) + j*elvec(2,3) + k*elvec(3,3) - zavg
      xx = i*elvec(1,1) + j*elvec(2,1) + k*elvec(3,1) - xavg - elvec(1,1)/2.0
      yy = i*elvec(1,2) + j*elvec(2,2) + k*elvec(3,2) - yavg - elvec(2,2)/2.0
      zz = i*elvec(1,3) + j*elvec(2,3) + k*elvec(3,3) - zavg - elvec(3,3)/2.0
      rnorm = sqrt(xx**2+yy**2+zz**2) !normalize grid point vector to unit sphere.
!  Do some counting on half the sphere (full sphere charge will be approx zero, so not useful)
      if (rnorm .ge. radius .and. zz .ge. 0.0 ) then !came from point outside sphere. Check inner sum
	    dengridout =  dengridout + dgridc(index)
	  elseif (rnorm .lt. radius .and. zz .ge. 0.0 ) then !came from point inside sphere
	    dengridin =  dengridin + dgridc(index)
	  end if
	  xx = xx/rnorm
	  yy = yy/rnorm
	  zz = zz/rnorm
      do kv = 1,Nvert !dalpha is approximate angular difference between grid and vertex vector
        dalpha = sqrt((xx-rvert(kv,1))**2+(yy-rvert(kv,2))**2+(zz-rvert(kv,3))**2)
        if (dalpha .lt. cutoff) then
           if (rnorm .ge. radius) then !came from point outside sphere
            denout(kv) =  denout(kv) + dgridc(index)*exp(-(dalpha/angwidth)**2) !put angular-gaussian point spread
		   else !came from point inside sphere
            denin(kv) =  denin(kv) + dgridc(index)*exp(-(dalpha/angwidth)**2)
		   end if
        end if !dalpha
      end do  ! kv
      end if !abs(dgridc(index))
      index = index + 1
      if (mod(index,100000) .eq. 0 ) then
	    write(*,*), 'Projecting grid point', index
	  end if

     enddo ! i
    enddo ! j
  enddo ! k
  write(*,*) 'Total grid outer charge, for z>0',dengridout
  write(*,*) 'Total grid inner charge, for z>0',dengridin
!
!=============== symmetrize ============
!   allocate(rrotvert(Nvert,3))
!   allocate(denintemp(Nvert),denouttemp(Nvert) )
!! average out small numerical noise and instabilities
!! rotate positions of all vertices by 2pi/5 about z axis, and add to form an average
!		! z rotation by 2pi/5 is our symmetry operation
!		theta = 2.0d0*pi/5.0d0
!		!define epsilon
!		eps2 = 0.0d0
!		eps2(1,1) = cos(theta)
!		eps2(1,2) = -sin(theta)
!		eps2(2,1) = sin(theta)
!		eps2(2,2) = cos(theta)
!		eps2(3,3) = 1.0d0
!
!      denintemp = cmplx(0.0d0,0.0d0)
!      denouttemp = cmplx(0.0d0,0.0d0)
!      do jv = 1,Nvert
!       if (mod(jv,1000) .eq. 0 ) then
!	    write(*,*), 'Averging rotated vertex point', jv
!	   end if
!       !rotate all vertex point by 2pi/5, z axis
!        rrotvert(jv,1) = eps2(1,1)*rrotvert(jv,1) + eps2(1,2)*rrotvert(jv,2)
!		rrotvert(jv,2) = eps2(2,1)*rrotvert(jv,1) + eps2(2,2)*rrotvert(jv,2)
!		rrotvert(jv,3) = rvert(jv,3)
!        !find old vertex closest to
!        dalphamin = 2*pi !Make it way too big
!       do kv = 1,Nvert !dalpha is approximate angular difference between grid and vertex vector
!        dalpha = sqrt((rrotvert(jv,1)-rvert(kv,1))**2+(rrotvert(jv,2)-rvert(kv,2))**2+(rrotvert(jv,3)-rvert(kv,3))**2)
!        if (dalpha .lt. dalphamin) then
!          dalphamin = dalpha
!          kvmin = kv
!        end if
!       end do !kv
!      !Add to vertex point
!      denouttemp(jv) = denouttemp(jv) + (denout(kvmin) + denout(jv))/2.0d0 ! rotated and nonrotated are averaged
!      denintemp(jv) = denintemp(jv) + (denin(kvmin) + denin(jv))/2.0d0
!     end do !jv
!     denout = denouttemp
!     denin = denintemp
!    deallocate(rrotvert)
!    deallocate(denintemp,denouttemp)

     !============== done with rotation

  ! Add up all the charge
  denvertout = 0.0d0
  denvertin = 0.0d0
  do kv = 1,Nvert
   if (rvert(kv,3) .ge. 0) then
    denvertout = denvertout + denout(kv)
    denvertin = denvertin + denin(kv)
   end if
  end do
  write(*,*) 'Total vertex outer charge',denvertout
  write(*,*) 'Total vertex inner charge',denvertin

! Write density on geodesic to file
! Make first part of file identical to geo.out file
  write(*,*) 'Writing to file'
  header = 0
  open ( unit = 302, file = filename, status = 'unknown' )
  write (302,*) header
  write (302,*) Nvert, Ntri, header !Number geodesic vertices
   do i = 1,Nvert !
	   write (302,*) (rvert(i,j), j = 1,3)
   end do
! Write triangles
   do i = 1,Ntri ! 4 for "3" then indices
	   write (302,*) (itri(i,j), j = 1,4)
   end do
! Write complex density on each vertex: first real, then imaginary outside then inside
   do i = 1,Nvert !
	   write (302,'(4f18.7)') real(denout(i)), aimag(denout(i)), real(denin(i)), aimag(denin(i))
   end do

! Write atomic coordinates
   write (302,*) 'ATOMS'
   do iatom = 1,natoms
!   write (302,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
    write (302,'(i2,3f14.8)') nzx(imass(iatom)),(ratom(i,iatom),i=1,3)
   enddo
! close file
  close (302)
! Deallocate Arrays
  deallocate( rvert, itri, denin, denout)
! ===========================================================================

! Format Statements
! ===========================================================================
100     format (70('='))
200     format (e14.6)
700		format (10000f14.8)

        return
      end subroutine writeout_drhoc_geo
