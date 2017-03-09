! copyright info:
!
!                             @Copyright 2002
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega

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


! laplace_fdm.f90
! Program Description
! ===========================================================================
!       Solve Laplace equation using finite difference method.
!
! ===========================================================================
! Code written by:
! ===========================================================================
!
! Program Declaration
! ===========================================================================
! subroutine laplace_fdm (uHdcc)
 subroutine laplace_fdm (icluster)

   use grid
   use charges
   use interactions
   use constants_fireball
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in) :: icluster

! Output
!   real, intent (inout) :: uHdcc 

! Local Parameters and Data Declaration
! ===========================================================================
   real, parameter :: fourpi = -1.0d0/(4.0d0*3.141592653589793238462643)
   real, parameter :: abohr3 = abohr**3.0d0 
   real, parameter :: eq = -1.0d0*sqrt(eq2)
! Local Variable Declaration and Description
! ===========================================================================

   integer i
   integer j
   integer jext
   integer ineigh
   integer i0, j0, k0
   integer iatom

   integer k,index
   
   real dens
   real vpot
   real, allocatable, dimension(:) :: vpot_old


! Procedure
! ===========================================================================
! Initialize
   allocate (vpot_old (0:nrm-1))
! copy hartree potential & convert it into a.u.   
   vpot_old(:) = vcaG(:)
   
! the density is in 1/Ang^3 units


! Initialize variables
   do i = 0, (nrm-1)     
! add the density term & convert it into into a.u.
     vpot = fourpi*eq*drhoG(i) 
     do ineigh = 1,nneighij-1
       jext = n2e(i) + neighij(ineigh)
       j = e2n(jext)
!       write (*,200) i,ineigh,neighij(ineigh),n2e(i),jext,j
       vpot = vpot + vpot_old(j)*d2f(ineigh)
     enddo ! do ineigh
! multiply by the factor and convert it from Hartree into eV
     vcaG(i) = vpot/d2f(nneighij)

   enddo ! do i

! deallocate
   deallocate (vpot_old)


! write out  hartree potential into dv.dat file
   open ( unit = 300, file = 'fdmpot.xsf', status = 'unknown' )
! print the list of atoms
  if (icluster .eq. 1) then

   write (300,*) 'ATOMS'
   do iatom = 1,natoms
    write (300,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
   enddo

  else

   write (300,*) 'CRYSTAL'
   write (300,*) 'PRIMVEC'
   write (300,*) (a1vec(i),i=1,3)
   write (300,*) (a2vec(i),i=1,3)
   write (300,*) (a3vec(i),i=1,3)

   write (300,*) 'PRIMCOORD'
   write (300,*) natoms,1
   do iatom = 1,natoms
    write (300,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
   enddo

  endif

  write (300,*)
  write (300,*) 'BEGIN_BLOCK_DATAGRID_3D'
  write (300,*) 'density_3D'
  write (300,*) 'DATAGRID_3D_DENSITY'
  write (300,*) rm1+1, rm2+1, rm3+1
! print origin of the grid
  write (300,*) 0.0d0, 0.0d0, 0.0d0
! print lattice vector
  write (300,*) (a1vec(i),i=1,3)
  write (300,*) (a2vec(i),i=1,3)
  write (300,*) (a3vec(i),i=1,3)

! print values of the grid point
  do k = 0, rm3
    if (k .eq. rm3) then
     k0 = 0
    else
     k0 = k
    endif
    do j = 0, rm2
      if (j .eq. rm2) then
       j0 = 0
      else
       j0 = j
      endif
      do i = 0, rm1
        if (i .eq. rm1) then
         i0 = 0
        else
         i0 = i
        endif
! mapping index within the regular mesh
        index = i0 + rm1*j0 + rm1*rm2*k0
        write (300,300) vcaG(index)
      enddo ! do i
    enddo ! do j
  enddo ! do k
  write (300,*) 'END_DATAGRID_3D'
  write (300,*) 'END_BLOCK_DATAGRID_3D'

! close file
  close (300)


! Format Statements
! ===========================================================================
100 format (2x, 'Hartree Energy = ',f14.7,' [eV]')
200  format('index = ',6i7)
300 format (e14.6)

   return
 end subroutine laplace_fdm
