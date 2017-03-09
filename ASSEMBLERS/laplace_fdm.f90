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
 subroutine laplace_fdm (natoms, icluster, a1vec, a2vec, a3vec)

   use grid
   use charges
   use interactions
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in) :: icluster
   integer, intent (in) :: natoms

   real, intent (inout), dimension (3)         :: a1vec
   real, intent (inout), dimension (3)         :: a2vec
   real, intent (inout), dimension (3)         :: a3vec

! Output
!   real, intent (inout) :: uHdcc 

! Local Parameters and Data Declaration
! ===========================================================================
   real, parameter :: fourpi = 1.0d0/(4.0d0*3.141592653589793238462643) 
! Local Variable Declaration and Description
! ===========================================================================

   integer i
   integer j
   integer jext
   integer ineigh
   integer i0, j0, k0
   integer iatom

   integer k,index

   real vpot
   real Etot_H

! Procedure
! ===========================================================================
   vcaG = 0.0d0
   Etot_H = 0.0d0
   do i = 0, (nrm-1)
     vpot = 0.0d0 
     do ineigh = 1,nneighij
       jext = n2e(i) + neighij(ineigh)
       j = e2n(jext)
!       write (*,200) i,ineigh,neighij(ineigh),n2e(i),jext,j
       vpot = vpot + drhoG(j)*drij(ineigh)
     enddo ! do ineigh
     vcaG(i) = vpot*fourpi
! Hartree energy int dVH*n_atm + int dVH*dn
     Etot_H =  Etot_H - rhoG0(i)*vcaG(i)*dvol - vcaG(i)*drhoG(i)*dvol
   enddo ! do i

! Get final value of Hartree double counting
!   uHdcc = 0.5d0*Etot_H
   write (*,100) Etot_H



! write out  hartree potential into dv.dat file
   open ( unit = 300, file = 'dpot.xsf', status = 'unknown' )
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
