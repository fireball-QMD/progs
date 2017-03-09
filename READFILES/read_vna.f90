! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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


! read_vna.f90
! Program Description
! ===========================================================================
!       This routine reads the neutral atomic poyential of each specie.
! The potential will be u sed to solve density problem in frame of numerical 
! grid method. 
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine read_vna ()

   use charges
   use grid
   use vnneutral
   use dimensions
   use constants_fireball
   use interactions
   use configuration 
   use integrals

   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
 
! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


   character (len = 100) filename
   character (len = 40) fileinvnn 

   integer issh
   integer ispec

   integer nzxvnn                   
   integer mesh
   integer inum
   integer iremainder
   integer ipoint

   real rcutoffvna
   real rcmax
   real xnoccwf
   real r
   real dr
   real sum

   real, dimension (:), allocatable :: psitemp

! test
   real, dimension (3):: vec
   real psi, dpsi
 
 
! Allocate Arrays
! ===========================================================================

 
! Procedure
! ===========================================================================
        
! root string

! loop over species
   do ispec = 1,nspecies

! append filename
      filename = trim(fdataLocation)//napot(0,ispec)

! get rcmax
      rcmax = 0.0d0
      do issh = 1,nssh(ispec) 
       if (rcutoff(ispec,issh) .gt. rcmax) rcmax = rcutoff(ispec,issh)
      enddo
      rcmax = rcmax/abohr


! open w.f. file
      write (*,200) filename 
      open (unit = 15, file = filename, status = 'unknown')

! read file header
      read (15,100) fileinvnn               ! filename (trash)
      read (15,*) nzxvnn                   ! atomic number (control)
      read (15,*) rcutoffvna               ! rcutoff
      read (15,*) mesh                     ! # mesh points
      read (15,*) sum                      ! trash

! Perform some checks
      if (mesh .gt. max_vna_points) then
         write (*,*) ' max_points_na = ', max_vna_points,' mesh = ', mesh
         write (*,*) ' The number of mesh points in your file is '
         write (*,*) ' greater than the dimensioned number of points.'
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readvnn'
      end if

      if (nzxvnn .ne. nzx(ispec)) then
         write (*,*) ' nzxvnn = ', nzxvnn, ' nzx = ', nzx(ispec)
         write (*,*) ' The cutoff radius in the wavefunction file, '
         write (*,*) ' for this shell, does not match the cutoff '
         write (*,*) ' radius that you put into the create.input file. '
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readvnn'
      end if

      if (abs(rcutoffvna - rcmax) .gt. 1.0d-5) then
         write (*,*) ' rcutoffvna = ', rcutoffvna,' rcutoff = ', rcmax
         write (*,*) ' The cutoff radius in the wavefunction file, for '
         write (*,*) ' this shell, does not match the cutoff radius '
         write (*,*) ' that you put into your create.input file. '
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readvnn'
      end if


! Read in the points
      mesh_na(ispec) = mesh
      do ipoint = 1, mesh
         read (15,300) rr_na(ipoint,ispec), vnna(ipoint,ispec)
      end do
      vnna(1,ispec) = vnna(2,ispec)
! set vna interval 
      drr_na(ispec) = rr_na(2,ispec) - rr_na(1,ispec)
! set max distance
      rmax_na(ispec) = rr_na(mesh,ispec)

! build spline       
      call buildspline2_1d (rr_na(1,ispec), vnna(1,ispec), mesh,         &
    &                       vnna_spline(1,ispec))

      close (unit = 15)

   enddo ! ispec

! write out vna potential into fort.* file for testing purposes
!   do ispec = 1, nspecies
!      mesh = 91 
!      rcmax = 0.0d0
!      do issh = 1,nssh(ispec) 
!       if (rcutoff(ispec,issh) .gt. rcmax) rcmax = rcutoff(ispec,issh)
!      enddo
!      rcmax = rcmax/abohr
!      dr = rcmax*abohr/dfloat(mesh - 1)
!      r = -dr
!      do ipoint = 1, mesh
!         r = r + dr
!         call getvna (ispec, r, psi, dpsi)
!         write (ispec*10,*) r, psi
!      end do
!   enddo

   write (*,*) '  '
   write (*,*) ' *---------------------  END READ_VNA  ----------------------*'

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format ( a25 )
200     format ( 4x,'Read neutral atomic potential from file: ', a40 )
300     format ( 2d24.16 )

  return
end subroutine read_vna
      
