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


! read_wf.f90
! Program Description
! ===========================================================================
!       This routine reads the radial wave functions of each specie.
! The radial wave functions are used to recontruct denisty on numerical grid.
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
 subroutine read_wf ()

   use charges
   use dimensions
   use constants_fireball
   use wavefunction
   use interactions
   use configuration
   use integrals
   use grid
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input

! Local Parameters and Data Declaration
! ===========================================================================
!   integer, parameter :: max_mesh = 5000 

! Local Variable Declaration and Description
! ===========================================================================


   character (len = 225) filename
   character (len = 40) fileinwf 

   integer issh
   integer ispec

   integer nzxwf                   
   integer mesh
   integer lqnwf      

   integer inum
   integer iremainder
   integer ipoint

   real rcutoffwf
   real rcmax
   real xnoccwf
   real r
   real sum
   real sum1
   real dr

   real, dimension (wfmax_points)  :: psitemp 
! test
   real psi
   real dpsi

 
 
! Allocate Arrays
! ===========================================================================

 
! Procedure
! ===========================================================================
! Read wave function
        
! loop over species
   do ispec = 1,nspecies


! loop over shell
      do issh = 1,nssh(ispec)

! open w.f. file
         filename = trim(fdataLocation)//wavefxn(issh,ispec)
         write (*,200) filename
         open (unit = 17, file = filename, status = 'unknown')

! read header
         read (17,100) fileinwf                 ! filename (trash)
         read (17,*) nzxwf                      ! atomic number (control)
         read (17,*) mesh                    ! # mesh points
         read (17,*) rcutoffwf, rcmax, xnoccwf ! rcutoff, occupancy
         read (17,*) lqnwf                      ! l-value

! Perform some checks
         if (nzxwf .ne. nzx(ispec)) then
            write (*,*) ' nzxwf = ', nzxwf, ' nzx = ', nzx
            write (*,*) ' The Z number the wavefunction file, for this '
            write (*,*) ' shell, does not match the cutoff radius '
            write (*,*) ' that you put into the create.input file. '
            write (*,*) ' Double check everything and rerun creator.'
            stop 'error in readpsi'
         end if

         if (rcutoffwf .le. (rcutoff(ispec,issh)/abohr - 1.0d-2) .or.     &
     &       rcutoffwf .ge. (rcutoff(ispec,issh)/abohr + 1.0d-2)) then
            write (*,*) ' rcutoffwf = ', rcutoffwf, ' rcutoff = ',        &
     &                                             rcutoff(ispec,issh)/abohr
            write (*,*) ' The cutoff radius in the wavefunction file, for '
            write (*,*) ' this shell, does not match the cutoff radius '
            write (*,*) ' that you put into your create.input file. '
            write (*,*) ' Double check everything and rerun creator.'
            stop 'error in readpsi'
         end if

         if (xnoccwf .ne. Qneutral(issh,ispec)) then
            write (*,*) ' xnoccwf = ', xnoccwf, ' xnoccin = ',             &
     &                                                Qneutral(issh,ispec)
            write (*,*) ' The occupation number in the wavefunction file, '
            write (*,*) ' for this shell, does not match the occupation '
            write (*,*) ' number that you put into your create.input file.'
            write (*,*) ' Double check everything and rerun creator.'
            stop 'error in readpsi'
         end if

         if (lqnwf .ne. lssh(issh,ispec)) then
            write (*,*) ' lqnwf = ', lqnwf, ' lqn = ', lssh(issh,ispec)
            write (*,*) ' The l quantum number in the wavefunction file, '
            write (*,*) ' for this shell, does not match the l quantum '
            write (*,*) ' number that you put into your create.input file.'
            write (*,*) ' Double check everything and rerun creator.'
            stop 'error in readpsi'
         end if

         if(mesh .gt. wfmax_points) then
            write (*,*) ' Error error ***** in read_wf. '
            write (*,*) ' Dimension of wavefunction = ', wfmax_points
            write (*,*) ' We are asking for mesh = ', mesh
            write (*,*) ' Redimension wfmax_points in MODULES/wavefunction.f90. '
            stop 'error in readpsi'
         end if

! Read in the points
         inum = idint(dfloat(mesh)/4)
         iremainder = mesh - (inum*4)
         do ipoint = 1, mesh - iremainder, 4
            read (17,300) psitemp(ipoint), psitemp(ipoint+1),             &
     &                     psitemp(ipoint+2), psitemp(ipoint+3)
         end do
! read final line of the data set
         if (iremainder .eq. 1) then
            read (17,300) psitemp(mesh)
         else if (iremainder .eq. 2) then
            read (17,300) psitemp(mesh-1), psitemp(mesh)
         else if (iremainder .eq. 3) then
            read (17,300) psitemp(mesh-2), psitemp(mesh-1), psitemp(mesh)
         end if
! close w.f. file
         close (17)

! Write psitemp to the wavefunction psi
         do ipoint = 1, mesh
            wf(ipoint,issh,ispec) = psitemp(ipoint)
         end do

! Set parameters 
         drr_wf(issh,ispec) = abohr*rcutoffwf/dfloat(mesh - 1)
         write (*,500) issh,ispec,drr_wf(issh,ispec)
         r = - drr_wf(issh,ispec)
         do ipoint = 1, mesh
            r = r + drr_wf(issh,ispec)
            rr_wf(ipoint,issh,ispec) = r
!            wf(ipoint,issh,ispec) = 2*exp(-2*r**2)
         end do
         mesh_wf (issh,ispec) = mesh
         rmax_wf (issh,ispec) = r

! Check normalization
         write (*,*) '      Checking normalization [NORM(l) should be 1]'

         sum = 0.0d0
         do ipoint = 1, mesh
            if (ipoint .ne. 1 .or. ipoint .ne. mesh) then
               sum = sum + drr_wf(issh,ispec)*                              &
     &             rr_wf(ipoint,issh,ispec)**2 * wf(ipoint,issh,ispec)**2
            else
               sum = sum + 0.5d0*drr_wf(issh,ispec)*                        &
     &             rr_wf(ipoint,issh,ispec)**2 * wf(ipoint,issh,ispec)**2
            end if
         end do ! do ipoint
         write (*,400) issh, sum
         write (*,*) '   --------------------------------------------------- ' 
! build spline 
         call buildspline2_1d (rr_wf(1,issh,ispec),wf(1,issh,ispec), mesh, &
                               wf_spline(1,issh,ispec))
! test orthogonalization with spline
         mesh = 40
         sum = 0.0d0
         write (*,*) '  Checking normalization for mesh: ',mesh
         dr = abohr*rcutoffwf/dfloat(mesh - 1)
         r = -dr
         do ipoint = 1, mesh
            r = r + dr
            call getpsi (ispec, issh, r, psi, dpsi)
            if (ipoint .ne. 1 .or. ipoint .ne. mesh) then
               sum = sum + dr*r**2 * psi**2
            else
               sum = sum + 0.5d0*dr*r**2 * psi**2        
            end if
         end do ! do ipoint
         write (*,400) issh, sum
         write (*,*) '   --------------------------------------------------- ' 
         
      enddo ! do nssh
   enddo ! ispec

   
! write out wavefunction into fort.* file for testing purposes
!   do ispec = 1, nspecies
!      do issh = 1,nssh(ispec)
!         write (*,*) 'ispec=',ispec,issh,'rmax = ',rmax_wf (issh,ispec)
!         mesh = 71 
!         dr = abohr*rcutoffwf/dfloat(mesh - 1)
!         r = -dr
!         do ipoint = 1, mesh
!            r = r + dr
!            call getpsi (ispec, issh, r, psi, dpsi)
!            write (ispec*100+issh, 350) r, psi
!         end do
!      enddo
!   enddo

   write (*,*) '  '
   write (*,*) ' *---------------------  END READ_WF  ----------------------*'

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format ( a25 )
200     format ( 4x,'Read wave function from file: ', a40 )
300     format ( 4d18.10 )       
350     format ( 2f12.6 )       
400     format ( 6x,' NORM (shell = ', i1, ') = ', f16.12 )
500     format ( 4x,' shell: ',i2,' ispec:',i2,' dr = ',f16.10)
  return
end subroutine read_wf
      
