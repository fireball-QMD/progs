! iopyright info:
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


! initdenmat.f90
! Program Description
! ===========================================================================
!    The subrotine initializes density matrix corresponding to atomic density
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! Program Declaration
! ===========================================================================
 subroutine initdenmat ( natoms )

   use density
   use neighbor_map
   use interactions
   use charges
   use options
   use grid
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent(in)    :: natoms

! Output

! Local Parameters and Data Declaration
! ===========================================================================
!   real, parameter :: pi = 3.141592653589793238462643

! Local Variable Declaration and Description
! ===========================================================================

   integer imu
   integer inu
   integer i
   integer issh
   integer iatom
   integer jatom
   integer mbeta
   integer matom
   integer in1
   integer in2
   integer ineigh
   integer a1
   integer a2
   integer a3

   real qmu
   logical isfile

! Allocate arrays
! ===========================================================================

! Procedure
! ===========================================================================
   write (*,*) '   Initialize density matrix (using atomic configuration) '
! reset variables
   rho = 0.0d0
   rhoA = 0.0d0
   vnaG = 0.0d0
   vcaG = 0.0d0
   vxcG = 0.0d0
   drhoG = 0.d00

! Loop over all atoms iatom in the unit cell
   do iatom = 1, natoms
      in1 = imass(iatom)
! find atom intself in list of neighbors
      matom = -99
      do ineigh = 1, neighn(iatom)
         mbeta = neigh_b(ineigh,iatom)
         jatom = neigh_j(ineigh,iatom)
         if (iatom .eq. jatom .and. mbeta .eq. 0) matom = ineigh
      end do
      imu = 1
      do issh = 1, nssh(in1)
! evaluate neutral charge per orbital
         qmu = Qneutral(issh,in1) / real(2*lssh(issh,in1)+1)
         do i = 1, (2*lssh(issh,in1)+1)
! set diagonal part of density matrix
            rhoA(imu,iatom) = qmu
            rho(imu,imu,matom,iatom) = qmu
            imu = imu + 1
         enddo ! do i
      enddo ! do imu
   end do ! do iatom

! if denmat.dat file exists
   inquire (file = 'denmat.dat', exist = isfile)
   if (isfile .and. iks .eq. 1) then
     write (*,*) ' Read denmat.dat '
     rho = 0.0d0
     rho_old = 0.0d0
     open (unit = 60, file = 'denmat.dat', status = 'old')
! loop over atoms
     do iatom = 1, natoms
       in1 = imass(iatom)
       read (60,*) a1,a2,a3
       if ( (neighn(iatom) .ne. a2) .and. (num_orb(in1) .ne. a3) ) then
        write (*,*) 'Input data are inconsistent !!'
        write (*,*) ' iatom =',iatom
        write (*,*) ' neighn =',a2, neighn(iatom)
        write (*,*) ' num_orb =',a3, num_orb(in1)
        write (*,*) ''
        stop
       endif
! loop over neighbors
       do ineigh = 1, neighn(iatom)
         jatom = neigh_j(ineigh,iatom)
         in2 = imass(jatom)
         read (60,*) a1,a2,a3
         if ( (jatom .ne. a2) .and. (num_orb(in2) .ne. a3) ) then
           write (*,*) 'Input data are inconsistent !!'
           write (*,*) ' ineigh =',ineigh
           write (*,*) ' jatom =',a2, jatom
           write (*,*) ' num_orb =',a3, num_orb(in1)
           write (*,*) ''
           stop
         endif
         do imu = 1, num_orb(in1)
           do inu = 1, num_orb(in2)
            read (60,*) rho(imu,inu,ineigh,iatom)
           end do ! do inu
         end do ! do imu
       end do ! do ineigh
     end do ! do iatom
! close file
     close (60)

!     write (*,*) ' Map psi onto the mesh'
!     call psi2mesh ()
!     write (*,*) ' Map psi*psi onto the mesh'
!     call psi22mesh ()
! assemble atomic density on the grid
     write (*,*) ' Assemble atomic density '
     call assemble_KS_den0 ()
! assemble density on the grid
     write (*,*) ' Assemble density'
     call assemble_KS_den ()

   endif

! copy rho to rho_old
   rho_old = rho

! Format Statements
! ===========================================================================
100     format (2x, i5, 2x, a40, 2x, i2)
200     format (2x, 70('='))

   return
 end subroutine initdenmat
