! copyright info:
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


! readdata.f90
! Program Description
! ===========================================================================
!       This routine perform post-scf tasks like DOS
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine postscf ()

        use options
        use transport
        use outputs

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================

!        write (*,*) '  '
!        write (*,100)
!        write (*,*) ' Now we are performing post-scf analysis. '
!        write (*,*) '  '
!

! write out the dos if iwrtdos is greatter than 1 CGP
         if (iwrtdos .ge. 1) call writeout_dos ( )

! write out the hoppings if iwrthop is greatter than 1
         if (iwrthop .ge. 1) call writeout_hop

! write out the Atomo_i files if iwrtatom is greatter than 1
         if (iwrtatom .ge. 1) call writeout_atom ()

! write down non-orthogonal hamiltoninan and overlap for each atom
! and his neighbors
         if (iwrtatom .eq. 2) call hamtrans ()

! itrans
         if (itrans .eq. 1) call calcG ()

!jel-grid
         if (iwrtden .eq. 1) then
          write (*,*) ' Assemble atomic density '
!          call psi2mesh ()
          call assemble_KS_den0 ()
          write (*,*) ' Assemble differential density '
          call den2mesh (icluster)
          write (*,*) ' Call Poisson solver'
          call laplace_fft (icluster, ibias)
         endif


! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine postscf

