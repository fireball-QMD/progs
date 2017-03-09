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


! ete_loop.f90
! Program Description
! ===========================================================================
!       This routine perform loop on electron time evolution
! This is variant where the electron wf are propagated using eigenstates
! so H,S are considered as entirely time-independent during MD time step
! at 1th step we have to obtain ground-state SCF solution H*psi = e*psi
! then psi is used as our adiabatic basis set where actual wfs are projected
! we also need get forces at 1-step to move ions at the end
! ===========================================================================
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ete_loop (itime_step)

        use options
        use tdse
        use outputs
        use scf

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
       integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================

       integer  ietime

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================

! be sure we update atomic neutral part of Hamiltonian at first step
        Kscf = 1

! Loop moving electrons
        do ietime = 1, netime

         write (*,*) '  '
         write (*,*) ' ****************************************************** '
         write (*,100) itime_step,ietime
         write (*,*) ' ****************************************************** '
! assemble H
         call assemble_h ()

! orthogonalize H
         call ortho_H (iqout, icluster)

! move psi(t) -> psi(t+dt)
         call propTpsi ()
         call tddenmat (iqout, icluster)
! write out psi
         if (iwrtpsit .eq. 1) then
! density matrix
          call tddenmat (iqout, icluster)
          call get_QMul (ifixcharge)
          call wrtout_psiT (itime_step,ietime)
         endif
! calculate rho, Qin
         if (iqout .eq. 2) then
! density matrix
          call tddenmat (iqout, icluster)
! Mulliken charge
          call get_QMul (ifixcharge)
         else
! Note: working in MO-space+McWeda we don't need to calculate rho
! Lowdin charge
          call get_QLow (ifixcharge)
         endif
! skip next time neutral atomic part of Hamiltonian
         Kscf = 2
        enddo ! ietime

        Kscf = 1
! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, ' MD time step =',i6,' electron time step =',i6)


        return
        end subroutine ete_loop

