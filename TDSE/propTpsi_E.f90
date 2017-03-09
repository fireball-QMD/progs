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


! propT_psi.f90
! Program Description
! ===========================================================================
!       This routine subroutine propagates wave-function in time base on
! a simple formulae of exponetial phase factor
!  |i(t+dt)> = exp(-i/h*Ek*dt)*|A>*<A|i(t)>
!
! ===========================================================================
! Code rewritten by:
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
        subroutine propTpsi ()

        use tdse
        use interactions
        use options
        use kpoints
        use density

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
!       integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================
       real, parameter :: hbar = 0.65822d0


! Local Variable Declaration and Description
! ===========================================================================
       integer  ielec
       integer  imu
       integer  inu
       integer  ikpoint
       integer  iband

       real fact

       complex phase
       complex cphi



! Procedure
! ===========================================================================
! printout
        write (*,*) '        ---- Propagate Psi in time ----'
! reset wf coefficients
        psi = 0.0d0
        if (icluster .eq. 1) then
! K-loop
         do ikpoint = 1, nkpoints
! loop over bands
          do iband = 1, norbitals
! calculate exp factor (exp(-i/h*ek*dt)
           fact = eigen_k(iband,ikpoint)*dte/hbar
           phase = cmplx(cos(fact),-1.0d0*sin(fact))
! loop over electrons
           do ielec = 1, nelec
! loop over orbitals
            do imu = 1, norbitals
             cphi = cmplx(bbnkre(imu,iband,ikpoint),0.0d0)
             psi(imu,ielec,ikpoint) = psi(imu,ielec,ikpoint) +                  &
      &           phase*cphi*psi2phi(iband,ielec,ikpoint)
            end do ! imu
           end do ! ielec
          end do ! iband
         end do ! ikpoint

        else
! K-loop
         do ikpoint = 1, nkpoints
! loop over bands
          do iband = 1, norbitals
! calculate exp factor (exp(-i/h*ek*dt)
           fact = eigen_k(iband,ikpoint)*dte/hbar
           phase = cmplx(cos(fact),-1.0d0*sin(fact))
! loop over electrons
           do ielec = 1, nelec
! loop over orbitals
            do imu = 1, norbitals
             cphi = cmplx(bbnkre(imu,iband,ikpoint),bbnkim(imu,iband,ikpoint))
             psi(imu,ielec,ikpoint) = psi(imu,ielec,ikpoint) +                  &
      &           phase*cphi*psi2phi(iband,ielec,ikpoint)
            end do ! imu
           end do ! ielec
          end do ! iband
         end do ! ikpoint
        endif ! if(icluster)

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (8f12.6)


        return
        end subroutine propTpsi

