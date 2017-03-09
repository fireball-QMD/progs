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


! move_predic.f90
! Program Description
! ===========================================================================
!       This routine moves ions via MD according to forces
! JOM "predictor-like" part
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine move_predic (itime_step)

        use options
        use configuration
        use options
        use outputs
        use MD
        use neb
        use forces
        use fragments
        use constants_fireball
        use energy
        use barrier
        use optimization
        use interactions
        use charges


        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ix
        integer in1

        real vscale
        real ddxx

! Procedure
! ===========================================================================



           if(iensemble .eq. 0 .or. iensemble .eq. 1) then
! Now we do the prediction for the next step.  Fragment projection done inside
! this subroutine.
            call predictor (iquench, T_instantaneous, T_previous,      &
     &                    taurelax, T_want, etotold, etot, itime_step,    &
     &                    dt, ftot,.true.)
           end if
!---------------------------------------------------------------------------
! JOM : we should update vatom and ratom after writing out
! position/velocities (at least for velocity verlet, iensemble=3)
           if(iensemble .eq. 2 .or. iensemble .eq. 3) then
! use Martyna (Mol Phys 87, 1117-1157, 1996) with Nose Hoover (Chain)
! this is the first step
            if (iensemble .eq. 2) call NHCThermostat(dt,natoms,xmass,vatom)
            do iatom=1,natoms
             vatom(:,iatom)=vatom(:,iatom) + 0.5*fovermp*dt*ftot(:,iatom)/xmass(iatom)
            end do
! update particle positions
            ratom=ratom + dt*vatom
           end if
!---------------------------------------------------------------------------

! Scale velocities again after predictor step
           if(iensemble .eq. 1) then
            tkinetic = 0.0d0
            do iatom = 1, natoms
             tkinetic = tkinetic                                               &
     &        + (0.5d0/fovermp)*xmass(iatom)                                   &
     &        *(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
            end do
            T_instantaneous = (2.0d0/3.0d0)*tkinetic*kconvert/natoms
            vscale = sqrt(T_want/T_instantaneous)
            vatom(:,1:natoms) = vatom(:,1:natoms)*vscale
            xdot(1,:,1:natoms) = xdot(1,:,1:natoms)*vscale
            if (vscale .gt. 1.3d0 .or. vscale .lt. 0.7d0) then
             write (*,*) ' '
             write (*,*) ' Warning: significant velocity scaling needed after '
             write (*,*) ' predictor. Scaling factor = ', vscale
             write (*,*) ' '
            end if
           end if

! Shift atoms away from origin.
           if (ishiftO .eq. 1) then
            do iatom = 1, natoms
             ratom(:,iatom) = ratom(:,iatom) + shifter
             xdot(0,:,iatom) = xdot(0,:,iatom) + shifter
             if (ibarrier .eq. 1)                                             &
     &        ratom_final(:,iatom) = ratom_final(:,iatom) + shifter
            end do
           end if

! For the crude energy barrier calculation - exit if the final configuration
! is obtained.
           if (ibarrier .eq. 1 .and. barrier_achieved) stop

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
600     format (2x, ' average temperature       = ', f8.2)
601     format (2x, '                    nuclear kinetic energy = ', f18.8)
602     format (2x, ' Grand Total = nuclear kinetic + potential = ', f18.8)
603     format (2x, '                     total energy per atom = ', f18.8)
604     format (2x, '               grand total energy per atom = ', f18.8)
605     format (2x, '                        deltaE/atom  (meV) = ', f18.8)
700     format (2x, i2, 3(2x,f12.6))
701     format (2x, a2, 3(2x,f12.6))
1889    format (2x, ' deltae, etotnew, etotold = ', 3f15.8)
1888    format (2x, ' ACVE ', 3f15.8)
2000    format (3x, ' Checking MD convergence : ')
2001    format (2x,' ++++ iter = ',i8,' Etot= ',f16.8,' Fi_max= ',f14.6)
2010    format (2x,'  +++ Etot  RES =', f16.8,'  TOL = ',f16.8,'   CONVERGED ')
2011    format (2x,'  +++ Etot  RES =', f16.8,'  TOL = ',f16.8,'   NOT CONVERGED ')
2020    format (2x,'  +++ Fmax  RES =', f16.8,'  TOL = ',f16.8,'   CONVERGED ')
2021    format (2x,'  +++ Fmax  RES =', f16.8,'  TOL = ',f16.8,'   NOT CONVERGED ')
907     format (2x, ' ETOT = ', f15.6,'      T_instantaneous =  ', f12.4 )
908     format (2x, ' ETOT = ', f15.6,' (eV)  T =  ', f12.4,' (K)   Time = ', f20.2,' (fs)')
        return
        end subroutine move_predic

