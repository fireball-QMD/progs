! copyright info:
!
! @Copyright 2005
! Fireball Enterprise Center, BYU
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

!                                 ==================== March 2010 ================
! changes - Zdenka Chromcova March 2010- because of classic MD:
! 1) outputs are not in each step but in intervals defined in Cdata/usePotential.in
! 2) some do-loops joined together (bacause of comp. time)
! 3) added vmax, if vmax*dt>2 then the program give warning that time step is too big
! 4) outputs
!                                 ================================================
! move_ions.f90
! Program Description
! ===========================================================================
! This routine moves ions via MD according to forces
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine move_ions (itime_step)

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
        use classicMD, only: freq_of_outputs

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom,k
        integer ix
        integer in1

        real vscale,vmax,rr(3)
        real ddxx
        logical :: writeOutput
   
! Procedures declared in this file
        external :: Move_ions_testTimeStep_cMD,Move_ions_writeOut_TempEtot
        external :: Move_ions_WriteOut_answerBasXyz, Move_ions_WriteOut_deltaEandF

! ===========================================================================
! JOM-test
! write(*,*)'JOM-TEST, move_ions CHANGED, dynamics NOT VALID!!!!!'
! ddxx = 0.0025
! ddxx = 0.010
! ratom(3,3) = ratom(3,3) - ddxx
! ratom(1,3) = ratom(1,3) - 0.014/0.8*ddxx
! vatom(3,3) = - ddxx/dt
! vatom(1,3) = - 0.014/0.8*ddxx/dt
! return
! JOM-test-end

! ============================================================================
! If we want internal fragment temperature to be zero, then we project
! out the forces, otherwise we wait until later (after T is calculated) to
! do the projections.

!in case of empirical potential the output is not written each cycle but with period freq_of_outputs
        if((iclassicMD == 0) .or. mod(itime_step,freq_of_outputs) == 0 )then
         writeOutput = .true.
        else
         writeOutput = .false.
        endif

        if (numfrags .ne. 0 .and. fragtemp .eq. 1)then
!         write(*,*) ' Projecting out inner-fragment forces. '
         call fixfrags2 (ftot)
        endif

! NVT isokinetic thermostat
        if (iensemble .eq. 1) call gaussT (natoms, vatom, xmass, T_want, ftot)

! 2nd step of velocity verlet algorithms
        if ((iensemble .eq. 2 .or. iensemble .eq. 3) .and. time .gt. 0.0) then
         vmax=0
         do iatom=1,natoms
          do k=1,3
           vatom(k,iatom)=vatom(k,iatom)+ 0.5*fovermp*dt*ftot(k,iatom)/xmass(iatom)
           if( abs(vatom(k,iatom)) > vmax ) vmax=abs(vatom(k,iatom))
          enddo
         enddo
         if (iensemble.eq.2) call NHCThermostat(dt,natoms,xmass,vatom)
         iwrtNHC=0 ! maybe include this in options.output someday
         if (iwrtNHC.eq.1) call writeHNose(time,dt,natoms,xmass,ratom,vatom,etot)
        endif

! Call umbrella sampling routine if iumbrella = 1
!        if (iumbrella .eq. 1)call get_umbrella (nstepi, nstepf, itime_step, time, natoms,&
!                                                 & iwrtfpieces, ratom, etot, ftot)

! Call steered routine if iumbrella = 2
!        if (iumbrella .eq. 2)call get_steered (nstepi, nstepf, itime_step,time,natoms,&
!                                                & iwrtfpieces, ratom, etot,ftot)


! Calculate crude energy barrier option.
! Push atoms along a direction "towards" the final configuation.
         if (ibarrier .eq. 1) then
         call push_atoms (natoms, ratom, etotold, etotnew, ftot, iquench,xmass)
         endif

! ----------------------------------------------------------------------------
! E Q U A T I O N S         O F         M O T I O N
! ----------------------------------------------------------------------------

        if(iensemble .eq. 0 .or. iensemble .eq. 1) then
! Now we have to do corrector, since we have all the forces. For the
! no-pressure run, call corrector (dynamics-corrector).
! NVT with velocity rescaling does this:
! 1. predict new position, velocity, etc. (in dynoPRED above)
! 2. calculate energy, density, forces, etc. at this new position
! 3. correct position, velocity, etc. (in corrector below)
! Repeat
! Better energy conservation can be obtained by additional corrections such
! as 1,2,3,2,3 and then repeat (or even more corrections).Simply
! correcting twice would double the computational work.It is fairly
! generally accepted to be a better idea to reduce the time step instead.
! Also, coding up multiple correction steps would be a pain.


!         if(writeOutput)then
!          write (*,*) ''
!          write (*,*) ' Predictor-Corrector: correct positions.'
!         endif

         call corrector (itime_step, dt, ftot)
         vmax=0
         do iatom = 1, natoms
          do k=1,3
           ratom(k,iatom) = xdot(0,k,iatom)
           vatom(k,iatom) = xdot(1,k,iatom)
           if( abs(vatom(k,iatom)) > vmax ) vmax=abs(vatom(k,iatom))
          enddo
         enddo
        endif

        if (numfrags .ne. 0 .and. fragtemp .eq. 1)then
         !write(*,*) ' Projecting out inner-fragment velocities. '
         call fixfrags ()
         vmax=0
         do iatom = 1, natoms
          do k=1,3
           ratom(k,iatom) = xdot(0,k,iatom)
           vatom(k,iatom) = xdot(1,k,iatom)
           if( abs(vatom(k,iatom)) > vmax ) vmax=abs(vatom(k,iatom))
          enddo
         enddo
        endif

        tkinetic = 0.0d0
        do iatom = 1, natoms
         tkinetic = tkinetic + (0.5d0/fovermp)*xmass(iatom) &
                & *(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
        enddo
! JOM-info : we need to correct this taking into account the fixed atoms
! (nfragments)
        T_instantaneous = (2.0d0/3.0d0)*tkinetic*kconvert/(natoms - nfragments)
!T_instantaneous = (2.0d0/3.0d0)*tkinetic*kconvert/natoms

! For constant temperature: rescale the velocities, so that the desired
! temperature is obtained. This makes sure that the work done to the forces
! in GaussT is not slightly off in velocities do to roundoff in dynoC
! The temperature we now have (3/2 kb * temperature = tkinetic)

         if ((iensemble .eq. 1 .and. .not. T_instantaneous .le. 0) ) then
         
         vscale = sqrt(T_want/T_instantaneous)
         call Move_ions_testTimeStep_cMD(vmax,dt,iclassicMD)

         if( verbosity .ge. 4) then
          write (*,*) ' Constant temp. (iensemble=1): scaling velocities. '
          write (*,*) ' Scaling =', vscale, ' It should be near 1.0!'
         endif
         do iatom = 1, natoms
          do k = 1,3
           vatom(k,iatom) = vatom(k,iatom)*vscale
           xdot(1,k,iatom) = xdot(1,k,iatom)*vscale
          enddo
         enddo
        endif
        T_average = ((itime_step - 1)*T_average + T_instantaneous)/itime_step
        getot = etot + tkinetic
        getotper = getot/natoms
! Check energy conservation
        deltaE = 1000.d0*(getotper - getot_initial)
        if (itime_step .eq. nstepi) getot_initial = getotper

! Write output on the screen
        call Move_ions_writeOut_TempEtot(writeOutput,T_instantaneous,T_average,tkinetic,getot,etotper,getotper,deltaE)

!---------------------------------------------------------------------------
! JOM : we should update vatom and ratom after writing out
! position/velocities (at least for velocity verlet, iensemble=3)
! if(iensemble .eq. 2 .or. iensemble .eq. 3) then
!! use Martyna (Mol Phys 87, 1117-1157, 1996) with Nose Hoover (Chain)
!! this is the first step
!if (iensemble .eq. 2) call NHCThermostat(dt,natoms,xmass,vatom)
!do iatom=1,natoms
! vatom(:,iatom)=vatom(:,iatom) + 0.5*fovermp*dt*ftot(:,iatom)/xmass(iatom)
!end do
!! update particle positions
!ratom=ratom + dt*vatom
! end if
!---------------------------------------------------------------------------

! ===========================================================================
! Write out positions/velocities
! ===========================================================================
! Project out FRAGMENTS before writing out
        if (numfrags .ne. 0 .and. fragtemp .eq. 0) then
         call fixfrags ()
         vmax = 0.d0
         do iatom = 1, natoms
          do k = 1,3
           ratom(k,iatom) = xdot(0,k,iatom)
           vatom(k,iatom) = xdot(1,k,iatom)
           if( abs(vatom(k,iatom)) > vmax ) vmax=abs(vatom(k,iatom))
          enddo
         enddo
        endif

! Possible center of mass constraint at each time step
        if (fixCenOfMass) then
         rcmNew = 0.0d0
         do iatom = 1, natoms
          rcmNew = rcmNew + xmass(iatom)*ratom(:,iatom)
         enddo
         rcmNew = rcmNew/xmasstot
         rcmDiff = rcmNew - rcmOld
         print *,'rcmDiff',rcmDiff
         rcmDiffMag=sqrt(rcmDiff(1)**2 + rcmDiff(2)**2 + rcmDiff(3)**2)
         if (rcmDiffMag .gt. 1.0d-6) then
          print *,'rcmDiffMag ',rcmDiffMag
          do iatom = 1, natoms
           ratom(:,iatom) = ratom(:,iatom) - rcmDiff
          enddo
         endif
        endif
! Shift atoms back towards origin.
        if (ishiftO .eq. 1) then
         do iatom = 1, natoms
          do k = 1,3
           ratom(k,iatom) = ratom(k,iatom) - shifter(k)
           xdot(0,k,iatom) = xdot(0,k,iatom) - shifter(k)
           if (ibarrier .eq. 1) ratom_final(k,iatom) = ratom_final(k,iatom) - shifter(k)
          enddo
         enddo
        endif
! write acceleration, coordinates, answer.bas, answer.xyz + outout on the screen

        if ( (MOD(itime_step,ntpr) .eq. 0 .and. itime_step .gt. 1 ) .or. (itime_step .eq. 1 ) .or. (itime_step .eq. nstepf)) then
           call Move_ions_WriteOut_answerBasXyz(writeOutput,time,itime_step,etot,T_instantaneous,dt)
        end if

! Update time
        time = time + dt

! control MD loop for quenching
        if (iquench .eq. -3 .or. iquench .eq. -2 .or. iquench .eq. -1 .or. &
                        & iquench .gt. 0) then
         deltaE = abs(etotnew - etotold)*natoms
! search max component of forces on free atoms
         deltaFmax = 0.0d0
         do iatom = 1,natoms
          do ix = 1,3
           deltaFmax = max(deltaFmax, abs(ftot(ix,iatom)*mask(ix,iatom)))
          enddo
         enddo

         call Move_ions_WriteOut_deltaEandF(writeOutput, itime_step, etotnew, etotold, &
                                        deltaFmax, deltaE, energy_tol, force_tol)

! check the criteria of convergence
         if ( (deltaE .lt. energy_tol) .and. (deltaFmax .lt. force_tol)) then
! dump info in a special format
          write (*,*) ' +++++Quenching process converged+++++ '
          write (*,*) 'That`sall for now, bye ..'
          stop
         endif

        endif

! Call TclMDUpdate routine to update 3D-visualization coordinates.
! if (itclmd .eq. 1) then
!if (mod(itime_step,itclmd_step) .eq. 0) then
! call tclmdupdate (natoms)
!end if
! end if

        if(iensemble .eq. 0 .or. iensemble .eq. 1) then
! Now we do the prediction for the next step.Fragment projection done inside
! this subroutine.
         call predictor (iquench, T_instantaneous, T_previous,&
                &taurelax, T_want, etotold, etot, itime_step,dt, ftot,writeOutput)
        endif
!---------------------------------------------------------------------------
! JOM : we should update vatom and ratom after writing out
! position/velocities (at least for velocity verlet, iensemble=3)
        if(iensemble .eq. 2 .or. iensemble .eq. 3) then
! use Martyna (Mol Phys 87, 1117-1157, 1996) with Nose Hoover (Chain)
! this is the first step
         if (iensemble .eq. 2) call NHCThermostat(dt,natoms,xmass,vatom)
          do iatom=1,natoms
           vatom(:,iatom)=vatom(:,iatom) + 0.5*fovermp*dt*ftot(:,iatom)/xmass(iatom)
          enddo
! update particle positions
          ratom = ratom + dt*vatom
          xdot(0,:,:) = ratom(:,:)
          xdot(1,:,:) = vatom(:,:) 
        endif
!---------------------------------------------------------------------------

 
! Scale velocities again after predictor step
        if(iensemble .eq. 1) then
         tkinetic = 0.0d0
         do iatom = 1, natoms
          tkinetic = tkinetic + (0.5d0/fovermp)*xmass(iatom) &
                &*(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
         enddo
         T_instantaneous = (2.0d0/3.0d0)*tkinetic*kconvert/natoms
         vscale = sqrt(T_want/T_instantaneous*(natoms-nfragments)/natoms)
!not take in acount the fix atoms ....  vscale = sqrt(T_want/T_instantaneous)
         do iatom=1,natoms
          do k = 1,3
           vatom(k,iatom) = vatom(k,iatom)*vscale
           xdot(1,k,iatom) = xdot(1,k,iatom)*vscale
          enddo
         enddo
         if (vscale .gt. 1.3d0 .or. vscale .lt. 0.7d0) then
          write (*,*) ' '
          write (*,*) ' Warning: significant velocity scaling needed after '
          write (*,*) ' predictor. Scaling factor = ', vscale
          write (*,*) ' '
         endif
        endif

! Shift atoms away from origin.
        if (ishiftO .eq. 1) then
         do iatom = 1, natoms
          do k = 1,3
           ratom(k,iatom) = ratom(k,iatom) + shifter(k)
           xdot(0,k,iatom) = xdot(0,k,iatom) + shifter(k)
           if (ibarrier .eq. 1)ratom_final(k,iatom) = ratom_final(k,iatom) + shifter(k)
          enddo
         enddo
        endif

         open (unit = 19, file = 'restart.xyz', status = 'unknown')
          write (19,*) natoms
          write (19,908) etot,T_instantaneous, itime_step*dt+init_time
          do iatom = 1, natoms
           in1 = imass(iatom)
            write (19,702) symbol(iatom), ratom(:,iatom) + ximage(:,iatom), vatom(:,iatom)
          enddo
         close (unit = 19)
! For the crude energy barrier calculation - exit if the final configuration
! is obtained.
        if (ibarrier .eq. 1 .and. barrier_achieved) stop

!         if(writeOutput)then
!         write(*,'(a,i7,a)') '=============================',itime_step,' ==================================='
!        endif


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
        100 format (2x, 70('='))
        702 format (2x, a2, 3(2x,f12.6), 3(2x,f16.9))
        908 format (2x, ' ETOT = ', f15.6,' eV; T =', f12.4,' K; Time = ', f12.1,' fs')

        return
end subroutine move_ions

subroutine Move_ions_testTimeStep_cMD(vmax,dt,iclassicMD)
        implicit none
        integer :: iclassicMD
        real :: vmax,dt

        if ( iclassicMD > 0 .and. vmax*dt > 2 ) then
         write(*,*)'=================================================================='
         write(*,*)'Warrning: The time step is probably too big: v*dt~',vmax*dt
         write(*,*)'Warrning: It should be better to set dt in fireball.in smaller.'
         write(*,*)'Warrning: acctual dt=',dt
         write(*,*)'=================================================================='
         stop
        endif
end subroutine

subroutine Move_ions_writeOut_TempEtot(writeOutput,T_instantaneous,T_average, &
                                        tkinetic,getot,etotper,getotper,deltaE)
        use options, only : verbosity
        implicit none
        logical :: writeOutput
        real :: T_instantaneous,T_average,tkinetic,getot,etotper,getotper,deltaE

        write (*,602) getot
        if( verbosity .ge. 1)then
         write (*,*)'T_instantaneous',T_instantaneous
          write (*,*) ''
         write(*,*) ' average temperature = ',T_average
         write (*,*) ''
          write (*,601) tkinetic
         write (*,603) etotper
          write (*,604) getotper
!                write (*,605) deltaE
         write (*,*)'deltaE/atom(meV) = ', deltaE
        endif

        600 format (2x, ' average temperature = ', f8.2)
        601 format (2x, 'nuclear kinetic energy = ', f18.8)
        602 format (2x, ' Grand Total = nuclear kinetic + potential = ', f18.8)
        603 format (2x, ' total energy per atom = ', f18.8)
        604 format (2x, ' grand total energy per atom = ', f18.8)
end subroutine

subroutine Move_ions_WriteOut_answerBasXyz(writeOutput,time,itime_step, &
                                            etot,T_instantaneous,dt)
        use outputs, only: iwrtvel,iwrtxyz
        use configuration, only: natoms,ratom,ximage,symbol,xyz2line,vatom,init_time
        use interactions, only: imass
        use charges, only: nzx
        use md, only: acfile,xvfile
        use options, only : restartxyz
        implicit none
        logical :: writeOutput
        integer :: itime_step,iatom,in1
        real :: etot,T_instantaneous,dt,time
!subroutines declaration
        external :: writeout_xv,writeout_ac

        if(writeOutput)then
         
!         call writeout_xv (xvfile, itime_step, time, imass, nzx, iwrtvel)
!         call writeout_ac (acfile, itime_step, time, imass, nzx)
! Open the answer file. This file has the input format and is overwritten
! every step. For empirical (potential) is written in freq_of_outputs cycle.
         open (unit = 17, file = 'answer.bas', status = 'unknown')
         if (iwrtxyz .eq. 1) then
          if (itime_step .eq. 1 .and. restartxyz .eq. 0) then
           open (unit = 18, file = 'answer.xyz', status = 'unknown')
          else
           open (unit = 18, file = 'answer.xyz', status = 'unknown',&
                                                 & position = 'append')
          endif
          endif

         write (17,*) natoms
         do iatom = 1, natoms
          in1 = imass(iatom)
          write (17,700) nzx(in1), ratom(:,iatom) + ximage(:,iatom)
         enddo
          close (unit = 17)

         if (iwrtxyz .eq. 1) then
          write (18,*) natoms
!d@ni
          if(xyz2line.eq.0) write (18,*) ''
          if(xyz2line.eq.1) write (18,907) etot,T_instantaneous
          if(xyz2line.eq.2) write (18,908) etot,T_instantaneous,itime_step*dt+init_time
!          if(xyz2line.eq.2 .and. restartxyz .eq. 1) write (18,908) etot,T_instantaneous,(itime_step-1)*dt+init_time
          if(xyz2line.eq.3) write (18,*) etot
          do iatom = 1, natoms
           in1 = imass(iatom)
            write (18,701) symbol(iatom), ratom(:,iatom) + ximage(:,iatom)
          enddo
          close (unit = 18)
         endif
        endif

        700 format (2x, i2, 3(2x,f12.6))
        701 format (2x, a2, 3(2x,f12.6))
        702 format (2x, a2, 3(2x,f12.6), 3(2x,f16.9))
        907 format (2x, ' ETOT = ', f15.6,'      T_instantaneous =', f12.4 )
        908 format (2x, ' ETOT = ', f15.6,' eV; T =', f12.4,' K; Time = ', f12.1,' fs')
end subroutine

subroutine Move_ions_WriteOut_deltaEandF(writeOutput,itime_step, &
                        etotnew, etotold, deltaFmax, deltaE, energy_tol, force_tol)
        implicit none
        logical :: writeOutput
        real :: etotnew, etotold, deltaFmax, deltaE, energy_tol, force_tol
        integer :: itime_step

         if(writeOutput)then
         write (*,1889) deltaE, etotnew, etotold
         write (*,*) 'deltaFmax = ', deltaFmax,'force_tol =', force_tol
! dump information on the geometric convergence
         write (*,2000)
         write (*,2001) itime_step, etotnew, deltaFmax
         if (deltaE .lt. energy_tol) then
          write (*,2010) deltaE, energy_tol
         else
          write (*,2011) deltaE, energy_tol
         endif
         if (deltaFmax .lt. force_tol) then
          write (*,2020) deltaFmax, force_tol
         else
          write (*,2021) deltaFmax, force_tol
         endif
        endif

        1889 format (2x, ' deltae, etotnew, etotold = ', 3f15.8)
        1888 format (2x, ' ACVE ', 3f15.8)
        2000 format (3x, ' Checking MD convergence : ')
        2001 format (2x,' ++++ iter = ',i8,' Etot= ',f16.8,' Fi_max= ',f14.6)
        2010 format (2x,'+++ EtotRES =', f16.8,'TOL = ',f16.8,' CONVERGED ')
        2011 format (2x,'+++ EtotRES =', f16.8,'TOL = ',f16.8,' NOT CONVERGED ')
        2020 format (2x,'+++ FmaxRES =', f16.8,'TOL = ',f16.8,' CONVERGED ')
        2021 format (2x,'+++ FmaxRES =', f16.8,'TOL = ',f16.8,' NOT CONVERGED ')
end subroutine
