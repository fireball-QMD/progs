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


! move_correc.f90
! Program Description
! ===========================================================================
!       This routine moves ions via MD according to forces
! JOM : "corrector-like" part
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine move_correc (itime_step)

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
! JOM-test
!       write(*,*)'JOM-TEST, move_ions CHANGED, dynamics NOT VALID!!!!!'
!       ddxx = 0.0025
!       ddxx = 0.010
!       ratom(3,3) = ratom(3,3) - ddxx 
!       ratom(1,3) = ratom(1,3) - 0.014/0.8*ddxx 
!       vatom(3,3) = - ddxx/dt
!       vatom(1,3) = - 0.014/0.8*ddxx/dt
!       return
! JOM-test-end


!#ifdef CCOMPILE
!-----------------------------------------------------------------------------
!thermodynamic integration server-client force interchange
!-----------------------------------------------------------------------------
!          if (ithermoint .eq. 1) then
!           open(unit=400,file='forceValues.txt',status='unknown')
!           do imarker=1,natoms
!            write(400,*)(ftot(isub,imarker),isub=1,3)
!           end do
!           close(400)
!           call cclient()
!           open(unit=401,file='newforce.txt',status='unknown')
!           write(*,*)'new forces are'
!           do imarker=1,natoms
!            read(401,*)(ftot(isub,imarker),isub=1,3)
!            write(*,*)(ftot(isub,imarker),isub=1,3)
!           end do
!           close(401)
!          endif
!----------------------------------------------------------------------------
!#endif

! ============================================================================
! If we want internal fragment temperature to be zero, then we project
! out the forces, otherwise we wait until later (after T is calculated) to
! do the projections.
           if (numfrags .ne. 0 .and. fragtemp .eq. 1)then
            !write(*,*) ' Projecting out inner-fragment forces. '
            call fixfrags2 (ftot)
           end if

! NVT isokinetic thermostat
           if (iensemble .eq. 1)           &
     &      call gaussT (natoms, vatom, xmass, T_want, ftot)

! 2nd step of velocity verlet algorithms
           if ((iensemble .eq. 2 .or. iensemble .eq. 3) .and. time .gt. 0.0) then
            do iatom=1,natoms
             vatom(:,iatom)=vatom(:,iatom)                                  &
              + 0.5*fovermp*dt*ftot(:,iatom)/xmass(iatom)
            end do
            if (iensemble.eq.2) call NHCThermostat(dt,natoms,xmass,vatom)
            iwrtNHC=0 ! maybe include this in options.output someday
            if (iwrtNHC.eq.1)                                                     &
     &       call writeHNose(time,dt,natoms,xmass,ratom,vatom,etot)
           end if

! Call umbrella sampling routine if iumbrella = 1
           if (iumbrella .eq. 1)                                              &
     &      call get_umbrella (nstepi, nstepf, itime_step, time, natoms,      &
     &                        iwrtfpieces, ratom, etot, ftot)


! Call steered routine if iumbrella = 2
           if (iumbrella .eq. 2)                                              &

     &      call get_steered (nstepi, nstepf, itime_step, time, natoms,       &
     &                        iwrtfpieces, ratom, etot, ftot)




! Calculate crude energy barrier option.
! Push atoms along a direction "towards" the final configuation.
           if (ibarrier .eq. 1) then
            call push_atoms (natoms, ratom, etotold, etotnew, ftot, iquench,  &
     &                       xmass)
           end if

! ----------------------------------------------------------------------------
!                   E Q U A T I O N S    O F    M O T I O N
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
! as 1,2,3,2,3 and then repeat (or even more corrections).  Simply
! correcting twice would double the computational work.  It is fairly
! generally accepted to be a better idea to reduce the time step instead.
! Also, coding up multiple correction steps would be a pain.
            !write (*,*) '  '
            !write (*,*) ' Predictor-Corrector: correct positions. '
            call corrector (itime_step, dt, ftot)
            ratom(:,1:natoms) = xdot(0,:,1:natoms)
            vatom(:,1:natoms) = xdot(1,:,1:natoms)
           end if

           if (numfrags .ne. 0 .and. fragtemp .eq. 1)then
            !write(*,*) ' Projecting out inner-fragment velocities. '
            call fixfrags ()
            ratom(:,1:natoms) = xdot(0,:,1:natoms)
            vatom(:,1:natoms) = xdot(1,:,1:natoms)
           end if

           tkinetic = 0.0d0
           do iatom = 1, natoms
            tkinetic = tkinetic                                               &
     &       + (0.5d0/fovermp)*xmass(iatom)                                   &
     &       *(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
           end do
! JOM-info : we need to correct this taking into account the fixed atoms
! (nfragments)
           T_instantaneous = (2.0d0/3.0d0)*tkinetic*kconvert/(natoms - nfragments)                   
!          T_instantaneous = (2.0d0/3.0d0)*tkinetic*kconvert/natoms
          ! write(*,*)'T_instantaneous',T_instantaneous

! For constant temperature: rescale the velocities, so that the desired
! temperature is obtained.   This makes sure that the work done to the forces
! in GaussT is not slightly off in velocities do to roundoff in dynoC
! The temperature we now have (3/2 kb * temperature = tkinetic)

           if (iensemble .eq. 1 .and. .not. T_instantaneous .le. 0) then
            vscale = sqrt(T_want/T_instantaneous)
!            write (*,*) ' Constant temp. (iensemble=1): scaling velocities. '
!            write (*,*) ' Scaling =', vscale, ' It should be near 1.0!'
            vatom(:,1:natoms) = vatom(:,1:natoms)*vscale
            xdot(1,:,1:natoms) = xdot(1,:,1:natoms)*vscale
           end if
           T_average = ((itime_step - 1)*T_average + T_instantaneous)/itime_step
          ! write (*,*) '  '
          ! write (*,600) T_average
          ! write (*,*) '  '
          ! write (*,601) tkinetic
           getot = etot + tkinetic
          ! write (*,602) getot
           getotper = getot/natoms
          ! write (*,603) etotper
          ! write (*,604) getotper
           if (itime_step .eq. nstepi) getot_initial = getotper
! Check energy conservation
           deltaE = 1000.d0*(getotper - getot_initial)
          ! write (*,605) deltaE
!---------------------------------------------------------------------------
! JOM : we should update vatom and ratom after writing out
! position/velocities (at least for velocity verlet, iensemble=3)
!           if(iensemble .eq. 2 .or. iensemble .eq. 3) then
!! use Martyna (Mol Phys 87, 1117-1157, 1996) with Nose Hoover (Chain)
!! this is the first step
!            if (iensemble .eq. 2) call NHCThermostat(dt,natoms,xmass,vatom)
!            do iatom=1,natoms
!             vatom(:,iatom)=vatom(:,iatom) + 0.5*fovermp*dt*ftot(:,iatom)/xmass(iatom)
!            end do
!! update particle positions
!            ratom=ratom + dt*vatom
!           end if
!---------------------------------------------------------------------------

! ===========================================================================
! Write out positions/velocities
! ===========================================================================
! Project out FRAGMENTS before writing out
           if (numfrags .ne. 0 .and. fragtemp .eq. 0) then
            call fixfrags ()
            ratom(:,1:natoms) = xdot(0,:,1:natoms)
            vatom(:,1:natoms) = xdot(1,:,1:natoms)
           end if

! Possible center of mass constraint at each time step
           if (fixCenOfMass) then
            rcmNew = 0.0d0
            do iatom = 1, natoms
             rcmNew = rcmNew + xmass(iatom)*ratom(:,iatom)
            end do
            rcmNew = rcmNew/xmasstot
            rcmDiff = rcmNew - rcmOld
            print *,'rcmDiff',rcmDiff
            rcmDiffMag=sqrt(rcmDiff(1)**2 + rcmDiff(2)**2 + rcmDiff(3)**2)
            if (rcmDiffMag .gt. 1.0d-6) then
             print *,'rcmDiffMag ',rcmDiffMag
             do iatom = 1, natoms
              ratom(:,iatom) = ratom(:,iatom) - rcmDiff
             end do
            end if
           end if

! Shift atoms back towards origin.
           if (ishiftO .eq. 1) then
            do iatom = 1, natoms
             ratom(:,iatom) = ratom(:,iatom) - shifter
             xdot(0,:,iatom) = xdot(0,:,iatom) - shifter
             if (ibarrier .eq. 1)                                             &
     &        ratom_final(:,iatom) = ratom_final(:,iatom) - shifter
            end do
           end if

           call writeout_xv (xvfile, itime_step, time, imass, nzx, iwrtvel)
           call writeout_ac (acfile, itime_step, time, imass, nzx)
!vlada begin
	do iatom =1, natoms
           write (211,101) iatom,vatom(1,iatom),vatom(2,iatom),vatom(3,iatom),sqrt((vatom(1,iatom)**2)+(vatom(2,iatom)**2)+(vatom(3,iatom)**2))
        end do
! vlada end
! Update time
           time = time + dt

! Open the answer file. This file has the input format and is overwritten
! every step.
           open (unit = 17, file = 'answer.bas', status = 'unknown')
           if (iwrtxyz .eq. 1) then
            if (itime_step .eq. 1) then
             open (unit = 18, file = 'answer.xyz', status = 'unknown')
            else
             open (unit = 18, file = 'answer.xyz', status = 'unknown',        &
     &             position = 'append')
            end if
           end if
           write (17,*) natoms
           do iatom = 1, natoms
            in1 = imass(iatom)
            write (17,700) nzx(in1), ratom(:,iatom) + ximage(:,iatom)
           end do
           close (unit = 17)
           if (iwrtxyz .eq. 1) then
            write (18,*) natoms
!d@ni       
            if(xyz2line.eq.0) write (18,*) '  '
            if(xyz2line.eq.1) write (18,907) etot,T_instantaneous
            if(xyz2line.eq.2) write (18,908) etot,T_instantaneous, itime_step*dt
            if(xyz2line.eq.3) write (18,*) etot
            do iatom = 1, natoms
             in1 = imass(iatom)
             write (18,701) symbol(iatom), ratom(:,iatom) + ximage(:,iatom)
            end do
            close (unit = 18)
           end if

! control MD loop for quenching
           if (iquench .eq. -3 .or. iquench .eq. -2 .or. iquench .eq. -1 .or. &
     &         iquench .gt. 0) then
            deltaE = abs(etotnew - etotold)*natoms
! search max component of forces on free atoms
            deltaFmax = 0.0d0
            do iatom = 1,natoms
             do ix = 1,3
              deltaFmax = max(deltaFmax, abs(ftot(ix,iatom)*mask(ix,iatom)))
             end do
            end do
            write (*,1889) deltaE, etotnew, etotold
            write (*,*) '  deltaFmax = ', deltaFmax,'  force_tol =', force_tol

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
! check the criteria of convergence
            if ( (deltaE .lt. energy_tol) .and. (deltaFmax .lt. force_tol)) then
! dump info in a special format
             write (*,*) '   +++++  Quenching  process converged  +++++ '
             write (*,*) '            That`s  all for now, bye ..'
! old format
!             write (*,1888) Vouc, etot, (etot - atomic_energy)/natoms
!             write (*,*) ' deltaFmax = ',deltaFmax,' force_tol =',force_tol
! write out volume and final energy (for e.o.s. and b.m.)
!             open (unit = 215, file = 'eos.intermediate', status = 'unknown')
!             write (215,*) Vouc
!             write (215,*) etotnew
!             close (unit=215)
             stop
            endif
!            if (itime_step .eq. nstepf .and. ((deltaE .gt. energy_tol) .or.    &
!     &         (deltaFmax .gt. force_tol))) then
!             write(*,*) ' Failed, ', ' Vouc =',Vouc,' etot =', etot
!             write(*,*) ' deltaFmax = ', deltaFmax,' force_tol =', force_tol
!            endif
           end if

! Call TclMDUpdate routine to update 3D-visualization coordinates.
!         if (itclmd .eq. 1) then
!          if (mod(itime_step,itclmd_step) .eq. 0) then
!           call tclmdupdate (natoms)
!          end if
!         end if


! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format ( 1i4, 4f14.8)
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
        end subroutine move_correc

