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


! predictor.f90
! Program Description
! ===========================================================================
!       The predictor part of dynamics.
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine predictor (iquench, T_instantaneous, T_previous,  &
     &                        taurelax, T_want, etotold, etot, istep, &
     &                        dt, ftot,classOutput)
        use configuration
        use dimensions
        use fragments
        use constants_fireball
        use options, only: iclassicMD, verbosity
		use classicMD, only: freq_of_outputs
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
	logical, intent(in) :: classOutput
        integer, intent(in) :: iquench
        integer, intent(in) :: istep

        real, intent(in) :: dt
        real, intent(in) :: etot
        real, intent(in) :: taurelax
        real, intent(in) :: T_want

! Output
        real, intent (out) :: etotold
        real, intent(in), dimension (3, natoms) :: ftot
        real, intent(inout) :: T_instantaneous
        real, intent(inout) :: T_previous
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ifactor
        integer iorder
        integer isum
        integer ix
        integer kquench

        real dot
        real, external :: factorial
        real rfactor
        real xmasstot

        real, dimension (3) :: rcm 
        real, dimension (3) :: vcm 
        real, dimension (3) :: xlcm 

! Procedure
! ===========================================================================
! Calculate the center of mass position and velocity.
        rcm = 0.0d0
        vcm = 0.0d0
        xmasstot = sum(xmass) 
        do iatom = 1, natoms
         rcm(:) = rcm(:) + xmass(iatom)*(ratom(:,iatom)+ximage(:,iatom))
         vcm(:) = vcm(:) + xmass(iatom)*vatom(:,iatom)
        end do
        rcm = rcm/xmasstot
        vcm = vcm/xmasstot

! Calculate the center of mass angular momentum.
        xlcm = 0.0d0
        do iatom = 1, natoms
         xlcm(1) = xlcm(1) + ((ratom(2,iatom) + ximage(2,iatom))*xmass(iatom)*vatom(3,iatom) -   &
     &                        (ratom(3,iatom)+ximage(3,iatom))*xmass(iatom)*vatom(2,iatom))
         xlcm(2) = xlcm(2) + ((ratom(3,iatom)+ximage(3,iatom))*xmass(iatom)*vatom(1,iatom) -   &
     &                        (ratom(1,iatom)+ximage(1,iatom))*xmass(iatom)*vatom(3,iatom))
         xlcm(3) = xlcm(3) + ((ratom(1,iatom)+ximage(1,iatom))*xmass(iatom)*vatom(2,iatom) -   &
     &                        (ratom(2,iatom)+ximage(2,iatom))*xmass(iatom)*vatom(1,iatom))
        end do
!CHROM  - classic interaction
	if(verbosity .ge. 4)then
!END CHROM            
         write (*, 102) xlcm
         write (*, *) '  '
         write (*, 100) rcm
         write (*, 101) vcm
!CHROM
	endif
!END CHROM
 
! Move the atoms
! For the gear algorithm update both the positions and the velocities after 
! each prediction and correction
! Actual prediction step - usually 5th order Gear algorithm.
        do iatom = 1, natoms
         do iorder = 0, gear_order - 1
          do isum = iorder + 1, gear_order
           ifactor = isum - iorder
           xdot(iorder,:,iatom) = xdot(iorder,:,iatom)                       &
     &      + xdot(isum,:,iatom)*(dt**ifactor)/factorial(ifactor)
          end do
         end do
        end do

! Project out some forces if using FRAGMENTS
        if (numfrags .ne. 0) then
         call fixfrags (xmass)
        end if
        ratom(:,1:natoms) = xdot(0,:,1:natoms)
        vatom(:,1:natoms) = xdot(1,:,1:natoms)

        if (iquench .eq. 0) then
         write (*, 200) T_instantaneous, T_previous
         T_previous = T_instantaneous
         etotold = etot
         return
        end if
 
! Completely quench the velocities if necessary.  Quench the velocities on 
! every n`th step if iquench = +n or whenever the instantaneous temperature 
! (T_instantaneous) is lower than the instantaneous temperature on the previous ! step (T_previous) if iquench = -1.
        kquench = 0
        if (iquench .gt. 0 .and. mod(istep,iquench) .eq. 0) kquench = 1
        if (iquench .eq. -1 .and. T_instantaneous .lt. T_previous) kquench = 1
 
        if (kquench .eq. 1) then
!CHROM  - classic interaction
	if(classOutput)then
!END CHROM            
 	 write (*, 201)  T_instantaneous, T_previous
!CHROM  - classic interaction
	endif
!END CHROM            
         xdot(1,:,1:natoms) = 0.0d0
         vatom(:,1:natoms) = 0.0d0
         T_instantaneous = 0.0d0  ! quenched temperature
!CHROM  - classic interaction
        elseif(classOutput)then
!END CHROM            
         write (*, 202) T_instantaneous, T_previous
        end if

! To really "anneal" a cell we need to do something like constant temperature 
! molecular dynamics.  The simplest way to proceed is to rescale velocities at 
! each time step get some desired temperature.  If you set iquench = -2, then 
! this is accomplished by the  following snippet of code. The input parameters
! are the anneal temperature ("T_want") and taurelax, which specifies how 
! rapidly the cell is forced to T_want.  The variable taurelax is essentially 
! the relaxation time to make the cell recover T_want. See the code!
        if (iquench .eq. -2) then
         rfactor = sqrt((1.0d0 + dt/taurelax*(T_want/T_instantaneous - 1.0d0)))
         xdot(1,:,1:natoms) = xdot(1,:,1:natoms)*rfactor
         vatom(:,1:natoms) = xdot(1,:,1:natoms)
        endif
 
! Coordinate power quench
! Quench velocities one coordinate at a time.
        if (iquench .eq. -3) then
         do iatom = 1, natoms
          do ix = 1, 3
           dot = ftot(ix,iatom)*vatom(ix,iatom)
           if (dot .lt. 0.0d0) then
            vatom(ix,iatom) = 0.0d0
            xdot(1,ix,iatom) = 0.0d0
           end if
          end do
         end do
        end if
 
! Project out necessary forces if using FRAGMENTS
        if (numfrags .ne. 0) then
         call fixfrags (xmass)
         ratom(:,1:natoms) = xdot(0,:,1:natoms)
         vatom(:,1:natoms) = xdot(1,:,1:natoms)
        end if

        T_previous = T_instantaneous
        etotold = etot
 
! Format Statements
! ===========================================================================
100     format (2x, '         center of mass position = ', 3d12.4)
101     format (2x, '         center of mass velocity = ', 3d12.4)
102     format (2x, ' center of mass angular momentum = ', 3d12.4)
200     format (2x, ' Free Dynamics: T_instantaneous =  ', f12.4,      & 
      &         ' T_previous = ', f12.4)
201     format (2x, ' Quenching !! T_instantaneous =  ', f12.4,      & 
      &         ' T_previous = ', f12.4)
202     format (2x, ' No quenching. T_instantaneous =  ', f12.4,      & 
      &         ' T_previous = ', f12.4)

        return
        end
