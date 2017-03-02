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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! initconstraints.f90
! Program Description
! ===========================================================================
!       This routines sets up the basis vector, ratom(3,natoms).  Any 
! constraints that are desired are then used to re-evaluate and initialize
! starting  velocites. 
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initconstraints (iconstraints, iensemble,         &
     &                              T_initial, ibarrier, ratom_final, imass, &
     &                              fixCenOfMass, rcmOld, xmassTot)
        use dimensions
        use configuration
        use constants_fireball
        use fragments
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ibarrier
        integer, intent (in) :: iensemble

        integer, intent (in), dimension (4) :: iconstraints
        integer, intent (in), dimension (natoms) :: imass

        real, intent (in) :: T_initial

        logical, intent(in) :: fixCenOfMass

! Output
        real, intent (inout), dimension (3, natoms) :: ratom_final
        real, intent (out), dimension (3) :: rcmOld
        real, intent (out) :: xmassTot
 
! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ix
        integer num_atoms

        real T_instantaneous
        real tkinetic
        real vscale
        real vscale1, vscale2
        real xmasstot2
 
        real, dimension (3) :: rcm
        real, dimension (3) :: vcm

        logical velocityfile
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize the velocities to zero
        vatom = 0.0d0

! Initialize velocites with a random uniform or Gaussian distribution.
!         call random_seed
        if ((iensemble .eq. 0) .or. (iensemble .eq. 3)) then  ! NVE ensembles or quenching
!----------------------------------------------------------------------------
! JOM-info : initialize with a normal distribution corresponding to
! T_initial
         do iatom = 1, natoms
          do ix = 1, 3
           call random_number (vscale)
! JOM-info random uniform
!          vatom(ix,iatom) =                                                 &
!    &      sqrt(fovermp*T_initial/(kconvert*xmass(iatom)))*(2.0d0*vscale - 1.0d0)
!----------------------------------------------------------------------------
! JOM-info random normal
           vscale1 = sqrt(-2.0d0*log(vscale))
           vatom(ix,iatom) =                                                 &
     &      sqrt(fovermp*T_initial/(kconvert*xmass(iatom)))*vscale1
           call random_number (vscale)
           vscale2 = cos(2.0d0*pi*vscale)
           vatom(ix,iatom) = vatom(ix,iatom)*vscale2                  
!----------------------------------------------------------------------------
          end do
         end do
!        do iatom = 1, natoms
!         do ix = 1, 3
!          call random_number (vscale)
!          vatom(ix,iatom) =                                                 &
!    &      sqrt(fovermp*600.0d0/(kconvert*xmass(iatom)))*(2.0d0*vscale - 1.0d0)
!         end do
!        end do
!----------------------------------------------------------------------------
        else   ! NVT ensembles--gaussian
         do iatom = 1, natoms
          do ix = 1, 3
           call random_number (vscale)
           vscale = 1.0d0 - vscale
           vatom(ix,iatom) =                                                 &
     &      7.4181d-04*sqrt(T_initial*abs(log(vscale))/xmass(iatom)) 
           call random_number (vscale)
           if (vscale .lt. 0.5d0) vatom(ix,iatom) = - vatom(ix,iatom)
          end do
         end do
        end if

! Read velocities from a velocities file. Note: if this is done, then it will
! wipe out the velocities originally initialized from a random temperature
! distribution. 
        inquire (file = 'VELOCITIES', exist = velocityfile)
        if (velocityfile) then
         if (iensemble .gt. 0 .or. T_initial .gt. 0.0d0) then 
          write (*,*) ' WARNING! Wiping out velcities that were initialized ' 
          write (*,*) ' from a constant temperature option or T_initial was '
          write (*,*) ' set greater than zero! '
         end if
         write (*,*) '  '
         write (*,*) ' We are reading from a velocity file. '
         open (unit = 12, file = 'VELOCITIES', status = 'old')
         read (12,*) num_atoms
         if (num_atoms .ne. natoms) then
          write (*,*) ' The velocity file that you are using must not belong '
          write (*,*) ' to the basis file that you are now calculating. '
          write (*,*) ' The number of atoms differs between the two. '
         end if
         do iatom = 1, natoms
          read (12,*) vatom(:,iatom) 
         end do
         close (unit = 12)
        end if
         
! Calculate the center-of-mass position.
        rcm = 0.0d0
        xmasstot2 = sum(xmass)
        do iatom = 1, natoms
         rcm = rcm + xmass(iatom)*ratom(:,iatom)
        end do
        rcm = rcm/xmasstot2
        write (*,*) '  '
        write (*,100) rcm
        write (*,*) '  '
        rcmOld = rcm
        xmassTot = xmasstot2

! constraint #1
! ----------------
! Shift new ratom so they are measured from the center of mass.
! in this way, rcmmol = 0.
        if (iconstraints(1) .eq. 1) then
         write (*,*) ' Constraining the positions about the center-of-mass '
         do iatom = 1, natoms
          ratom(:,iatom) = ratom(:,iatom) - rcm
         end do
        end if

! constraint #4
! ----------------        
        do iatom = 1,natoms
        end do
        if (iconstraints(4) .eq. 1) then
         write (*,*) ' Constraining angular momentum!!!! '
         call zero_ang_mom ()
        end if

! constraint #2
! -----------------
! Now adjust the velocities to get velocity of vcm = 0
        if (iconstraints(2) .eq. 1) then
         write (*,*) ' Constraining the velocities about the center-of-mass '
         vcm = 0.0d0
         do iatom = 1, natoms 
          vcm = vcm + xmass(iatom)*vatom(:,iatom)
         end do
         vcm = vcm/xmasstot
         do iatom = 1, natoms
          vatom(:,iatom) = vatom(:,iatom) - vcm
         end do
        end if

! constraint # 3
! ----------------
! Finally rescale the velocities to get the average temp = temperature_want
! tkinetic = average kinetic energy per particle in ev.
        if (iconstraints(3) .eq. 1) then
         write (*,*) ' Rescaling the velocities based on desired temperature. '
         tkinetic = 0.0d0
         do iatom = 1, natoms
          tkinetic = tkinetic + 0.5d0*xmass(iatom)/fovermp                   &
     &     *(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
         end do
! JOM-info : if "nfragments" atoms are fixed
! JOM-info : if only x,y or x coordiantes are fixed, we have to change
! this
!        tkinetic = tkinetic/natoms
         tkinetic = tkinetic/(natoms - nfragments)

! The temperature we now have (3/2 kb * T_instantaneous = tkinetic )
         T_instantaneous = (2.0d0/3.0d0)*tkinetic*kconvert
!--------------------------------------------------------------------------
! JOM-test
         write(*,*)'T_instantaneous',T_instantaneous
         write(*,*)'T_initial',T_initial
! JOM-test
!--------------------------------------------------------------------------
         if (T_instantaneous .gt. 0.0d0) then
          vscale = sqrt(T_initial/T_instantaneous)
         else
          vscale = 0.0d0
         end if
         vatom(:,1:natoms) = vatom(:,1:natoms)*vscale
        end if

! ***************************************************************************
! Now write out the velocities.
        write (*,*) '  '
        write (*,*) ' Atom Velocities: '
        write (*,200)
        write (*,201)
        write (*,200)
        do iatom = 1, natoms
         write (*,202) iatom, symbol(iatom), vatom(:,iatom), imass(iatom)
        end do

! ***************************************************************************
!
!       Shift the coordinates
! ***************************************************************************
! Initialize the shifting constants_fireball. This is to make sure that none of the
! atoms fall on (0.0, 0.0, 0.0), because the three-center integrals does
! not like it when this occurs.
! Need to avoid getting r1vector-cross-r2vector = 0 in nanlxc routine.
        shifter(1) = 4.0d0*atan(1.0d0)    ! pi
        shifter(2) = 1.0/exp(1.0d0)       ! 1/e
        shifter(3) = sqrt(2.0d0)          ! square root of 2

! Now shift coordinates to form computerese basis vectors.
        if (ishiftO .eq. 1) then
         do iatom = 1, natoms
          ratom(:,iatom) = ratom(:,iatom) + shifter
          if (ibarrier .eq. 1)                                               &
     &     ratom_final(:,iatom) = ratom_final(:,iatom) + shifter
         end do
         if (fixCenOfMass) rcmOld = rcmOld + shifter 
        end if
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, ' Calculated position of the Center-of-Mass: ', 3(2x,f7.3))
200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 6x,   &
     &              ' x ', 9x, ' y ', 9x, ' z ', 6x, ' Species # ')
202     format (3x, i5, 7x, a2, 3(2x,f10.5), 7x, i2)
 
        return
        end
