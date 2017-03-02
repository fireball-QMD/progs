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

! push_atoms.f90
! Program Description
! ===========================================================================
!       This routine pushes the atoms from a given initial configuration to
! a final configuration.  It is used for calculating a crude energy 
! barrier between two such configurations. 
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine push_atoms (natoms, ratom, etotold, etotnew, ftot, iquench, &
     &                         xmass)
        use dimensions
        use barrier
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms

        real, intent (in) :: etotnew
        real, intent (in) :: etotold

        real, intent (in), dimension (3, natoms) :: ratom
        real, intent (in), dimension (natoms) :: xmass
 
! Output
        real, intent (inout), dimension (3, natoms) :: ftot

        integer, intent (out) :: iquench

! Local Parameters and Data Declaration
! ===========================================================================
! use_mass determines whether or not you mass weight things
! Since we are talking about forces and not velocities, probably a bad idea
        logical, parameter :: use_mass = .false.
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom

        real distance
        real dotp
        real rms
        real mass

        real stn, dstn

        real, dimension (3) :: rhat
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! First step - calculate unit vectors between current and final configurations.
! In addition, calculate the rms difference between the current and final
! configurations.  If the rms difference is less than some tolerance (as
! determined in the barrier.optional file), then end the simulation.
        rms = 0.0d0
        do iatom = 1, natoms
         mass = 1.0d0
         if (use_mass) mass = xmass(iatom)
         rhat(:) = ratom_final(:,iatom) - ratom(:,iatom)
         distance = sqrt(rhat(1)**2 + rhat(2)**2 + rhat(3)**2)

! skip a case when actual and final position of atom are practically identical
         if (distance**2 .gt. 0.00001d0) then          
          rms = rms + distance**2
          rhat = rhat/distance

! Next step - check ftot*rhat.  If ftot*rhat is negative, then "push" the
! atom by strength barrier_push.  If ftot*rhat is positive, then do nothing,
! since this atom is going to the correct final configuration on its own accord
          dotp = ftot(1,iatom)*rhat(1) + ftot(2,iatom)*rhat(2)                &
     &         + ftot(3,iatom)*rhat(3)

! Wrong direction correct this (bar_sav is how much of the bad force to keep)
          if (dotp .lt. 0.0d0) then  
           ftot (:,iatom)=ftot(:,iatom)*bar_sav +                              &
     &          (barrier_push*mass - dotp*bar_sav)*rhat(:)
           dotp = ftot(1,iatom)*rhat(1) + ftot(2,iatom)*rhat(2)                &
     &         + ftot(3,iatom)*rhat(3)
          endif

! Going in the right direction, but with too much force
          if (dotp/mass .gt. bar_too_much) then
           ftot (:,iatom) = ftot (:,iatom)*bar_too_much/(dotp/mass)
           dotp = ftot(1,iatom)*rhat(1) + ftot(2,iatom)*rhat(2)                &
     &         + ftot(3,iatom)*rhat(3)
          endif

! Not completly wrong, but still bad (we smoothly go from pushing to no pushing)
! Requires the length of the force vector
          distance = sqrt(ftot(1,iatom)**2 + ftot(2,iatom)**2 + ftot(3,iatom)**2)
          if (dotp/distance .lt. bar_stop_push) then
           call smoother (dotp/distance, bar_stop_push, 0.8, stn, dstn)
           ftot (:,iatom) = ftot (:,iatom)*(1.0d0-stn) +                      &
     &            (ftot(:,iatom)*bar_sav + barrier_push*mass*rhat(:))*stn
          endif
         else
! the actual and the final position are identical, we reset forces on the atom
           ftot (:,iatom) = 0.0d0
         end if ! if (distance**2)

        end do ! do iatom

        rms = sqrt(rms/(3*natoms))

! Stop the simulation if the current and final configurations are the same.
        barrier_achieved = .false.
        if (rms .lt. barrier_tol) barrier_achieved = .true.
        write (*,*) '  '
        write (*,*) ' The current and final configurations are within each '
        write (*,*) ' other by the the following rms deviation = ', rms
        write (*,*) '  '

! Once we have gotten over the barrier, then go back to dynamically quenching
! in order to stop the simulation from oscillating wildly. 
        iquench = bar_how_often
        if (etotold .gt. etotnew) iquench = -1
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
