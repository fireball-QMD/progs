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


! readbarrier.f90
! Program Description
! ===========================================================================
!       This routine reads in the barrier.optional file.  This is an option
! that exists for calculating a crude energy barrier.  

!     The additions described here are designed to "push" the system of 
! interest from a given initial configuration to a given final configuration 
! via a "path of least resistance".  Here's how it's supposed to work:

!    At each time step we take the vector difference, for each atom, between
! the current position (at time step #1 this is the initial configuration)
! and the final, or desired, configuration.  These difference vectors are
! then transformed into unit vectors.  These unit vectors then describe the 
! direction that each atom must move in order to reach the final configuration.

!    Next we examine the calculated forces for each atom.  We ask the question,
! "Does the force on iatom have a nonnegative component in the desired 
! direction?" (Is force dot unit direction .ge. 0?).  If the answer is yes, 
! then we are satisfied that iatom is not being pushed away from its desired 
! final position and we do nothing.  But, if the answer is no, we add enough 
! force IN THE DESIRED DIRECTION to make the component of the force in the 
! desired direction equal to a given value.

!    This procedure insures that iatom always has some component of
! acceleration towards the desired final postion of iatom.  Given enough time
! iatom will arrive at the desired final position.  (Two rare cases are 
! ignored here, the force on iatom may be exactly zero, and the dot product 
! may be exactly zero, but these two cases could be handled if necessary.)
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
        subroutine readbarrier (natoms, nspecies, imass, nzx) 
        use dimensions
        use barrier
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies

        integer, intent (in), dimension (natoms) :: imass
        integer, intent (in), dimension (nspec_max) :: nzx

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ispec
        integer natoms_in
        integer nucz

        integer, dimension (natoms) :: imass_in

        character (len=30) file_barrier

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize all of the positions of the final structure to zero.
        ratom_final = 0.0d0

        if (ibarrier .eq. 1) then 
         write (*,*) '  '
         write (*,100) 
         write (*,*) ' You have chosen to calculate a crude energy barrier '
         write (*,*) ' between an initial configuration and a final '
         write (*,*) ' configuration. '

         open (unit = 79, file = 'barrier.optional', status = 'old')

         write (*,*) '  '
         write (*,*) ' An amount of force is added to the calculated forces '
         write (*,*) ' to push the atoms from the initial configuration to ' 
         write (*,*) ' the final configuration. '
         write (*,*) ' Insert the strength of the "push": '
         read (79,*) barrier_push
         write (*,*) ' barrier_push = ', barrier_push
        
         write (*,*) '  ' 
         write (*,*) ' The simulation continues until the rms deviation ' 
         write (*,*) ' between the current and final configurations is within '
         write (*,*) ' some accepted tolerance. ' 
         write (*,*) ' Insert the tolerance desired: '
         read (79,*) barrier_tol
         write (*,*) ' barrier_tol = ', barrier_tol
        
         write (*,*) '  '
         write (*,*) ' There is an option to quench velocities every X '
         write (*,*) ' steps.  You might need this to remove off excess heat '
         read (79,*) bar_how_often
         write (*,*) ' bar_how_often = ', bar_how_often

         write (*,*) '  '
         write (*,*) ' There is an option to push atoms that are technically'
         write (*,*) ' heading in the right direction, but only barely. '
         write (*,*) ' Acceptable values are 0(no push) to 1(push all) '
         read (79,*) bar_stop_push
         write (*,*) ' bar_stop_push = ', bar_stop_push

         write (*,*) '  '
         write (*,*) ' There is an option to slow down atoms that are '
         write (*,*) ' heading in the right direction, but very rapidly. '
         write (*,*) ' This is to avoid overshooting the target. '
         write (*,*) ' Acceptable values are 0(set forces to zero-dumb) to '
         write (*,*) ' a huge number (32000, never quench)  '
         read (79,*) bar_too_much
         write (*,*) ' bar_too_much = ', bar_too_much 

         write (*,*) '  '
         write (*,*) ' There is an option to save some of the original force '
         write (*,*) ' of the atoms that are heading in the wrong direction. '
         write (*,*) ' Acceptable values are 0(just push, save nothing) to '
         write (*,*) ' 1(save all and push) '
         read (79,*) bar_sav
         write (*,*) ' bar_sav = ', bar_sav

         write (*,*) '  '
         write (*,*) ' Insert the filename with the final configuration: '
         read (79,200) file_barrier
         write (*,201) file_barrier

         open (unit = 80, file = file_barrier, status = 'old')
         read (80, *) natoms_in
         if (natoms_in .ne. natoms) then
          write (*,*) ' Sorry natoms_in = ', natoms_in, ' and natoms = ', natoms
          write (*,*) ' The initial and final configuration files are not the '
          write (*,*) ' same structures - please check! '
          stop
         end if

! Loop over the number of atoms
         do iatom = 1, natoms
          read (80,*) nucz, ratom_final(:,iatom)
          do ispec = 1, nspecies
           if (nucz .eq. nzx(ispec)) then
            imass_in(iatom) = ispec
           end if
          end do
          if (imass_in(iatom) .ne. imass(iatom)) then
           write (*,*) ' iatom = ', iatom
           write (*,*) ' Sorry imass_in(iatom) = ', imass_in(iatom)
           write (*,*) ' and imass(iatom) = ', imass(iatom)
           write (*,*) ' The initial and final configuration files are not the '
           write (*,*) ' same structures - please check! '
           write (*,*) ' The same ordering must be used in the two files. '
           stop
          end if
         end do
        end if
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (a30)
201     format (2x, a30)
 
        return
        end
