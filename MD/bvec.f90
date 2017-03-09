! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University of Utah - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Motorola, Physical Sciences Research Labs - Alex Demkov
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio University - Dave Drabold

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


! bvec.f90
! Program Description
! ===========================================================================
! This changes b ---> b, and b---> b+/-u0disp (with b=atom positions)
! depending on the index itime. The +/- is + if iplumin=1, and - if iplumin=-1.
! In a given time step, we call it once at the beginning to add u0disp to a
! single atom, and then at the end of the time step, we call it to subtract
! off the u0disp to get back to the undisplaced situation. We do this so
! that at the end of the loop, all the atoms are in undisplaced positions
! because thats the way we wish to start each time step. It makes it really
! easy to start the program at some intermediate place using this method.
! We displace one atom at a time and in one direction.
! For each itime, we displace the atom given by index iatom, and in the ix
! direction.
! The sequence for which atom and for which direction gets the displacement
! is as follows:
!
!     itime       ix       iatom
!       1          1         1
!       2          2         1
!    nspace      nspace      1
!    nspace+1      1         2
!     .            .         .
!    2*nspace    nspace      2
!    2*nspace+1    1         3
!      .           .         .
!      .           .         .
!      .           .         .
!  natoms*nspace nspace   natoms
!
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
!
! (rewritten by P. Jelinek February 2005)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
  subroutine bvec (itime, natoms, nstepf, ratom)

    use dimensions
    use dynamo

    implicit none

! Argument Declaration and Description
! ===========================================================================
    integer, intent (in)  :: itime
    integer, intent (in)  :: natoms
    integer, intent (in)  :: nstepf
    real, dimension (3,natoms), intent (inout) :: ratom

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
    integer jatom
    integer iatom
    integer jx
    real    dir

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================

    write (*,*) '  '
    write (*,*) '  Welcome to bvec subroutine.  '


    if ( ldynamo ) then
       dir = 1.0d0
    else
       dir = -1.0d0
    endif
! ========================================================================
! First we restore the initial atomic position.
    do iatom = 1, natoms
       ratom(:,iatom) = ratom0(:,iatom)
    enddo

! ========================================================================
! Second, we displace i-th atom of the list for the next step.
! We skip this part at last step.
    if (itime .lt. nstepf ) then
! Now we find l-th atom and a related direction in which we move the atom.
       if (mod(ltime,ndx) .ne. 0) then
          iatom = 1 + ltime / ndx
          jx = mod (ltime, ndx)
       else
          iatom = ltime / ndx
          jx = ndx
       end if

! find index of the atom signed to displace
       jatom = jatoms_dm(iatom)
       write (*,*) '  '
       write (*,200) jatom,jx

       if (jatom .gt. natoms) then
          write (*,*) 'error: jatom to big!!'
          stop
       end if

! we now displace atom jatom in jx direction
       ratom(jx,jatom) = ratom(jx,jatom) + dir*u0disp


    endif



! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100    format ('    Reset atom: ',i4,' in direction ix: ',i2)
200    format ('    Displace atom: ',i4,' in direction ix: ',i2)
300    format (' ====================================================== ')
!301    format ( i4,3f14.6)
    return
  end subroutine bvec
