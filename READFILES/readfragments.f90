! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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


! readfragments.f90
! Program Description
! ===========================================================================
!       Reads in the Fragments that are to be internally fixed
!
! Now first number in file FRAGMENTS is 0 or 1 and has following meaning 
!   0 .. position of selected atoms is frozen 
!   1 .. center of mass of selected atoms is frozen
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
! Modified by  P.Jelinek (2003)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
         subroutine readfragments ()
         use dimensions
         use configuration
         use fragments
         implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
         integer ifrag     ! loop over fragments
         integer imu       ! loop index
         integer iatom

         logical readsup   ! does the FRAGMENTS file exist

! Procedure
! ===========================================================================
! Initialize
         numfrags = 0
         fragtemp = 0
! dani.JOM
         nfragments = 0

! Check to see if we have anything to do
         inquire (file = 'FRAGMENTS', exist = readsup)
         if(.not. readsup) return

! Read in how many fragments
         allocate (fraggots(natoms))
         fraggots = 0
         open (unit = 71, file = 'FRAGMENTS', status = 'old')
         write (*,*) ' '
         write (*,*) ' Reading FRAGMENTS file'
         write (*,*) ' '

! Type of fragments 0 .. 1
         read (71,*) ifrags
         write(*,*) ' ifrags = ',ifrags
         if (ifrags .eq. 0) then 
          write (*,*) ' Listed atoms will be fixed in the initial positions.'
          fragtemp = 1
         else
          write (*,*) ' Block of selected atoms will be fixed.'
         end if

! Number of fragments
         read (71,*) numfrags
         if (numfrags .lt. 0) then
          numfrags = - numfrags
          fragtemp = 1
          write (*,*) ' The variable numfrags is negative, we will project '
          write (*,*) ' forces this means that temperatures will not include '
          write (*,*) ' inner fragment forces. '
          stop
         else
          write (*,*) ' The variable numfrags is positive, we will not project '
          write (*,*) ' forces this means that temperatures will include '
          write (*,*) ' inner fragment forces'
         end if
         if (numfrags .gt. natoms)then
          write (*,*) ' TOO MANY FRAGMENTS! numfrags = ', numfrags
          stop
         end if
         write (*,*) ' Number of individual fragments : ',numfrags

         allocate (fragsize (numfrags))
         allocate (fragatm (natoms,numfrags))
! define partially fixed atoms (0..free; 1..fixed)
         allocate (fragxyz (3,natoms))
         fragxyz = 0

! Loop over fragments
! first read in the size of the fragment and then the atoms in the fragment
         do ifrag = 1, numfrags
          write (*,*) ' Fragment # ', ifrag
          read (71,*) fragsize(ifrag)
! dani.JOM
          if (ifrags .eq. 0) then
           nfragments = nfragments + fragsize(ifrag)
          end if

          do imu = 1, fragsize(ifrag)
           read (71,*) iatom,fragxyz(1,iatom),fragxyz(2,iatom),fragxyz(3,iatom)
           fragatm(imu,ifrag) = iatom
           write (*,100) iatom,fragxyz(1,iatom),fragxyz(2,iatom),     &
      &      fragxyz(3,iatom)
           if (fraggots(iatom) .ne. 0) then
            write (*,*) ' atom # ', iatom, ' is in both fragment',           &
      &      fraggots(iatom), 'and ', ifrag
            stop
           end if
           fraggots(iatom) = ifrag
          end do
          write (*,*) ' '
         end do
         close (unit = 71)

! Store a copy of the original stucture (for later reference)
         allocate (ratom_frag(3,natoms))
         allocate (ratom_frag_save(3,natoms))
         ratom_frag_save (1:3,1:natoms) = ratom(1:3,1:natoms)

! Format Statements
! ===========================================================================a
100   format ('  atom: ',i6,3i6)
         return
         end

