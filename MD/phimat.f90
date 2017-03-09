! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
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


! phimat.f90
! Program Description
! ===========================================================================
!       Computes one column of the force constant matrix 
! phi(ix,iatom,jx,jatom) by displacing atom jatom in the jx direction by a 
! displacement u0disp. The row is stored in phirow(ix,iatom) and is written to 
! disk file filephi.
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ( old version from Kurt R. Glaesemann )
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine phimat( itime, natoms, ftot, iephc)

   use dimensions
   use dynamo
   use density
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
   integer,  intent (in) :: itime
   integer,  intent (in) :: natoms
   real,  dimension (3,natoms), intent (in) :: ftot
   integer,  intent (in) :: iephc

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
   real phirow(3, natoms_dm)
   integer iatom
   integer ix
   integer jx
   integer jatom
   integer jtime
   integer matom
   integer mx
   integer iprev
   integer iline
   integer ndim
   character (len=70) message
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================

   ndim = natoms_dm*ndx

!   do iatom = 1,natoms
!      write (100,301) iatom,(ftot(ix,iatom),ix=1,3)
!   enddo

! first (+) displacement, we only save forces vector
   if ( ldynamo ) then 

      ftot1  = ftot
      ldynamo = .false.

      return 

   endif

! Here we're only when second (-) displacement is done, it means 
! ldynamo must be '.false.'. We will evaluate a row of dynamical matrix.
 
! Set local time step

   if (mod(itime,2) .ne. 0) then 
      write (*,*) '   ---- Wrong time step for (-) displacement  ---', itime
      stop
   endif


! Determin which atom do we displace. See the ordering described in bvec.
   if (mod(ltime,ndx) .ne. 0) then
      matom = 1 + ltime / ndx
      mx = mod ( ltime, ndx )
   else
      matom = ltime / ndx
      mx = ndx
   end if

   
   if (ltime .gt. 1) then
! Determin how many lines belongs to one DM record.
      if (mod(ndim,4) .ne. 0) then
         iline = 1 + (ndim) / 4
      else
         iline = (ndim) / 4
      end if

! reopen the file 
      open ( unit=72, file=filephi, status='old' )

! read in previous results to get to end of file
! skip header of the record
      read (72,*) message
! set number of previous tracks
      iprev = ltime - 1
      do jtime = 1, iprev
! skip header of a individual track
         read (72,*) message
! skip data of the track
         do ix = 1, iline
            read (72,*) message
         enddo
      end do

! Write to a "fresh" file for istep=1
   else if (ltime .eq. 1) then
      open ( unit=72, file=filephi, status='unknown' )
! write header containing general information about dimension of dyn. mat.
      write (72,200) ndim, natoms_dm, ndx
   end if

   write(*,*) '  '
   write(*,400) itime, matom, mx

! Now construct phirow(ixx,iatom)
   do iatom = 1,natoms_dm
      ix = 1
      do jx = 1,3
         if (ndvec(jx) .eq. 1) then 
            jatom = jatoms_dm(iatom) 
            phirow(ix,iatom) = -(ftot1(jx,jatom) - ftot(jx,jatom)) / (2*u0disp)
            ix = ix + 1
         endif
      end do
   end do
! Now write to disk file filephi. We add to the end of the existing file.
   write (72,200) ltime, mx, matom
   write (72,300) ((phirow(ix,iatom), ix=1,ndx), iatom=1,natoms_dm)

! Assign the row to the global dynamical matrix
   jx = 1
   do iatom = 1,natoms_dm
      do ix = 1,ndx 
         phidm(ltime,jx) = phirow(ix,iatom)
         jx = jx + 1
      enddo
   enddo

! store derivatives of eigenvalues for e-ph couplings
   if (iephc .eq. 1) then
     do ix = 1, nephc
      deigen (ix,ltime) = eigsave(eiglist(ix),1)
      write(*,*) 'ooo',ix,eiglist(ix),ltime,deigen (ix,ltime)
     enddo 
   end if

! close the output file
   close (unit = 72)

! set flag 
   ldynamo = .true.
! increase local step
   ltime = ltime + 1

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format ('===========================================================')
101     format (' in phimat at time: ',i4,' doing ground state ..')
200     format ( 3i4)
300     format ( 4f16.8 )
!301     format ( i4,3f16.8 )
400     format ( ' in phimat at time: ',i4,' doing matom,mx=',i4,i2 )
 
   return
 end subroutine phimat
