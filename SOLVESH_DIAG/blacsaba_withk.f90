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


! blacsaba.f90
! Program Description
! ===========================================================================
!       This subroutine does a(ia,ja)=alpha*a(ia,ja) for a BLACS distributed
! matrix.
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
      SUBROUTINE blacsaba( A, IA, JA, DESCA, ALPHA,  &
     &                     MYCOL, MYROW, NPCOL, NPROW )

! Argument Declaration and Description
! ===========================================================================
! Input
      implicit none
      INTEGER, intent (in) :: IA, JA, MYCOL, MYROW, NPCOL, NPROW
      INTEGER, intent (in) :: DESCA( * )
      REAL   , intent (in) :: ALPHA

! Input/Output
      COMPLEX, intent (inout) :: A( * )

! Local Parameters and Data Declaration
! ===========================================================================
      INTEGER             LLD_
      PARAMETER          ( LLD_ = 9 )

! Local Variable Declaration and Description
! ===========================================================================
      INTEGER       IACOL, IAROW, IIA, JJA

! Procedure
! ===========================================================================

!     Get grid parameters.

      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,  &
     &              IAROW, IACOL )
! Do it!
      IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
           A( IIA+(JJA-1)*DESCA( LLD_ ) ) =  &
     &     A( IIA+(JJA-1)*DESCA( LLD_ ) ) * ALPHA
      END IF
      RETURN
      END

