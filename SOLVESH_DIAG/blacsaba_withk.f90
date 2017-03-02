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

