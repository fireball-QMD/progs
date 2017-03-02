      SUBROUTINE PCLAPUTTER( A, DESCA, ALPHA, NORBITALS)
      implicit none
!
!
!     Written by Kurt R. Glaesemann
!
!     EXTENSION OF, BUT NOT PART OF
!  -- ScaLAPACK tools routine (version 1.5) --
!
      real*8, intent (out), dimension (*) :: a
      integer, intent (in), dimension (*) :: desca
      integer, intent (in) :: norbitals
      real*8, intent (in), dimension (norbitals,norbitals) :: alpha

!     .. Things that used to be Scalar Arguments but we fix ..
      INTEGER            IA, ICPRNT, IRPRNT, JA, M, N
!     ..
!
!  Purpose
!  =======
!
!  This stores values of local array alpha (*) in A on
!  the the process of coordinates (IRPRNT, ICPRNT)
!
!  Notes
!  =====
!
!  Each global data object is described by an associated description
!  vector.  This vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  Arguments
!  =========
!
!  M       (global input) INTEGER
!          The number of rows to be operated on i.e the number of rows
!          of the distributed submatrix sub( A ). M >= 0.
!
!  N       (global input) INTEGER
!          The number of columns to be operated on i.e the number of
!          columns of the distributed submatrix sub( A ). N >= 0.
!
!  A       (local input) COMPLEX pointer into the local memory to a
!          local array of dimension (LLD_A, LOCc(JA+N-1) ) containing
!          the local pieces of the distributed matrix sub( A ).
!
!  IA      (global input) INTEGER
!          The row index in the global array A indicating the first
!          row of sub( A ).
!
!  JA      (global input) INTEGER
!          The column index in the global array A indicating the
!          first column of sub( A ).
!
!  DESCA   (global and local input) INTEGER array
!          The array descriptor for the distributed matrix A.
!
!  IRPRNT  (global input) INTEGER
!          The row index of the printing process.
!
!  ICPRNT  (global input) INTEGER
!          The column index of the printing process.
!
!  Norbitals (local input) INTEGER size of alpha matrix
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            CTXT_, MB_, NB_, LLD_
      PARAMETER          ( CTXT_ = 2, MB_ = 5, NB_ = 6, LLD_ = 9 )
!     ..
!     .. Local Scalars ..
      INTEGER            KK, I, IACOL, IAROW, IB, ICTXT, ICURCOL,  &
     &                   ICURROW, II, IIA, IN, J, JB, JJ, JJA, JN, K,  &
     &                   LDA, NB_A, MB_A, MYCOL, MYROW, NPCOL, NPROW
!     ..
!     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
!     ..
!     .. Executable Statements ..
!
!     These are usually inputs, but we will require square matrices
!     and that A begins at (0,0), and that the master get the data
!     This simplyfies the subroutine call.
!
      IA=1
      JA=1
      IRPRNT=0
      ICPRNT=0
      N=NORBITALS
      M=NORBITALS
!
!     Get grid parameters
!
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
!
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,  &
     &              IIA, JJA, IAROW, IACOL )
      ICURROW = IAROW
      ICURCOL = IACOL
      II = IIA
      JJ = JJA
      LDA = DESCA( LLD_ )
      NB_A = DESCA( NB_ )
      MB_A = DESCA( MB_ )
!
!     Handle the first block of column separately
!
      JN = MIN( ICEIL( JA, NB_A ) * NB_A, JA+N-1 )
      IN = MIN( ICEIL( IA, MB_A ) * MB_A, IA+M-1 )
      JB = JN-JA+1
      IB = IN-IA+1
      J = JA
      I = IA
      IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
         IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
            DO KK = 0, JB-1
               DO K = 0, IB-1
                  A(II+K+(JJ+KK-1)*LDA)=alpha (I+K, J+KK)
               END DO
            END DO
         END IF
      ELSE
         IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
            CALL DGERV2D( ICTXT, IB, JB, A( II+(JJ-1)*LDA ), LDA,  &
     &                    IRPRNT, ICPRNT )
         ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
            CALL DGESD2D( ICTXT, IB, JB, alpha (I, J), norbitals,  &
     &                    ICURROW, ICURCOL )
         END IF
      END IF
      IF( MYROW.EQ.ICURROW ) II = II + IB
      ICURROW = MOD( ICURROW+1, NPROW )
      CALL BLACS_BARRIER( ICTXT, 'All' )
!
!     Loop over remaining block of rows (in first column)
!
      DO 50 I = IN+1, IA+M-1, MB_A
         IB = MIN( MB_A, IA+M-I )
         IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
            IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               DO KK = 0, JB-1
                  DO K = 0, IB-1
                     A( II+K+(JJ+KK-1)*LDA )=alpha(I+K, J+KK)
                  END DO
               END DO
            END IF
         ELSE
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL DGERV2D( ICTXT, IB, JB, A( II+(JJ-1)*LDA ),  &
     &                       LDA, IRPRNT, ICPRNT )
            ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               CALL DGESD2D( ICTXT, IB, JB, alpha(I, J), norbitals,  &
     &                       ICURROW, ICURCOL )
            END IF
         END IF
         IF( MYROW.EQ.ICURROW ) II = II + IB
         ICURROW = MOD( ICURROW+1, NPROW )
         CALL BLACS_BARRIER( ICTXT, 'All' )
   50 CONTINUE
!
      II = IIA
      ICURROW = IAROW
!
      IF( MYCOL.EQ.ICURCOL ) JJ = JJ + JB
      ICURCOL = MOD( ICURCOL+1, NPCOL )
      CALL BLACS_BARRIER( ICTXT, 'All' )
!
!     Loop over remaining column blocks
!
      DO 130 J = JN+1, JA+N-1, NB_A
         JB = MIN(  NB_A, JA+N-J )
         IN = MIN( ICEIL( IA, MB_A ) * MB_A, IA+M-1 )
         IB = IN-IA+1
         I = IA
         IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
            IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               DO KK = 0, JB-1
                  DO K = 0, IB-1
                     A( II+K+(JJ+KK-1)*LDA )=alpha(I+K, J+KK)
                  END DO
               END DO
            END IF
         ELSE
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL DGERV2D( ICTXT, IB, JB, A( II+(JJ-1)*LDA ),  &
     &                       LDA, IRPRNT, ICPRNT )
            ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               CALL DGESD2D( ICTXT, IB, JB, alpha(I, J), norbitals,  &
     &                       ICURROW, ICURCOL )
            END IF
         END IF
         IF( MYROW.EQ.ICURROW ) II = II + IB
         ICURROW = MOD( ICURROW+1, NPROW )
         CALL BLACS_BARRIER( ICTXT, 'All' )
!
!        Loop over remaining block of rows
!
         DO 110 I = IN+1, IA+M-1, MB_A
            IB = MIN( MB_A, IA+M-I )
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO KK = 0, JB-1
                     DO K = 0, IB-1
                        A( II+K+(JJ+KK-1)*LDA )=alpha(I+K, J+KK)
                     END DO
                  END DO
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL DGERV2D( ICTXT, IB, JB, A( II+(JJ-1)*LDA ),  &
     &                          LDA, IRPRNT, ICPRNT )
               ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  CALL DGESD2D( ICTXT, IB, JB, alpha(I, J),  &
     &                          norbitals, ICURROW, ICURCOL )
               END IF
            END IF
            IF( MYROW.EQ.ICURROW ) II = II + IB
            ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'All' )
  110    CONTINUE
!
         II = IIA
         ICURROW = IAROW
!
         IF( MYCOL.EQ.ICURCOL ) JJ = JJ + JB
         ICURCOL = MOD( ICURCOL+1, NPCOL )
         CALL BLACS_BARRIER( ICTXT, 'All' )
!
  130 CONTINUE
!
      RETURN
      END
