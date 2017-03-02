MODULE MatMult

INTERFACE OPERATOR(.x.)
 MODULE PROCEDURE MatTimesMat,  MatTimesVector,  VectorTimesVector, &
                  CMatTimesCMat,CMatTimesCVector,CVectorTimesCVector
END INTERFACE

CONTAINS

FUNCTION MatTimesMat( A, B )
 REAL, DIMENSION(:,:), INTENT(in) :: A, B
 REAL, DIMENSION( SIZE(A,1), SIZE(B,2) ) :: MatTimesMat
 MatTimesMat = MATMUL( A, B )
END FUNCTION MatTimesMat

FUNCTION MatTimesVector( A, X )
 REAL, DIMENSION(:,:), INTENT(in) :: A
 REAL, DIMENSION(:),  INTENT(in)  :: X
 REAL, DIMENSION( SIZE(A,1) ) :: MatTimesVector
 MatTimesVector = MATMUL( A, X ) 
END FUNCTION MatTimesVector

FUNCTION VectorTimesVector( Y, X )
 REAL, DIMENSION(:), INTENT(in) :: Y, X
 REAL                           :: VectorTimesVector
 VectorTimesVector = DOT_PRODUCT( Y, X ) 
END FUNCTION VectorTimesVector

FUNCTION CMatTimesCMat( A, B )
 COMPLEX, DIMENSION(:,:), INTENT(in) :: A, B
 COMPLEX, DIMENSION( SIZE(A,1), SIZE(B,2) ) :: CMatTimesCMat
 CMatTimesCMat = MATMUL( A, B )
END FUNCTION CMatTimesCMat

FUNCTION CMatTimesCVector( A, X )
 COMPLEX, DIMENSION(:,:), INTENT(in) :: A
 COMPLEX, DIMENSION(:),  INTENT(in)  :: X
 COMPLEX, DIMENSION( SIZE(A,1) ) :: CMatTimesCVector
 CMatTimesCVector = MATMUL( A, X ) 
END FUNCTION CMatTimesCVector

FUNCTION CVectorTimesCVector( Y, X )
 COMPLEX, DIMENSION(:), INTENT(in) :: Y, X
 COMPLEX                           :: CVectorTimesCVector
 CVectorTimesCVector = DOT_PRODUCT( Y, X ) 
END FUNCTION CVectorTimesCVector

END MODULE MatMult



MODULE PosMult

INTERFACE OPERATOR(.p.)
 MODULE PROCEDURE PosicTimesPosic
END INTERFACE

CONTAINS

FUNCTION PosicTimesPosic( A, B )

 REAL, DIMENSION( : , : ), INTENT(in) :: A, B
 REAL, DIMENSION( SIZE(A,1) ) :: PosicTimesPosic
 INTEGER I, EM
 EM = SIZE(A,2)
 DO I = 1, SIZE(A,1)  
   PosicTimesPosic(I) =  SUM(A(I,1:EM)  * B(I,1:EM))! scalar product
 END DO
END FUNCTION PosicTimesPosic

END MODULE PosMult

