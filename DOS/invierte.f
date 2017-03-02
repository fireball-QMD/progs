      SUBROUTINE INV(B,A,N,NDIM)
      COMPLEX A(NDIM,NDIM),B(NDIM,NDIM),X
      DO 5 I=1,N
      DO 5 J=1,N
5      A(I,J)=B(I,J)
      DO 4 I=1,N
      X=A(I,I)
      A(I,I)=1.
      DO 1 J=1,N
1      A(J,I)=A(J,I)/X
      DO 4 K=1,N
      IF(K.EQ.I) GO TO 4
      X=A(I,K)
      A(I,K)=0.
      DO 3 J=1,N
3      A(J,K)=A(J,K)-A(J,I)*X
4      CONTINUE
      RETURN
      END
