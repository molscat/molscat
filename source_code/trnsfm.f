      SUBROUTINE TRNSFM(T,W,A,N,ISTOP,ISYM)
C  This subroutine is part of the MOLSCAT, BOUND and FIELD suite of programs
C-------------------------------------------------------------------
C  WRITTEN BY G. A. PARKER.
C  MODIFIED TO USE BLAS BY J. M. HUTSON
C  THIS ROUTINE TRANSFORMS THE MATRIX W INTO A NEW BASIS SET
C  ISTOP=.TRUE.   ==>  RETURN AFTER A = TRANSPOSE(W) * T
C  ISTOP=.FALSE.  ==>  CONTINUE TO FORM W = TRANSPOSE(T) * W * T
C  ISYM =.TRUE.   ==>  FORCE THE RESULTING MATRIX TO BE SYMMETRIC.
C  N IS THE DIMENSION OF THE MATRICES.
C-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ISTOP,ISYM
      DIMENSION T(1),W(1),A(1)
      DATA ZERO/0.D0/,HALF/0.5D0/,ONE/1.D0/
C-------------------------------------------------------------------
C  MULTIPLY THE TRANSPOSE OF THE MATRIX W TIMES T AND
C  STORE THE RESULT INTO MATRIX A.
C-------------------------------------------------------------------
      IF (N.EQ.1) GOTO 300
      IF (ISYM) GOTO 140
      CALL DGEMUL(W,N,'T',T,N,'N',A,N,
     1            N,N,N)
      IF (ISTOP) RETURN
C-------------------------------------------------------------------
C  MULTIPLY THE TRANSPOSE OF MATRIX A TIMES MATRIX
C  T AND STORE THE RESULT INTO MATRIX W
C-------------------------------------------------------------------
      CALL DGEMUL(A,N,'T',T,N,'N',W,N,
     1            N,N,N)
      RETURN
C-------------------------------------------------------------------
C  THIS IS REACHED ONLY WHEN W AND THE RESULT MATRIX ARE SYMMETRIC,
C  SO THAT ONLY HALF THE MATRIX NEED BE COMPUTED
C  AND THE OTHER HALF STORED BY SYMMETRY.
C-------------------------------------------------------------------
  140 CALL DSYMM('L','L',N,N,ONE,W,N,T,N,ZERO,A,N)
      IF (ISTOP) RETURN
      CALL DSYR2K('L','T',N,N,HALF,A,N,T,N,ZERO,W,N)
      CALL DSYFIL('U',N,W,N)
      RETURN
C
  300 A(1)=W(1)*T(1)
      IF (ISTOP) RETURN
      W(1)=A(1)*T(1)
      RETURN
      END
