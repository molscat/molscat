      SUBROUTINE DSYNEG(UPLO,A,KPVT,N,INERT)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  SUBROUTINE BY J. M. HUTSON, APRIL 1994
C  BASED ON LINPACK ROUTINE DSIDI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UPPER,LSAME
      CHARACTER(1) UPLO
      DIMENSION A(N,N), KPVT(N)
C
C  SUBROUTINE TO FIND THE NUMBER OF NEGATIVE EIGENVALUES OF A
C  SYMMETRIC INDEFINITE MATRIX AFTER BUNCH-KAUFMAN DECOMPOSITION
C  USING LAPACK ROUTINE DSYTRF
C
      UPPER = LSAME (UPLO , 'U')
      IF (UPPER) THEN
         INERT=0
         T = 0.0D0
         DO 190 K = 1, N
            D = A(K,K)
C
C  CHECK IF 1 BY 1
C
            IF (KPVT(K).GT.0) GOTO 150
C
C  2 BY 2 BLOCK
C  USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S)
C          (S  C)
C  TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C  TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (T.NE.0.0D0) GOTO 130
                  T = ABS(A(K,K+1))
                  D = (D/T)*A(K+1,K+1) - T
               GOTO 140
  130          CONTINUE
                  D = T
                  T = 0.0D0
  140          CONTINUE
  150       CONTINUE
C
               IF (D.LT.0.0D0) INERT = INERT + 1
  190    CONTINUE
      ELSE
         INERT=0
         T = 0.0D0
         DO 290 K = N, 1, -1
            D = A(K,K)
C
C  CHECK IF 1 BY 1
C
            IF (KPVT(K).GT.0) GOTO 250
C
C  2 BY 2 BLOCK
C  USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S)
C          (S  C)
C  TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C  TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (T.NE.0.0D0) GOTO 230
                  T = ABS(A(K,K-1))
                  D = (D/T)*A(K-1,K-1) - T
               GOTO 240
  230          CONTINUE
                  D = T
                  T = 0.0D0
  240          CONTINUE
  250       CONTINUE
C
               IF (D.LT.0.0D0) INERT = INERT + 1
  290    CONTINUE
      ENDIF
C
      RETURN
      END
