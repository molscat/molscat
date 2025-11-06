      SUBROUTINE SYMINV(A, IA, N, IFAIL)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,N)
C
C  ON A TYPICAL SCALAR MACHINE, A LAPACK-BASED SYMINL IS FASTER
C  THAN THE INLINE CODE HERE FOR MATRICES BIGGER THAN ABOUT 35*35.
C  THESE 4 LINES COULD BE COMMENTED OUT IF LAPACK IS NOT AVAILABLE
C
      IF (N.GT.30) THEN
        CALL SYMINL(A, IA, N, IFAIL)
        RETURN
      ENDIF
C
C  IN SITU INVERSION OF A REAL SYMMETRIC MATRIX BY L.D.L(T)
C  DECOMPOSITION. KOUNT = NUMBER OF NEGATIVE EIGENVALUES.
C  ROUTINE FAILS (IFAIL = N+1) IF THE DECOMPOSITION IS UNDEFINED.
C  J.H.WILKINSON AND C.REINSCH, LINEAR ALGEBRA, I/1
C
      KOUNT = 0
      X = A(1,1)
      IF (X.LT.0.D0) KOUNT = KOUNT + 1
      IF (X.EQ.0.D0) GOTO 320
      A(1,1) = 1.D0/X
      IF (N.EQ.1) GOTO 300
C  FORM L
      DO I = 2,N
         I1 = I - 1
         IF (I1.LT.2) GOTO 60
         DO J = 2,I1
            J1 = J - 1
            X = A(I,J)
            DO K = 1,J1
               X = X - A(I,K)*A(J,K)
            ENDDO
            A(I,J) = X
         ENDDO
   60    X = A(I,I)
         DO K = 1,I1
            Y = A(I,K)
            Z = Y*A(K,K)
            A(I,K) = Z
            X = X - Y*Z
         ENDDO
         IF (X.LT.0.D0) KOUNT = KOUNT + 1
         IF (X.EQ.0.D0) GOTO 320
         A(I,I) = 1.D0/X
      ENDDO
C  INVERT L
      N1 = N - 1
      DO I = 1,N1
         I1 = I + 1
         A(I1,I) = - A(I1,I)
         IF (I1.GT.N1) CYCLE
         DO J = I1,N1
            J1 = J + 1
            X = - A(J1,I)
            DO K = I1,J
               X = X - A(J1,K)*A(K,I)
            ENDDO
            A(J1,I) = X
         ENDDO
      ENDDO
C  FORM INVERSE OF A
      DO I = 1,N1
         I1 = I + 1
         X = A(I,I)
         DO K = I1,N
            Y = A(K,I)
            Z = Y*A(K,K)
            A(K,I) = Z
            X = X + Y*Z
         ENDDO
         A(I,I) = X
         IF (I1.GT.N1) CYCLE
         DO J = I1,N1
            J1 = J + 1
            X = A(J,I)
            DO K = J1,N
               X = X + A(K,J)*A(K,I)
            ENDDO
            A(J,I) = X
         ENDDO
      ENDDO
C  COPY INVERSE INTO UPPER TRIANGLE
      DO I = 2,N
         I1 = I - 1
         DO J = 1,I1
            A(J,I) = A(I,J)
         ENDDO
      ENDDO
  300 IFAIL = KOUNT
      RETURN
  320 IFAIL = N + 1
      RETURN
      END
