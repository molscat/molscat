      FUNCTION THREEJ(J1,J2,J3)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  COMPUTATION OF SPECIAL WIGNER 3J COEFFICIENT WITH
C  VANISHING PROJECTIONS.  SEE EDMONDS, P. 50.
C
C  THIS VERSION EVALUATES BINOM AND PARITY IN-LINE
C  SHOULD IMPROVE EFFICIENCY, ESPECIALLY ON CRAY;
C  ALSO GIVES IMPROVEMENT ON AMDAHL  (S Green: 20 DEC 92)
C
C  STATEMENT FUNCTION FOR DELTA ASSOCIATED W/ RACAH AND SIXJ SYMBOLS
C     DELTA(I,J,K)= SQRT(1.D0/ ( BINOM(I+J+K+1,I+J-K) *
C    1                 BINOM(K+K+1,I-J+K) * DBLE(K+J-I+1) )  )
C
      I1=J1+J2+J3
      IF (I1-2*(I1/2).NE.0) GOTO 8
    1 I2=J1-J2+J3
      IF (I2) 8,2,2
    2 I3=J1+J2-J3
      IF (I3) 8,3,3
    3 I4=-J1+J2+J3
      IF (I4) 8,4,4
    4 I5=I1/2
      I6=I2/2

      SIGN3J=1.D0
      IF (I5-2*(I5/2).NE.0) SIGN3J=-SIGN3J
C   7 THREEJ=SIGN3J*DELTA(J1,J2,J3)*BINOM(I5,J1)*BINOM(J1,I6)
C     B1,B2 ARE BINOM ASSOCIATED W/ DELTA

      N=J1+J2+J3+1
      M=J1+J2-J3
      NM = N-M
      MNM = MIN(NM,M)
      IF (MNM.LE.0) THEN
        B1=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 101 I = 1,MNM
          F = F+1.D0
          C = (FN-F)*B
  101     B = C/F
        B1 = B
      ENDIF

      N=J3+J3+1
      M=J1-J2+J3
      NM = N-M
      MNM = MIN(NM,M)
      IF (MNM.LE.0) THEN
        B2=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 102 I = 1,MNM
          F = F+1.D0
          C = (FN-F)*B
  102     B = C/F
        B2 = B
      ENDIF
      DELTA=SQRT(1.D0/(B1*B2*(J3+J2-J1+1)))
C     B3=BINOM(I5,J1),  B4=BINOM(J1,I6)

      N=I5
      M=J1
      NM = N-M
      MNM = MIN(NM,M)
      IF (MNM.LE.0) THEN
        B3=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 103 I = 1,MNM
          F = F+1.D0
          C = (FN-F)*B
  103     B = C/F
        B3 = B
      ENDIF

      N=J1
      M=I6
      NM = N-M
      MNM = MIN(NM,M)
      IF (MNM.LE.0) THEN
        B4=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 104 I = 1,MNM
          F = F+1.D0
          C = (FN-F)*B
  104     B = C/F
        B4 = B
      ENDIF

      THREEJ=SIGN3J*DELTA*B3*B4
      RETURN

    8 THREEJ=0.D0
      RETURN
      END
