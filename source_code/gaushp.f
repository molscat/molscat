      SUBROUTINE GAUSHP(NN,X,A)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  CALCULATES THE ZEROS, X(I), AND WEIGHTS, A(I), I=1,NN, FOR
C  GAUSS-HERMITE QUADRATURE IN ORDER TO
C  APPROXIMATE THE INTEGRAL FROM -INFINITY TO INFINITY F(X)*EXP(-X**2)
C  BY THE SUM(I=1,NN) W(I)*F(X(I)).
C  ADAPTED BY S. GREEN FROM STROUD AND SECREST GAUSSIAN QUADRATURE FORMULAS.
C  VERSION OF 18 APRIL 94; FIXED NN=1 BUG 10 MAR 95 (SG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXIT=15)
      DIMENSION X(NN),A(NN)
      DATA EPS/1.D-15/
C
      GAM(Y)=(((((((.035868343D0*Y-.193527818D0)*Y+.482199394D0)*Y-
     1             .756704078D0)*Y+.918206857D0)*Y-.897056937D0)*Y+
     2             .988205891D0)*Y-.577191652D0)*Y+1.D0
C
      IF (NN.LE.0) THEN
        WRITE(6,*) ' *** GAUSHP CALLED FOR ILLEGAL NPT=',NN
        STOP
      ELSEIF (NN.EQ.1) THEN
        WRITE(6,*) ' *** GAUSHP. WARNING, SINGLE POINT REQUESTED.'
        X(1)=0.D0
        A(1)=SQRT(ACOS(-1.D0))
        RETURN
      ELSE
        FN=NN
        N1=NN-1
        N2=(NN+1)/2
C  COMPUTE GAMMA FN BY HASTINGS APPROX; 0.LE.X.LE.70.
        Z=FN
        IF (Z.LE.0.D0 .OR. Z.GE.7.D1) THEN
          WRITE(6,600) Z
  600     FORMAT('  *** GAUSHP. CANNOT GET GAMMA FUNCTION FOR',F10.2)
          STOP
        ENDIF
        IF (Z.EQ.1.D0) THEN
          GAMMA1=1.D0
          GOTO 20
        ELSEIF (Z.LT.1.D0) THEN
          GAMMA1=GAM(Z)/Z
          GOTO 20
        ELSE
          ZA=1.D0
   10     Z=Z-1.D0
          IF (Z-1.D0) 13,11,12
   11     GAMMA1=ZA
          GOTO 20
   12     ZA=ZA*Z
          GOTO 10
   13     GAMMA1=ZA*GAM(Z)
          GOTO 20
        ENDIF
   20   CC=1.7724538509D0*GAMMA1*(2.D0**(-N1))
        S=(2.D0*FN+1.D0)**(1.D0/6.D0)
        DO 100 I=1,N2
        IF (I.EQ.1) THEN
C  LARGEST ZERO
          XT=S**3-1.85575D0/S
          GOTO 50
        ELSEIF (I.EQ.2) THEN
C  SECOND ZERO
          XT=XT-1.14D0*FN**.426D0/XT
          GOTO 50
        ELSEIF (I.EQ.3) THEN
C  THIRD ZERO
          XT=1.86D0*XT-0.86D0*X(1)
          GOTO 50
        ELSEIF (I.EQ.4) THEN
C  FOURTH ZERO
          XT=1.91D0*XT-0.91D0*X(2)
          GOTO 50
        ELSE
C  ALL HIGHER ZERO'S
          XT=2.D0*XT-X(I-2)
        ENDIF
C
C  IMPROVE THE APPROXIMATE ROOT XT AND OBTAIN
C  DPN = DERIVATIVE OF H(N) AT XT;  PN1 = VALUE OF H(N-1) AT XT
   50   IT=0
   60   IT=IT+1
        IF (IT.GT.MXIT) THEN
          WRITE(6,*) ' *** GAUSHP FAILED TO CONVERGE. ITERATIONS ='
     1               ,MXIT
          STOP
        ENDIF
        CALL HRECUR(P,DP,PN1,XT,NN)
        D=P/DP
        XT=XT-D
        IF (ABS(D).GT.EPS) GOTO 60
        DPN=DP
        X(I)=XT
        A(I)=CC/(DPN*PN1)
        NI=NN-I+1
        X(NI)=-XT
  100   A(NI)=A(I)
      ENDIF
      RETURN
      END
