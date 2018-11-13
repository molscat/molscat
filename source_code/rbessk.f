      SUBROUTINE RBESSK (ELL,X,VRATIO)
C  This subroutine is part of the MOLSCAT, BOUND and FIELD suite of programs
C------------------------------------------------------------------------------
C  THIS ROUTINE ORIGINATES FROM THE ABC CODE BY MANOLOPOULOS ET AL.
C  COMPUTER PHYSICS COMMUNICATIONS 133 (2000) 128-135.
C------------------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  -----------------------------------------------------------------
C  MODIFIED RICCATI-BESSEL FUNCTION OF THE THIRD KIND
C  AND ITS FIRST DERIVATIVE WITH RESPECT TO X:
C
C  K(ELL,X) = CK * EXP(EK)
C  D/DX K(ELL,X) = DK * EXP(EK)
C  -----------------------------------------------------------------
C

      IF (X.LE.0.0D0 .OR. ELL.LT.-0.5D0) STOP 'RBESSK 0'
      V = ELL+0.5D0
      CALL MBESSK (V,X,CK,DK,EK)
      PI = ACOS(-1.D0)
      EX = 0.5D0*DLOG(PI*X/2.D0)
      DK = DK+CK/(2.D0*X)
      SK = SQRT(CK*CK+DK*DK)
      CK = CK/SK
      DK = DK/SK
      EK = EK+DLOG(SK)+EX

C-------------------------------------------------------
C      WRITE(6,*) "K(X)_ELL",CK*DEXP(EK)
C      WRITE(6,*) "DK(X)_ELL",DK*DEXP(EK)
C      WRITE(6,*) "RATIO D/C(KOCCHI), C/D", DK/CK, CK/DK
       VRATIO=DK/CK
C-------------------------------------------------------

      RETURN
      END
C==============================================================================
      SUBROUTINE MBESSK (V,X,CK,DK,EK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  -----------------------------------------------------------------
C  THIS SUBROUTINE USES TEMME'S METHOD [ N.M.TEMME, J COMPUT PHYS
C  19 (1975) 324-337 ] TO CALCULATE THE MODIFIED BESSEL FUNCTION
C
C  K(V,X) = CK * EXP(EK)
C
C  AND ITS FIRST DERIVATIVE WITH RESPECT TO X
C
C  D/DX K(V,X) = DK * EXP(EK)
C
C  FOR A GIVEN REAL ORDER V >= 0 AND REAL ARGUMENT X > 0.
C  NOTE THE EXPONENTIAL SCALING, WHICH IS USED TO AVOID
C  OVERFLOW OF K(V,X) FOR V >> X AND UNDERFLOW FOR V << X.
C  -----------------------------------------------------------------
C
      PARAMETER (EPS = 1.D-15)! CONSISTENT WITH RGAMMA
      PARAMETER (MAXIT = 1000)
C
      IF (V.LT.0.D0 .OR. X.LE.0.D0) STOP 'MBESSK 0'
      PI = ACOS(-1.D0)
      XMIN = 1.D0
C
C  BEGIN BY CALCULATING K(A,X) AND K(A+1,X) FOR |A| <= 1/2
C
      NA = INT(V+0.5D0)
      A = V-NA
      IF (X.LT.XMIN) THEN
C
C  USING TEMME'S SERIES FOR SMALL X
C
         B = X/2.D0
         D = -DLOG(B)
         E = A*D
         C = A*PI
         IF (ABS(C).LT.EPS) THEN
            C = 1.D0
         ELSE
            C = C/SIN(C)
         ENDIF
         IF (ABS(E).LT.EPS) THEN
            S = 1.D0
         ELSE
            S = SINH(E)/E
         ENDIF
         E = EXP(E)
         G = E*RGAMMA(A,P,Q)
         E = (E+1.D0/E)/2.D0
         F = C*(P*E+Q*S*D)
         E = A*A
         P = 0.5D0*G*C
         Q = 0.5D0/G
         C = 1.D0
         D = B*B
         AK = F
         AK1 = P
         DO N = 1,MAXIT
            F = (F*N+P+Q)/(N*N-E)
            C = C*D/N
            P = P/(N-A)
            Q = Q/(N+A)
            G = C*(P-N*F)
            H = C*F
            AK = AK+H
            AK1 = AK1+G
            IF (H/AK+ABS(G)/AK1.LT.EPS) GOTO 1
         ENDDO
         STOP 'MBESSK 1'
   1     F = AK
         G = AK1/B
         EX = 0.D0
      ELSEIF (X.GE.XMIN) THEN
C
C  AND TEMME'S PQ METHOD FOR LARGE X
C
         C = 0.25D0-A*A
         G = 1.D0
         F = 0.D0
         E = X*COS(A*PI)/PI/EPS
         DO N = 1,MAXIT
            H = (2*(N+X)*G-(N-1+C/N)*F)/(N+1)
            F = G
            G = H
            IF (H*N.GT.E) GOTO 2
         ENDDO
         STOP 'MBESSK 2'
   2     P = F/G
         Q = P
         B = X+X
         E = B-2.D0
         DO M = N,1,-1
            P = (M-1+C/M)/(E+(M+1)*(2.D0-P))
            Q = P*(Q+1.D0)
         ENDDO
         F = SQRT(PI/B)/(1.D0+Q)
         G = F*(A+X+0.5D0-P)/X
         EX = X
      ENDIF
C
C  NOW RECUR UPWARDS FROM K(A,X) TO K(V,X),
C  SCALING TO AVOID OVERFLOW ALONG THE WAY
C
      P = 0.D0
      IF (NA.GT.0) THEN
         Y = 2.D0/X
         DO N = 1,NA
            H = Y*(A+N)*G+F
            F = G
            G = H
   3        IF (ABS(F).GT.4.D0) THEN
               P = P+1.D0
               F = 0.0625D0*F
               G = 0.0625D0*G
               GOTO 3
            ENDIF
         ENDDO
      ENDIF
      CK = F
      DK = (V/X)*F-G
      SK = SQRT(CK*CK+DK*DK)
      CK = CK/SK
      DK = DK/SK
      EK = DLOG(SK)+P*DLOG(16.D0)-EX
      RETURN
      END
