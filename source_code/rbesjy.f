      SUBROUTINE RBESJY (ELL,X, V_U, DV_U, V_Y, DV_Y )
C  This subroutine is part of the MOLSCAT, BOUND and FIELD suite of programs
C
C------------------------------------------------------------------------
C  THIS ROUTINE ORIGINATES FROM THE ABC CODE BY MANOLOPOULOS ET AL.
C  COMPUTER PHYSICS COMMUNICATIONS 133 (2000) 128-135.
C------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C-----------------------------------------------------------------
C  RICCATI-BESSEL FUNCTIONS OF FRACTIONAL ORDER
C  AND THEIR FIRST DERIVATIVES WITH RESPECT TO X:
C
C  J(ELL,X) = CJ * EXP(EJ)
C  Y(ELL,X) = CY * EXP(EY)
C  D/DX J(ELL,X) = DJ * EXP(EJ)
C  D/DX Y(ELL,X) = DY * EXP(EY)
C-----------------------------------------------------------------
C

      IF (X.LE.0.0D0 .OR. ELL.LT.-0.5D0) STOP 'RBESJY 0'
      V = ELL+0.5D0
      CALL BESSJY (V,X,CJ,DJ,EJ,CY,DY,EY)
      PI = ACOS(-1.D0)
      EX = 0.5D0*DLOG(PI*X/2.D0)
      DJ = DJ+CJ/(2.D0*X)
      SJ = SQRT(CJ*CJ+DJ*DJ)
      CJ = CJ/SJ
      DJ = DJ/SJ
      EJ = EJ+DLOG(SJ)+EX
      DY = DY+CY/(2.D0*X)
      SY = SQRT(CY*CY+DY*DY)
      CY = CY/SY
      DY = DY/SY
      EY = EY+DLOG(SY)+EX

C-------------------------------------
C     WRITE(6,*) "J(X)_ELL",CJ*DEXP(EJ)
C     WRITE(6,*) "Y(X)_ELL",CY*DEXP(EY)
C     WRITE(6,*) "DJ(X)_ELL",DJ*DEXP(EJ)
C     WRITE(6,*) "DY(X)_ELL",DY*DEXP(EY)

      V_U=CJ*DEXP(EJ)    ! UJ
      DV_U=DJ*DEXP(EJ)   ! UJP
      V_Y=CY*DEXP(EY)    ! UN
      DV_Y=DY*DEXP(EY)   ! UNP

      RETURN
      END
C==============================================================================
      SUBROUTINE BESSJY (V,X,CJ,DJ,EJ,CY,DY,EY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  THIS SUBROUTINE USES A COMBINATION OF METHODS (MOSTLY DUE
C  TO TEMME) TO CALCULATE THE ORDINARY BESSEL FUNCTIONS
C
C  J(V,X) = CJ * EXP(EJ)
C  Y(V,X) = CY * EXP(EY)
C
C  AND THEIR FIRST DERIVATIVES WITH RESPECT TO X
C
C  D/DX J(V,X) = DJ * EXP(EJ)
C  D/DX Y(V,X) = DY * EXP(EY)
C
C  FOR A GIVEN REAL ORDER V >= 0 AND REAL ARGUMENT X > 0.
C  NOTE THE EXPONENTIAL SCALING, WHICH IS USED TO AVOID
C  OVERFLOW OF Y(V,X) AND UNDERFLOW OF J(V,X) FOR V >> X.
C  -----------------------------------------------------------------
C
      PARAMETER (EPS = 1.D-15)! CONSISTENT WITH RGAMMA
      PARAMETER (MAXIT = 1000)
C
      IF (V.LT.0.D0 .OR. X.LE.0.D0) STOP 'BESSJY 0'
      PI = ACOS(-1.D0)
      XMIN = 3.D0
      XMAX = 5.D0-DLOG10(EPS)
C
C  BEGIN BY CALCULATING Y(A,X) AND Y(A+1,X) FOR |A| <= 1/2
C
      NA = INT(V+0.5D0)
      A = V-NA
      IF (X.LT.XMIN) THEN
C
C  USING TEMME'S SERIES (BESSYA) FOR SMALL X
C  [ N.M.TEMME, J COMPUT PHYS 21 (1976) 343-350 ]
C
         B = X/2.D0
         D = -DLOG(B)
         E = A*D
         IF (ABS(A).LT.EPS) THEN
            C = 1.D0/PI
         ELSE
            C = A/SIN(A*PI)
         ENDIF
         IF (ABS(E).LT.EPS) THEN
            S = 1.D0
         ELSE
            S = SINH(E)/E
         ENDIF
         E = EXP(E)
         G = E*RGAMMA(A,P,Q)
         E = (E+1.D0/E)/2.D0
         F = 2*C*(P*E+Q*S*D)
         E = A*A
         P = G*C
         Q = 1.D0/G/PI
         C = A*PI/2.D0
         IF (ABS(C).LT.EPS) THEN
            R = 1.D0
         ELSE
            R = SIN(C)/C
         ENDIF
         R = PI*C*R*R
         C = 1.D0
         D = -B*B
         YA = F+R*Q
         YA1 = P
         DO N = 1,MAXIT
            F = (F*N+P+Q)/(N*N-E)
            C = C*D/N
            P = P/(N-A)
            Q = Q/(N+A)
            G = C*(F+R*Q)
            H = C*P-N*G
            YA = YA+G
            YA1 = YA1+H
            DEL = ABS(G)/(1.D0+ABS(YA))
            DEL1 = ABS(H)/(1.D0+ABS(YA1))
            IF (DEL+DEL1.LT.EPS) GOTO 1
         ENDDO
         STOP 'BESSJY 1'
   1     F = -YA
         G = -YA1/B
      ELSEIF (X.GE.XMIN .AND. X.LT.XMAX) THEN
C
C  TEMME'S PQ METHOD (BESSPQA) FOR INTERMEDIATE X
C  [ N.M.TEMME, J COMPUT PHYS 21 (1976) 343-350 ]
C
         C = 0.25D0-A*A
         B = X+X
         P = PI
         E = (X*COS(A*PI)/PI/EPS)**2
         P = 1.D0
         Q = -X
         R = 1.D0+X*X
         S = R
         DO N = 2,MAXIT
            D = (N-1+C/N)/S
            P = (2*N-P*D)/(N+1)
            Q = (-B+Q*D)/(N+1)
            S = P*P+Q*Q
            R = R*S
            IF (R*N*N.GT.E) GOTO 2
         ENDDO
         STOP 'BESSJY 2'
   2     P = P/S
         F = P
         Q = -Q/S
         G = Q
         DO M = N,1,-1
            R = (M+1)*(2.D0-P)-2.D0
            S = B+(M+1)*Q
            D = (M-1+C/M)/(R*R+S*S)
            P = D*R
            Q = D*S
            E = F+1.D0
            F = P*E-G*Q
            G = Q*E+P*G
         ENDDO
         F = 1.D0+F
         D = F*F+G*G
         PA = F/D
         QA = -G/D
         D = A+0.5D0-P
         Q = Q+X
         PA1 = (PA*Q-QA*D)/X
         QA1 = (QA*Q+PA*D)/X
         B = X-PI*(A+0.5D0)/2.D0
         C = COS(B)
         S = SIN(B)
         D = SQRT(2.D0/X/PI)
         F = D*(PA*S+QA*C)
         G = D*(QA1*S-PA1*C)
      ELSEIF (X.GE.XMAX) THEN
C
C  AND HANKEL'S ASYMPTOTIC EXPANSIONS FOR LARGE X
C  [ ABRAMOWITZ AND STEGUN, SECTION 9.2 ]
C
         P = 0.D0
         Q = 0.D0
         DO IA = 0,1
            PA = P
            QA = Q
            Y = 4.D0*(A+IA)**2
            Z = 8.D0*X
            D = 0.D0
            W = -1.D0
            P = 1.D0
            Q = 0.D0
            TP = 1.D0
            DO K = 1,MAXIT
               D = D+Z
               W = W+2.D0
               TQ = +TP*(Y-W*W)/D
               Q = Q+TQ
               D = D+Z
               W = W+2.D0
               TP = -TQ*(Y-W*W)/D
               P = P+TP
               IF (ABS(TP)+ABS(TQ).LT.EPS) GOTO 3
            ENDDO
            STOP 'BESSJY 3'
   3        P = P-0.5D0*TP
            Q = Q-0.5D0*TQ
         ENDDO
         PA1 = P
         QA1 = Q
         B = X-PI*(A+0.5D0)/2.D0
         C = COS(B)
         S = SIN(B)
         D = SQRT(2.D0/X/PI)
         F = D*(PA*S+QA*C)
         G = D*(QA1*S-PA1*C)
      ENDIF
C
C  NOW RECUR UPWARDS FROM Y(A,X) TO Y(V,X),
C  SCALING TO AVOID OVERFLOW ALONG THE WAY
C
      P = 0.D0
      IF (NA.GT.0) THEN
         Y = 2.D0/X
         DO N = 1,NA
            H = Y*(A+N)*G-F
            F = G
            G = H
   4        IF (ABS(F).GT.4.D0) THEN
               P = P+1.D0
               F = 0.0625D0*F
               G = 0.0625D0*G
               GOTO 4
            ENDIF
         ENDDO
      ENDIF
      CY = F
      DY = (V/X)*F-G
      SY = SQRT(CY*CY+DY*DY)
      CY = CY/SY
      DY = DY/SY
      EY = DLOG(SY)+P*DLOG(16.D0)
C
C  FINALLY, CALCULATE J(V,X) AND DJ(V,X)/DX
C
      VV = MAX(XMIN,V)
      IF (X.GE.VV) THEN
C
C  USING UPWARD RECURSION IN THE CLASSICALLY ALLOWED REGION
C
         F = D*(PA*C-QA*S)
         G = D*(QA1*C+PA1*S)
         IF (NA.GT.0) THEN
            Y = 2.D0/X
            DO N = 1,NA
               H = Y*(A+N)*G-F
               F = G
               G = H
            ENDDO
         ENDIF
         CJ = F
         DJ = (V/X)*F-G
         SJ = SQRT(CJ*CJ+DJ*DJ)
         CJ = CJ/SJ
         DJ = DJ/SJ
         EJ = DLOG(SJ)
      ELSE
C
C  AND CF1 IN THE CLASSICALLY FORBIDDEN REGION
C  [ NUMERICAL RECIPES, 2ND EDITION, SECTION 6.7 ]
C
         AP = 1.D0
         A = V/X
         BP = 0.D0
         B = 1.D0
         F = 0.D0
         G = 0.D0
         Y = 2.D0/X
         W = Y/PI
         DO N = 1,MAXIT
            AN = Y*(V+N)*A-AP
            AP = A
            A = AN
            BN = Y*(V+N)*B-BP
            BP = B
            B = BN
            IF (ABS(B).GT.ABS(A)) THEN
               AP = AP/B
               A = A/B
               BP = BP/B
               B = 1.D0
               IF (ABS(A-F).LT.EPS*ABS(F)) THEN
                  CJ = W/(DY-CY*A)
                  DJ = A*CJ
                  GOTO 5
               ENDIF
               F = A
            ELSE
               BP = BP/A
               B = B/A
               AP = AP/A
               A = 1.D0
               IF (ABS(B-G).LT.EPS*ABS(G)) THEN
                  DJ = W/(DY*B-CY)
                  CJ = B*DJ
                  GOTO 5
               ENDIF
               G = B
            ENDIF
         ENDDO
         STOP 'BESSJY 4'
   5     SJ = SQRT(CJ*CJ+DJ*DJ)
         CJ = CJ/SJ
         DJ = DJ/SJ
         EJ = DLOG(SJ)-EY
      ENDIF
      RETURN
      END
C==============================================================================
      FUNCTION RGAMMA(X,ODD,EVEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  -----------------------------------------------------------------
C  DIRECT FORTRAN TRANSLATION OF TEMME'S ALGOL ROUTINE FOR COMPUTING
C  RGAMMA = 1/GAMMA(1-X), ALONG WITH ITS ODD AND EVEN PARTS, FOR
C  ABS(X) .LE. 0.5. [ N.M.TEMME, J COMPUT PHYS 19 (1975) 324-337 ]
C  -----------------------------------------------------------------
C
      DIMENSION B(12)
      DATA B / -0.283876542276024D0, -0.076852840844786D0,
     *         +0.001706305071096D0, +0.001271927136655D0,
     *         +0.000076309597586D0, -0.000004971736704D0,
     *         -0.000000865920800D0, -0.000000033126120D0,
     *         +0.000000001745136D0, +0.000000000242310D0,
     *         +0.000000000009161D0, -0.000000000000170D0 /
      SAVE B
C
      X2 = X*X*8.D0
      ALFA = -0.000000000000001D0
      BETA = 0.D0
      DO I = 12,2,-2
         BETA = -(2*ALFA+BETA)
         ALFA = -BETA*X2-ALFA+B(I)
      ENDDO
      EVEN = (BETA/2.D0+ALFA)*X2-ALFA+0.921870293650453D0
      ALFA = -0.000000000000034D0
      BETA = 0.D0
      DO I = 11,1,-2
         BETA = -(2*ALFA+BETA)
         ALFA = -BETA*X2-ALFA+B(I)
      ENDDO
      ODD = 2*(ALFA+BETA)
      RGAMMA = ODD*X+EVEN
      RETURN
      END
