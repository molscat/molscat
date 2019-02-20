      SUBROUTINE PERT2(N,CX,SX,XSQ,X,H,V0,V1,V2,EYE11,
     1                 EYE12,EYE22,A1,A2,A1P,A2P)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------------------------
C  THIS ROUTINE CALCULATES AND STORES THE PERTURBATION
C  INTEGRALS IN EYE11, EYE12 AND EYE22 FOR THE INTERVAL.
C  ALSO ON RETURN V0, V1 AND V2 CONTAIN THE PERTURBATION
C  INTEGRALS FOR THE CURRENT STEP WHICH ARE USED FOR
C  THE STEP SIZE DETERMINATION ONLY.
C-------------------------------------------------------------------
      LOGICAL IFIRST
      DIMENSION CX(1),SX(1),XSQ(1),X(1)
      DIMENSION V0(N,N),V1(N,N),V2(N,N)
      DIMENSION EYE11(N,N),EYE22(N,N)
      DIMENSION EYE12(N,N),A1(N),A1P(N),A2(N),A2P(N)
      DATA TOL/1.D-3/
      DATA IFIRST/.FALSE./
      SAVE
C-------------------------------------------------------------------
C  STORE CONSTANTS ON THE FIRST CALL.
C-------------------------------------------------------------------
      IF (IFIRST) GOTO 100
      CON1 = 1.D0/24.D0
      CON2 = 1.D0/3.D0
      CON3 = 1.D0/120.D0
      CON4 = 1.D0/12.D0
      CON5 = 0.2D0
      CON6 = 1.D0/5040.D0
      CON7 = 1.D0/30240.D0
      CON8 = 1.D0/15.D0
      CON9 = 1.D0/840.D0
      CON11 = 1.D0/180.D0
      CON12 = 1.D0/6.D0
      CON13 = 1.D0/560.D0
      CON10 = 1.D0/20160.D0
      CON14 = 1.D0/42.D0
      CON15 = 1.D0/8.D0
      CON16 = 1.D0/48.D0
      CON17 = 1.D0/168.D0
      CON18 = 1.D0/112.D0
      CON19 = 1.D0/40.D0
      CON20 = 1.D0/36.D0
      CON21 = 1.D0/216.D0
      IFIRST = .TRUE.
  100 HSQ = H*H
      HINV = 1.D0/H
      HSQINV = 1.D0/HSQ
      HF = H**3
C-------------------------------------------------------------------
C  THE FOLLOWING IS USED WHEN CORRECTIONS TO THE SECOND
C  DERIVATIVE OF THE POTENTIAL ARE DESIRED.
C-------------------------------------------------------------------
      DO 190 I = 1,N
        XI = X(I)
        SXI = SX(I)
        CXI = CX(I)
        XSXI = XI*SXI
        XSQI = XSQ(I)
        G1 = XSXI/XSQI
        G2 = 8.D0*XSQI
        G3 = XSQI*XSQI
        G4 = 4.D0*XSQI
        G5 = 2.D0*XSQI
        G6 = 0.5D0*G1
        G7 = G3*XSQI
      DO 190 J = I,N
        XSQJ = XSQ(J)
        IF (I.EQ.J) GOTO 160

        XJ = X(J)
        SXJ = SX(J)
        CXJ = CX(J)
        XSXJ = XJ*SXJ
        F5 = XSQI*XSQJ
        F6 = XSXI*XSXJ
        F7 = CXI*CXJ
        SD = XSQI-XSQJ
        TOTAL = XSQI+XSQJ
        COEF = HF/(2.D0*F6)
        HHC = HSQ*COEF*V2(I,J)
        HC = H*COEF*V1(I,J)
        COEF = COEF*V0(I,J)
        D = XI-XJ
        IF (ABS(D).GT.TOL) GOTO 130

        IF (ABS(TOTAL).LT.TOL) GOTO 150

        S = XI+XJ
        SSQ = S*S
        DSQ = D*D
        SHALF = S*0.5D0
        IF (XSQI.LT.0.D0 .AND. XSQJ.LT.0.D0) GOTO 110

        IF (F5.LT.0.D0) GOTO 130
        C2 = COS(SHALF)
        S2 = SIN(SHALF)
        GOTO 120

  110   C2 = COSH(SHALF)
        S2 = -SINH(SHALF)
        SSQ = -SSQ
        DSQ = -DSQ
  120   SS2 = S*S2
        GOTO 140

C-------------------------------------------------------------------
C  GENERAL FORMULAS
C-------------------------------------------------------------------
  130   F25 = 1.D0/SD
        F8 = XSXJ*CXI
        F9 = XSXI*CXJ
        F10 = SD*SD
        F11 = 2.D0*XSQJ
        F12 = 2.D0*COEF*F25
        F13 = HC*F25
        F18 = 2.D0*TOTAL
        F19 = 0.25D0*F25
        F24 = -8.D0*(TOTAL+F11)
        F26 = -8.D0*(TOTAL+G5)
        F27 = 1.D0/F10
        F31 = HHC*F27
        R1V0 = F12*(F9-F8)
        R2V0 = F12*(XSXI-XSXJ)
        R1V1 = F13*((F18*(F7-1.D0)+4.D0*F6)*F25+F9-F8)
        R2V1 = F13*(F18*(CXI-CXJ)*F25+XSXI+XSXJ)
        R1V2 = F31*((((F10+F24)*F9-(F10+F26)*F8)*F19+TOTAL*(F7+1.D0)
     1                +2.D0*F6))
        R2V2 = F31*(((F10+F24)*XSXI-(F10+F26)*XSXJ)*F19+TOTAL*(CXI+CXJ))
        GOTO 170

C-------------------------------------------------------------------
C  THE FOLLOWING FORMULAS ARE VALID WHEN THE DIFFERENCE BETWEEN
C  THE WAVEVECTORS TIMES THE STEP SIZE IS SMALL.
C-------------------------------------------------------------------
  140   F14 = SS2/SSQ
        F16 = DSQ*CON1
        F17 = DSQ/80.D0
        F15 = SS2*C2
        F36 = DSQ*CON17
        F37 = DSQ*CON18
        F38 = DSQ*CON19
        F39 = 1.D0-DSQ*CON15*(1.D0-DSQ*CON16*(1.D0-DSQ*CON3))
        F40 = 1.D0/SSQ
        F41 = (SSQ-8.D0)*0.5D0*F14
        R1V0 = COEF*(1.D0-DSQ*CON12*(1.D0-DSQ/20.D0*(1.D0-DSQ*CON14))+
     X               2.D0*F15*F40)
        R2V0 = COEF*(F14*(F39)+0.5D0*C2*(1.D0-F16*(1.D0-F17*(1.D0-F36)))
     X                 )*2.D0
        R1V1 = HC*((-2.D0*SS2*F14+F15)*F40-F16*(1.D0-DSQ*CON8*(1.D0
     1                                                      -3.D0*F37)))
        R2V1 = HC*(-S2*CON4*(1.D0-F38*(1.D0-F37))+(0.5D0*C2-F14)/S
     1                                         *(1.D0-F16*(1.D0-F17)))*D
        R1V2 = HHC*(C2*(2.D0*C2+F41)*F40+(1.D0-DSQ*CON5*(1.D0-11.D0*F36*
     X                                    (1.D0-DSQ*CON20)))*CON4)*0.5D0
        R2V2 = HHC*((F39)*(C2+F41*0.5D0)*F40+(1.D0-F38*(3.D0-F37*(5.D0-
     X                                        7.D0*DSQ*CON21)))*C2*CON1)
        GOTO 170
C-------------------------------------------------------------------
C  THE FOLLOWING FORMULAS ARE VALID WHEN BOTH OF THE WAVEVECTORS
C  TIMES THE STEP SIZE ARE SMALL.
C-------------------------------------------------------------------
  150   F20 = TOTAL*TOTAL
        F21 = TOTAL*CON5
        F22 = G4*XSQJ
        F23 = XSQJ*XSQJ
        F29 = TOTAL*CON12
        F28 = 14.D0*F5
        F30 = COEF*2.D0
        R1V0 = F30*(1.D0-F29+(F20+F22)*CON3-TOTAL*(G3+F28+F23)*CON6)
        R2V0 = F30*(1.D0-F29+(XSQI*TOTAL+F23)*CON3-(G3+F23)*TOTAL*CON6)
        R1V1 = -HC*CON4*(TOTAL-(F20+F22)*CON8+TOTAL*(G3+F28+F23)*CON13)
        R2V1 = HC*(CON4+TOTAL*CON11-(3.D0*F20-G5*XSQJ)*CON10)*SD
        R1V2 = HHC*(1.D0-F21+11.D0*(F20+F22)*CON9
     1                              -11.D0*TOTAL*(G3+F28+F23)*CON7)*CON4
        R2V2 = HHC*(1.D0-F21+(11.D0*F20-19.D0*F5)*CON9-(11.D0*(G7+F23*
     X                                   XSQJ)+3.D0*F5*TOTAL)*CON7)*CON4
        GOTO 170
C-------------------------------------------------------------------
C  FORMULAS VALID FOR DIAGONAL ELEMENTS ONLY.
C-------------------------------------------------------------------
  160   F6 = XSXI*XSXI
        CXJ = CXI
        XSXJ = XSXI
        F7 = CXI*CXI
        COEF = HF/(2.D0*F6)
        HHC = HSQ*COEF*V2(I,J)
        HC = H*COEF*V1(I,J)
        COEF = COEF*V0(I,J)
        R1V0 = COEF*(1.D0+CXI*G1)
        R2V0 = COEF*(G1+CXI)
        R1V1 = HC*(CXI-G1)*G6
        R2V1 = 0.D0
        R1V2 = HHC*(CXI*(CXI+(XSQI-2.D0)*G6)/XSQI*0.25D0+CON1)
        R2V2 = HHC*((6.D0+XSQI)*CXI*CON2+(XSQI-2.D0)*G1)/G2
C-------------------------------------------------------------------
C  CC- COSINE-COSINE INTEGRALS
C  CS- COSINE-SINE INTEGRALS
C  SS- SINE-SINE INTEGRALS
C-------------------------------------------------------------------
  170   TERM1 = -CXJ*(R1V0+R1V1+R1V2)
        CCIJ = F6*(R1V0+R1V1+R1V2)*HSQINV
        CSIJ = XSXI*(R2V0+R2V1+R2V2+TERM1)*HINV
        SSIJ = (R1V0-R1V1+R1V2-CXJ*(R2V0-R2V1+R2V2)-CXI*(R2V0+R2V1+R2V2)
     X                                             +F7*(R1V0+R1V1+R1V2))
        CCJI = CCIJ
        CSJI = XSXJ*(R2V0-R2V1+R2V2-CXI*(R1V0+R1V1+R1V2))*HINV
        SSJI = SSIJ
        V0(I,J) = CCIJ
        V0(J,I) = CCJI
        V1(I,J) = CSIJ
        V1(J,I) = CSJI
        V2(I,J) = SSIJ
        V2(J,I) = SSJI
        TERM1 = A1P(J)*SSIJ+A1(J)*CSJI
        TERM2 = A1P(J)*CSIJ+A1(J)*CCIJ
        TERM3 = A2P(J)*SSIJ+A2(J)*CSJI
        TERM4 = A2P(J)*CSIJ+A2(J)*CCIJ
        EYE11(I,J) = A1P(I)*TERM1+A1(I)*TERM2+EYE11(I,J)
        EYE11(J,I) = EYE11(I,J)
        EYE22(I,J) = A2P(I)*TERM3+A2(I)*TERM4+EYE22(I,J)
        EYE22(J,I) = EYE22(I,J)
        EYE12(I,J) = A1P(I)*TERM3+A1(I)*TERM4+EYE12(I,J)
        IF (I.EQ.J) GOTO 190
        EYE12(J,I) = A2P(I)*TERM1+A2(I)*TERM2+EYE12(J,I)
  190 CONTINUE
      RETURN
C----------------***END-PERT2***-------------------------------------
      END
