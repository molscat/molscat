      SUBROUTINE H2H2
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  FARRAR-LEE H2-H2 MSV POTENTIAL PLUS ANISOTROPIC TERMS USED BY ZARUR AND
C                                                                RABITZ.
C
C  ENTRY VINIT(I,RM,EPSIL)  INITIALIZES ROUTINE FOR LAM(I)
C                           SETS RM, EPSIL
C  ENTRY VSTAR (I,R,V)  GIVES POTENTIAL FOR LAM(I) AT R.
C  ENTRY VSTAR1(I,R,V) GIVES DERIVATIVE FOR LAM(I) AT R.
C  ENTRY VSTAR2(I,R,V) GIVES 2ND DERIV. FOR LAM(I) AT R.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION B(4),A(4,4)
      DIMENSION CONST(4)
      PARAMETER (PI=ACOS(-1.D0))
C
C  DEFINE STATEMENT FUNCTIONS
      E(R)=EXP(BETA*(1D0-R))
      P4(C1,C2,C3,C4,R)=C1+R*(C2+R*(C3+R*C4))
      P3(C1,C2,C3,R)=C1+R*(C2+R*C3)
      P2(C1,C2,R)=C1+R*C2
C
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      ENTRY VINIT(I,RM,EPSIL)
C
C  SET CONSTANT FACTORS . . .
C
      RM=3.49D0
      EPSIL=24.17D0
      RMSAVE=RM
C     CONST(1)=(4D0*3.1415926D0)**( 1.5D0)
      CONST(1)=(4D0*PI)**(1.5D0)
      CONST(2)=(.14D0/5D0)*CONST(1)
      CONST(3)=CONST(2)
      CONST(4)=54836D0 / EPSIL
C
      IF (I.EQ.4)  GOTO 1400
      IF (I.LT.1  .OR. I.GT.4)  GOTO 9999
C
      WRITE(6,601)
  601 FORMAT(' ',10X,'FARRAR-LEE MSV H2-H2 POTENTIAL /     RM=3.49,',
     1        ' EPSIL=24.17   SET INTERNALLY.')
      IF (I.GT.1)  WRITE(6,606)
  606 FORMAT(16X,'ANISOTROPIC TERM IS .14/5.0 TIMES ISOTROPIC TERM')
C
C  THE FOLLOWING PARAMETERS MUST BE SPECIFIED INTERNALLY.
C        R1,R2 IN RM .   C6,C8 IN (CM-1 )/ANG**N        BETA
C
C  SCALE PARAMETERS WITH RM, EPS.
C
      R1=1.0754D0
      R2=1.4D0
      C6=( 57913.D0 /EPSIL)/  RM**6
      C8=( 156263.D0/EPSIL)/  RM**8
      BETA=6.5D0
C
      WRITE(6,602)  R1,BETA,R2,C6,C8
  602 FORMAT(20X,'FOR R LESS THAN',F8.4,', MORSE BETA =',F9.4/
     1  20X,'FOR R GREATER THAN',F8.4,', VAN DER WAALS C6, C8 =',
     2      2F10.5,' (EPSIL/RM**N).')
C
C  SET-UP SPLINE COEFFS. B(J,I), J=1,4  FOR LAM(I).
C  ** MAR 1979 CHANGED TO INCORPORATE PRECOMPUTED SPLINE COEFFS.
C
      B(1)=-.6829 4444 3121 D1
      B(2)= .7211 3009 9412 D1
      B(3)=-.7689 9426 2915 D0
      B(4)=-.7125 2209 5040 D0
C
 1100 WRITE(6,605)  (B(J  ),J=1,4)
  605 FORMAT(20X,'CUBIC SPLINE COEFFICIENTS, TO INCREASING POWERS OF',
     1       ' R, ARE AS FOLLOWS'/25X,4E16.8)
      RETURN
C
 1400 WRITE(6,607)
  607 FORMAT(' ',10X,'QUADRUPOLE-QUADRUPOLE LONG-RANGE W/ Q=.662',
     1       ' BUCKINGHAM.')
      AC=.0687D0
      RETURN
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      ENTRY VSTAR (I,R,V)
      IF (I.EQ.4)  GOTO 2400
      IF (R.LE.R1)  GOTO 2100
      IF (R.GE.R2)  GOTO 2200
      V=P4(B(1),B(2),B(3),B(4),R)
      GOTO 5000
 2100 T1=E(R)
      V=T1*(T1-2D0)
      GOTO 5000
 2200 RSQ=1D0/(R*R)
      R6=RSQ*RSQ*RSQ
      V=-       R6*(C6+C8*RSQ)
      GOTO 5000
 2400 RR=R*RMSAVE
      R12=R**(-12)
      RRT=RR+AC*R12
      V=RRT**(-5)
      GOTO 5000
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      ENTRY VSTAR1(I,R,V)
      IF (I.EQ.4)  GOTO 3400
      IF (R.LE.R1)  GOTO 3100
      IF (R.GE.R2)  GOTO 3200
      V=P3(B(2),2D0*B(3),3D0*B(4),R)
      GOTO 5000
 3100 T1=E(R)
      V=-2D0*BETA*       T1*(T1-1D0)
      GOTO 5000
 3200 RSQ=1D0/(R*R)
      R6=RSQ*RSQ*RSQ/R
      V=       R6*(6D0*C6+8D0*C8*RSQ)
      GOTO 5000
 3400 RR=R*RMSAVE
      R12=R**(-12)
      RRT=RR+AC*R12
      V=-5D0*RRT**(-6)*(1D0-12D0*AC*R12/RR)
      GOTO 5000
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      ENTRY VSTAR2(I,R,V)
      IF (I.EQ.4)  GOTO 4400
      IF (R.LE.R1)  GOTO 4100
      IF (R.GE.R2)  GOTO 4200
      V=P2(2D0*B(3),6D0*B(4),R)
      GOTO 5000
 4100 T1=E(R)
      V=T1*BETA*BETA*       (4D0*T1-2D0)
      GOTO 5000
 4200 RSQ=1D0/(R*R)
      R6=RSQ*RSQ
      R6=R6*R6
      V=-       R6*(42D0*C6+56D0*C8*RSQ)
      GOTO 5000
 4400 RR=R*RMSAVE
      R12=R**(-12)
      RRT=RR+AC*R12
      XXX=1D0-12D0*AC*R12/RR
      V=-5D0*RRT**(-6)*(156D0*AC*R12/(RR*RR)-6D0*XXX*XXX/RRT)
      GOTO 5000
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C  UNIFIED RETURN POINT
 5000 V=V*CONST(I)
      RETURN
C
C
 9999 WRITE(6,699)  I
  699 FORMAT('  * * * ERROR.  POTENTIAL NOT DEFINED FOR SYMMETRY=',I10)
      STOP
      END
