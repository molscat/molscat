      SUBROUTINE GAUSSP(A,B,NPT,XPT,WHT)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XPT(NPT),WHT(NPT)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                                                                     *
C   THIS ROUTINE SETS UP ABSCISSAE AND WEIGHTS FOR NPT-POINT          *
C   GAUSS-LEGENDRE INTEGRATION IN THE INTERVAL (A,B).                 *
C                                                                     *
C   ON RETURN, THE FUNCTION TO BE INTEGRATED SHOULD BE EVALUATED      *
C   AT THE POINTS XPT(I). INTEGRAL = SUM(I=1,NPT) F(XPT(I))*WHT(I)    *
C                                                                     *
C   THIS VERSION (SG 11/7/91) CALCULATES POINTS/WEIGHTS FROM          *
C   GASLEG/ZBES CODE OF AD VAN DER AVOIRD                             *
C   IT DOES ANY NUMBER OF PTS FROM 1 TO MXPT, WHERE LIMIT IS FROM     *
C   DIMENSION STATEMENTS IN GASLEG (P,PD AT LEAST (MXPT+1) )          *
C   AND HERE W,X DIMENSIONED AT LEAST  ((MXPT+1)/2)                   *
C                                                                     *
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DIMENSION X(128),W(128)
      DATA MXPT/256/
C
      T1=(B-A)/2.D0
      T2=(B+A)/2.D0
      IF (NPT-1) 9999,9998,9997

 9997 IF (NPT.LE.MXPT) GOTO 3100
      WRITE(6,601) NPT,MXPT
  601 FORMAT(/' * * * WARNING.  GAUSS-LEGENDRE NPT =',I6,
     1       '  REDUCED TO',I4)
      NPT=MXPT
 3100 CALL GASLEG(NPT,X,W)
      N2=(NPT+1)/2
      I1=1
      I2=NPT
      IC=1
      DO 2000 I=1,N2
        XPT(I1)=-X(IC)*T1+T2
        XPT(I2)=X(IC)*T1+T2
        WHT(I1)=W(IC)*T1
        WHT(I2)=WHT(I1)
        I1=I1+1
        I2=I2-1
 2000   IC=IC+1
C  N.B FOR NPT ODD, THE LAST (I.E. MIDDLE) TERM IS EVALUATED TWICE.
      RETURN

 9999 WRITE(6,610) NPT
  610 FORMAT(/' * * * WARNING.  GAUSS-LEGENDRE REQUESTED WITH NPT =',I6)
C  REPLACE WITH SINGLE-POINT AT (A+B)/2 * (B-A)
      NPT=1

 9998 XPT(1)=T2
      WHT(1)=2.D0*T1
      RETURN
      END
