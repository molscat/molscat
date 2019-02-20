      FUNCTION RSYMTP(J1,K1,J2,J1P,K1P,J2P,
     1                JJ,JJP,MU,P1,Q1,P2,PP)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS FUNCTION CALCULATES
C  3J( J1 P1 J1P)*3J(JJ PP  JJP)*3J(J2 P2 J2P)*9J(JJ PP JJP)
C    (-K1 Q1 K1P)   (MU  0 -MU )   ( 0  0   0)   (J1 P1 J1P)
C                                                (J2 P2 J2P)
C  *ROOT((2J1+1)(2J1P+1)(2J2+1)(2J2P+1)(2PP+1)(2P2+1)(2JJ+1)(2JJP+1))
C  *(-1)^(J1P+J2P+JJ+MU-K1)/4PI
C
C  ON ENTRY: ALL VARIABLES ARE INTEGER BASIS FUNCTION LABELS.
C            ALL UNCHANGED ON EXIT.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER P1,Q1,P2,PP
      DATA Z0 /0.D0/, PI /3.14159265358979289D0/
C  INTERNAL FUNCTION...
      Z(X)=2.D0*X+1.D0
C
      XJ1 = J1
      XK1 = K1
      XJ2 = J2
      XJ1P = J1P
      XK1P = K1P
      XJ2P = J2P
      XJJ = JJ
      XJJP = JJP
      XMU = MU
      XQ1 = Q1
      XP1 = P1
      XP2 = P2
      XPP = PP
      RSYMTP=0.D0
      F=THRJ(XJ1,XP1,XJ1P,-XK1,XQ1,XK1P)
      IF (ABS(F).LE.1.D-8) RETURN
      F=F*THRJ(XJJ,XPP,XJJP,XMU,Z0,-XMU)
      IF (ABS(F).LE.1.D-8) RETURN
      F = F*THREEJ(J2,P2,J2P)
      IF (ABS(F).LE.1.D-8) RETURN
      F = F*XNINEJ(JJ,PP,JJP,J1,P1,J1P,J2,P2,J2P)
      IF (ABS(F).LE.1.D-8) RETURN
      RSYMTP=F*SQRT(Z(XJ1)*Z(XJ1P)*Z(XJ2)*Z(XJ2P)*Z(XPP)*Z(XP2)
     1       *Z(XJJ)*Z(XJJP))*PARSGN(J1P+J2P+JJ+MU-K1)/(4.0D0*PI)
      RETURN
      END
