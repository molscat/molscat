      FUNCTION PASLEG(L,MM,X)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  CALCULATE NORMALISED ASSOCIATED LEGENDRE POLYNOMIALS
C
      DIMENSION P(450)

      NREQ=(L+1)*(L+2)/2

      IF (NREQ.GT.450) THEN
        WRITE(6,601) L
  601   FORMAT('  *** ERROR IN PASLEG - NOT ENOUGH STORAGE FOR L=',I3)
        STOP
      ENDIF

      CALL ASSLEG(P,L,X,450)

      M=ABS(MM)
      IND=L*(L+1)/2+M+1
      FAC=0.5D0*DBLE(L+L+1)

      DO 100 I=L-M+1,L+M
        FAC=FAC/DBLE(I)
  100 CONTINUE

      PASLEG=P(IND)*SQRT(FAC)

      RETURN
      END
