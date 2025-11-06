      FUNCTION PLEG(N,X)
C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  FUNCTION TO GENERATE LEGENDRE POLYNOMIALS
C
      PLEG=1.D0
      IF (N.EQ.0) RETURN

      P0=1.D0
      P1=X

      DO K=3,N+1
        TEMP=(DBLE(2*K-3)*X*P1 - DBLE(K-2)*P0) / DBLE(K-1)
        P0=P1
        P1=TEMP
      ENDDO

      PLEG=P1

      RETURN
      END
