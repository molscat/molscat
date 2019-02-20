      SUBROUTINE CHNSRT(NB,EINT,N)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  Extracted from other routines by CR Le Sueur 2015 ish
C  THIS ROUTINE ORDERS THE INTERNAL ENERGIES FROM LOWEST TO HIGHEST.
C  THE ORDERING IS RETURNED IN NB
      IMPLICIT NONE

      INTEGER, INTENT(IN)::N
      INTEGER, INTENT(OUT)::NB(N)
      INTEGER I,J,IT

      DOUBLE PRECISION, INTENT(IN):: EINT(N)

      DO I=1,N
        NB(I)=I
      ENDDO

      IF (N.EQ.1) RETURN

      DO I=1,N-1
        DO J=I+1,N
          IF (EINT(NB(I)).LE.EINT(NB(J))) CYCLE
          IF (ABS(EINT(NB(I))-EINT(NB(J))).LT.1D-11) CYCLE
          IT=NB(I)
          NB(I)=NB(J)
          NB(J)=IT
        ENDDO
      ENDDO

      RETURN
      END
