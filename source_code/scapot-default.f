      SUBROUTINE SCAPOT(P,MXLAM)
      USE efvs, ONLY: SCALAM
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  WRITTEN ON 26-10-2018 BY CRLS
C
C  THIS SUBROUTINE SCALES ALL THE POTENTIAL COEFFICIENTS P BY THE
C  SCALING FACTOR
C
      IMPLICIT NONE

      INTEGER,          INTENT(IN)    :: MXLAM
      DOUBLE PRECISION, INTENT(INOUT) :: P(MXLAM)

      INTEGER I

      DO I=1,MXLAM
        P(I)=P(I)*SCALAM
      ENDDO

      RETURN
      END
