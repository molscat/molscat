      DOUBLE PRECISION FUNCTION CALCK(EP,ERED,EINT,N)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3

C  CR Le Sueur Dec 2018
      IMPLICIT NONE

C  THIS ROUTINE CALCULATES THE HIGHEST K VECTOR (WITH AN OPTIONAL OFFSET EP)

      INTEGER,          INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(IN) :: EP,ERED,EINT(N)

      INTEGER I

      CALCK=0.D0
      DO I=1,N
        IF (EP+ERED-EINT(I).LE.0.D0) CYCLE
        CALCK=MAX(CALCK,SQRT(EP+ERED-EINT(I)))
      ENDDO

      RETURN
      END
