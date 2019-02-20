      SUBROUTINE VRTP(IDERIV,RM,P)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION P(*)

      IF (IDERIV.NE.-1) THEN
        WRITE(6,601) IDERIV
  601   FORMAT(/' *** ERROR.  DUMMY VERSION OF VRTP CALLED WITH ',
     1         'IDERIV =',I4/
     2         14X,'VRTP MUST BE PROVIDED IF NTERM(I) IS ZERO.')
        STOP
      ENDIF
      END
