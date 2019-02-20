      FUNCTION PARSGN(I)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARSGN=1.D0
      IF ((I/2)*2-I.NE.0) PARSGN=-1.D0
      RETURN
      END
