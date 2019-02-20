      module angles
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE sizes, ONLY: MXANG
      IMPLICIT NONE

C  DATA FOR ANGLES FOR QUADRATURE ETC
      INTEGER          :: IHOMO, ICNSYM, IHOMO2, ICNSY2

      DOUBLE PRECISION :: COSANG(MXANG), FACTOR
C  REPLACES THIS COMMON BLOCK
C     COMMON /ANGLES/ COSANG,FACTOR,IHOMO,ICNSYM,IHOMO2,ICNSY2

      end module angles
