      SUBROUTINE VINIT(I,RM,EPSIL)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      ENTRY VSTAR (I,RR,SUM)
      ENTRY VSTAR1(I,RR,SUM)
      ENTRY VSTAR2(I,RR,SUM)
      WRITE(6,601) I
  601 FORMAT(/' *** ERROR.  DUMMY VERSION OF VINIT CALLED WITH I =',
     1       I4/14X,'VINIT MUST BE PROVIDED IF NTERM(I) IS ZERO.')
      STOP
      END
