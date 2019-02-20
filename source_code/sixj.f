      FUNCTION SIXJ(J1,J2,J5,J4,J3,J6)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  CALCULATES 6-J SYMBOL:   _(J1 J2 J3 )_
C                            (J4 J5 J6 )
C  INTERFACE TO J6J ROUTINE.
C  MODIFIED BY S. GREEN 20 AUG 93; PASS DIMENSION OF XJ6J FOR CHECKING
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXDIM=200)
      DIMENSION XJ6J(MXDIM)
      IVAL=MXDIM
      CALL J6J(         DBLE(J2),DBLE(J3),
     1         DBLE(J4),DBLE(J5),DBLE(J6),
     3         IVAL,XJ1MIN,XJ6J)
      IND=1+J1-INT(XJ1MIN+0.1D0)
      SIXJ=0.D0
      IF (IND.GE.1 .AND. IND.LE.IVAL) SIXJ=XJ6J(IND)
      RETURN
      END
