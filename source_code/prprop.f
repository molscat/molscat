C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  This set of subroutines used to print out information relating to
C  the progress of a bound-state location
C================================================================ PRPROP
      SUBROUTINE PRPROP(SVNAME,SVVAL,SVUNIT)
C  PRINT ENERGY OR EFV FOR CURRENT PROPAGATION
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: SVVAL
      CHARACTER(*),     INTENT(IN) :: SVNAME,SVUNIT
      CHARACTER(38) SVTEMP !73-35

      SVTEMP=SVNAME
      SVTEMP=ADJUSTR(SVTEMP)

      WRITE(6,100) SVTEMP,SVVAL,SVUNIT
  100 FORMAT(/2X,'COUPLED EQUATIONS PROPAGATED AT ',A,' = ',1PG17.10, !35+A1
     1       1X,A)

      RETURN
      END
C================================================================ PREREF
      SUBROUTINE PREREF(EREF,EFACT,EUNIT)
C  PRINT REFERENCE ENERGY
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: EREF,EFACT
      CHARACTER(*),     INTENT(IN) :: EUNIT

      IF (EREF.EQ.0.D0) RETURN
      WRITE(6,200) EREF,EREF/EFACT,EUNIT
  200 FORMAT(2X,1P,'REFERENCE ENERGY IS ',25X,G17.10,1X,'CM-1',3X,
     1             ' = ',G17.10,1X,A)
      RETURN
      END
C================================================================ PREABS
      SUBROUTINE PREABS(E,EREF,EFACT,EUNIT)
C  PRINT ABSOLUTE ENERGY
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: E,EREF,EFACT
      CHARACTER(*),     INTENT(IN) :: EUNIT

      IF (EREF.EQ.0.D0) RETURN
      IF (EFACT.EQ.1.D0) THEN
        WRITE(6,100) E
      ELSE
        WRITE(6,200) E,E/EFACT,EUNIT
      ENDIF
  100 FORMAT(2X,1P,'ABSOLUTE ENERGY  IS ',53X,G17.10,1X,'CM-1')
  200 FORMAT(2X,1P,'ABSOLUTE ENERGY  IS ',25X,G17.10,1X,'CM-1',3X,
     1             ' = ',G17.10,1X,A)
      RETURN
      END
