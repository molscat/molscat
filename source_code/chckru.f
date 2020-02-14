      SUBROUTINE CHCKRU(RUNIT,RUNAME,RPUNIT,RSCALE,unset,IPRINT)
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE potential, ONLY: RMNAME
      USE physical_constants, ONLY: bohr_in_SI
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: RUNIT,RPUNIT,RSCALE
      DOUBLE PRECISION, INTENT(IN)    :: unset
      INTEGER,          INTENT(IN)    :: IPRINT
      CHARACTER(10),    INTENT(INOUT) :: RUNAME

      IF (RUNIT.EQ.unset) THEN
        IF (RPUNIT.EQ.unset) THEN
          RUNIT=1.D0
          RPUNIT=1.D0
        ELSE
          RUNIT=RPUNIT
          IF (RUNAME.EQ.'RUNIT') RUNAME=RMNAME
        ENDIF
      ELSE
        IF (RPUNIT.EQ.unset) THEN
          RPUNIT=RUNIT
          IF (RUNAME.EQ.'RUNIT') RUNAME=RMNAME
        ELSE
           IF (RUNIT.EQ.RPUNIT) THEN
            IF (RUNAME.EQ.'RUNIT') RUNAME=RMNAME
          ELSEIF (IPRINT.GE.6) THEN
            WRITE(6,100)
  100 FORMAT(/'  NOTE THAT POTENTIAL ROUTINE OPERATES IN DIFFERENT ',
     1        'LENGTH UNITS FROM REST OF PROGRAM')
          ENDIF
        ENDIF
      ENDIF
      RSCALE=RUNIT/RPUNIT
C
      IF (RUNIT.EQ.1.D0) THEN
        RUNAME='ANGSTROM'
      ELSEIF (ABS(RUNIT-bohr_in_SI*1.D10).LT.1.D-6) THEN
        RUNAME='BOHR    '
      ENDIF
C
      IF (IPRINT.GE.1) THEN
        IF (RUNIT.EQ.1.D0) THEN
          WRITE(6,200)
  200     FORMAT(/'  ALL LENGTHS ARE IN UNITS OF ANGSTROM UNLESS ',
     1            'OTHERWISE STATED'/)
        ELSE
          WRITE(6,210) TRIM(RUNAME),RUNIT
  210     FORMAT(/'  ALL LENGTHS ARE IN UNITS OF ',A,' (',F12.8,
     1            ' ANGSTROM ) UNLESS OTHERWISE STATED'/)
        ENDIF
      ENDIF

      RETURN
      END
