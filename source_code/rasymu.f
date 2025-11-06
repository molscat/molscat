      SUBROUTINE RASYMU(NSTATE,IASYMU,EASYM,EFACT,EUNAME,IPRINT,LSET6)
C  Copyright (C) 2025 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  WRITTEN BY CRLS 08-2022

C  COMMON CODE EXTRACTED FROM SET4/SET6 WITH ADDITIONAL REFINEMENTS
C
      USE basis_data, ONLY: ELEVEL, EMAX, JLEVEL, JMAX, JMIN,
     b                      JSTEP, NLEVEL, MXELVL, ROTI
      USE pair_state, ONLY: ATAU, JSTATE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION EASYM(*)
      CHARACTER(8) EUNAME
      INTEGER IASYMU
      LOGICAL LSET6

      ALLOCATABLE :: JTEMP(:),ATEMP(:),AT2(:)

      CHARACTER(1) PLUR(2)
      LOGICAL LEVIN,LKEEP,EIN,LIN

      DATA PLUR /' ','S'/

      NJTERM=3
      IF (LSET6) NJTERM=2
      LIN=NLEVEL.LT.0
      LEVIN=NLEVEL.GT.0
      EIN=.FALSE.

      IF (IASYMU.NE.5) THEN
        OPEN (IASYMU,STATUS='OLD',ACTION='READ')
        REWIND (IASYMU)
      ENDIF

      IF (IPRINT.GE.1) WRITE(6,602) IASYMU,TRIM(EUNAME),EFACT
  602 FORMAT(/'  ASYMMETRIC TOP BASIS WILL BE INPUT FROM UNIT ',
     1       'IASYMU =',I4//'  QUANTITIES IN IASYMU MUST BE IN UNITS',
     2       ' SPECIFIED BY ',A,' = ',G12.6,' CM-1'/)

!  Start of long IF block #1
      IF (NLEVEL.GT.0) THEN
        NREAD=NLEVEL
        IF (IPRINT.GE.1) WRITE(6,603) NLEVEL,PLUR(MIN(NLEVEL,2))
  603   FORMAT(' ',10X,I6,'  INPUT LEVEL',A1,' SPECIFIED BY NLEVEL.')
      ELSE
        IF (IASYMU.EQ.5) THEN
          WRITE(6,*) '  *** RASYMU. CANNOT READ FROM STD INPUT FOR',
     1               ' NLEVEL.LE.0'
          STOP
        ENDIF
        NREAD=1000000
        IF (LIN) THEN
          NLEVEL=ABS(NLEVEL)
          EINMAX=0.D0
          EIN=LSET6
          IF (EIN) THEN
            DO I=1,NLEVEL
            EINMAX=MAX(EINMAX,ABS(ELEVEL(I)))
            ENDDO
            EIN=EIN .AND. EINMAX.GT.0.D0
          ENDIF
          IF (IPRINT.GE.1) THEN
            WRITE(6,613) -ABS(NLEVEL)
  613       FORMAT(2X,'NEGATIVE NLEVEL =',I5,
     1             ' WILL SCREEN INPUT ON &BASIS JLEVEL()')
            IF (EIN) THEN
              WRITE(6,*) '      ENERGIES TAKEN FROM &BASIS ELEVEL'
            ELSEIF (LSET6) THEN
              WRITE(6,*) '      ENERGIES TAKEN FROM IASYMU'
            ENDIF
          ENDIF
        ENDIF
      ENDIF
!  End of long IF block #1
C
      NSTATE=0
      NKVAL=0
!  Start of long DO loop #1
      DO III=1,NREAD
        READ(IASYMU,*,END=9000) JI,ITAU,EINP
        IF (NSTATE.GT.MXELVL) THEN
          WRITE(6,*) '  *** RASYMU. DIMENSION OF ELEVEL EXCEEDED',NSTATE
          STOP
        ENDIF
        JI=ABS(JI)
        NK=2*JI+1
        IF (ALLOCATED(AT2)) DEALLOCATE (AT2)
        ALLOCATE (AT2(NK))
        READ(IASYMU,*,END=9100) AT2
!  Start of long IF block #2
        IF (LIN) THEN
C  CODE BELOW FILTERS IASYMU INPUT ON JLEVEL
          LKEEP=.FALSE.
          DO IND=1,NLEVEL
            IF (JLEVEL(NJTERM*(IND-1)+1).EQ.JI .AND.
     1          JLEVEL(NJTERM*(IND-1)+2).EQ.ITAU) THEN
              LKEEP=.TRUE.
              INDX=IND
              IF (EIN .AND. ELEVEL(INDX).NE.0.D0) EINP=ELEVEL(IND)/EFACT
            ENDIF
            IF (LKEEP) EXIT
          ENDDO
          IF (IPRINT.GE.1 .AND. .NOT.LKEEP)
     1      WRITE(6,683) JI,ITAU,EINP,TRIM(EUNAME)
  683       FORMAT(/'  SKIP INPUT LEVEL      J, TAU =',2I5,
     1             '  ENERGY =',F13.3,1X,A,'.  NOT IN JLEVEL LIST')
          IF (.NOT.LKEEP) CYCLE
        ELSE
          LKEEP=.TRUE.
          IF (EMAX.GT.0.D0) THEN
            IF (EINP.GT.EMAX) LKEEP=.FALSE.
            IF (IPRINT.GE.1 .AND. .NOT.LKEEP)
     1        WRITE(6,680) JI,ITAU,EINP,TRIM(EUNAME),EMAX,TRIM(EUNAME)
  680         FORMAT(/'  SKIP INPUT LEVEL      J, TAU =',2I5,
     1               '  ENERGY = ',F13.3,1X,A,' ABOVE EMAX =',
     2               F11.3,1X,A)
            IF (.NOT.LKEEP) CYCLE
          ENDIF
          IF (JMAX.GT.0) THEN
            IF (JI.LT.JMIN .OR. JI.GT.JMAX) THEN
              LKEEP=.FALSE.
              IF (IPRINT.GE.1 .AND. .NOT.LKEEP)
     1          WRITE(6,610) JI,ITAU,EINP,TRIM(EUNAME),JMIN,JMAX
  610           FORMAT(/'  SKIP INPUT LEVEL      J, TAU =',2I5,
     1                 '  ENERGY = ',F13.3,1X,A,'.  J NOT IN RANGE',
     2                 I4,' TO',I4)
              IF (.NOT.LKEEP) CYCLE
            ENDIF
          ENDIF

C  GET PARITY CODE FROM ATAU SYMMETRIES. . .
          IPAR=IPASYM(JI,NK,AT2)
C  IPAR=-1 IS ERROR RETURN FROM IPASYM.
          IF (IPAR.EQ.-1) THEN
            WRITE(6,699)
            STOP
          ENDIF

C  SEE IF WE SHOULD SKIP ON JMIN,JMAX,JSTEP OR EMAX
          IF (JMAX.GT.0) THEN
C  PRE-2023 CODE FOR JSTEP SELECTION ONLY ON J
C           JDIF=JI-JMIN
C           JDIF=MOD(JDIF,JSTEP)
C           IF (JDIF.NE.0 .OR. JI.LT.JMIN .OR. JI.GT.JMAX) THEN
C  NEW CODE TO SELECT ON (-1)**(J+K+IPAR)
            KK=IPAR/2
            IF (MOD(JI+KK+IPAR,JSTEP).NE.MOD(JMIN,JSTEP)) THEN
              LKEEP=.FALSE.
              IF (IPRINT.GE.1 .AND. .NOT.LKEEP)
     1          WRITE(6,611) JI,ITAU,EINP,TRIM(EUNAME),JSTEP
  611           FORMAT(/'  SKIP INPUT LEVEL      J, TAU =',2I5,
     1                 '  ENERGY = ',F13.3,1X,A,
     2                 '. J+K+PRTY NOT AS REQUIRED FOR ',
     3                 'JSTEP =',I4)
              IF (.NOT.LKEEP) CYCLE
            ENDIF
          ENDIF
        ENDIF
!  End of long IF block #2
        ETEMP=EINP*EFACT
C
C  REACH BELOW IF WE ARE INCLUDING THIS SET
C  CRLS 06-2022: EVERY TIME A NEW STATE IS INCLUDED, EXPAND THE SIZES OF
C                JSTATE AND ATAU BY COPYING THEM TO TEMPORARY ARRAYS,
C                REALLOCATING THEM AS LARGER ARRAYS AND COPYING THE CONTENTS
C                OF THE TEMPORARY ARRAYS BACK
        NSTATE=NSTATE+1
        IF (LSET6) THEN
          ELEVEL(NSTATE)=EINP*EFACT
        ELSE
          EASYM(NSTATE)=EINP*EFACT
        ENDIF

        IF (NSTATE.GT.1) THEN
          ALLOCATE (JTEMP(6*(NSTATE-1)),ATEMP(NKVAL))
          JTEMP=JSTATE
          ATEMP=ATAU
          DEALLOCATE (JSTATE,ATAU)
        ENDIF
        ALLOCATE (JSTATE(6*NSTATE),ATAU(NKVAL+NK))
        IF (NSTATE.GT.1) THEN
          IPOS=0
          DO ILAB=1,6
          DO ISTATE=1,NSTATE-1
            IPOS=IPOS+1
            INEW=(ILAB-1)*NSTATE+ISTATE
            JSTATE(INEW)=JTEMP(IPOS)
          ENDDO
          ENDDO
          ATAU(1:NKVAL)=ATEMP
          DEALLOCATE (JTEMP,ATEMP)
        ENDIF
        DO I=1,NK
          ATAU(NKVAL+I)=AT2(I)
        ENDDO
C  OUTPUT INFORMATION READ.
        IF (IPRINT.GE.1) WRITE(6,604) NSTATE,JI,ITAU,EINP,TRIM(EUNAME),
     1                                EINP*EFACT
  604   FORMAT(/'  KEEP INPUT LEVEL',I4,'  J, TAU =',2I5,
     1         '  ENERGY =',F15.5,1X,A,'  =  ',F15.5,' CM-1')
        MJI=-JI
        IF (IPRINT.GE.1) WRITE(6,605) (ATAU(NKVAL+I+JI+1),I, I=MJI,JI)
  605   FORMAT(10X,'INPUT COEFFICIENTS ARE (K VALUES IN PARENTHESES)'/
     1         (10X,6(F12.6,'(',I3,')')))
C
C  GET PARITY CODE FROM ATAU SYMMETRIES. . .
        IPAR=IPASYM(JI,NK,ATAU(NKVAL+1))
C  IPAR=-1 IS ERROR RETURN FROM IPASYM.
        IF (IPAR.EQ.-1) THEN
          WRITE(6,699)
          STOP
        ENDIF

C  REORDER JSTATE TO RECEIVE NEW ROW.
C
        JSTATE(  NSTATE)=JI
        JSTATE(2*NSTATE)=ITAU
        JSTATE(3*NSTATE)=IPAR
        JSTATE(4*NSTATE)=NKVAL
        JSTATE(5*NSTATE)=NK
        NKVAL=NKVAL+NK
        IF (LIN .AND. LSET6) THEN
          JSTATE(6*NSTATE)=INDX
        ELSE
          JSTATE(6*NSTATE)=NSTATE
        ENDIF
        CYCLE
C
C  * * * END OF FILE CONDITIONS * * *
 9000   IF (.NOT.LEVIN) THEN

          WRITE(6,606) IASYMU,NSTATE
  606     FORMAT(/'  END OF FILE ENCOUNTERED ON UNIT',I4,'.   ',I5,
     &           '  ROTOR FUNCTIONS KEPT.')
          EXIT
        ENDIF

c  somETHING WRONG ABOUT THIS...
        WRITE(6,607) IASYMU,NSTATE
  607   FORMAT(/'  PREMATURE E.O.F. ON UNIT',I4,
     &         '.  NLEVEL REDUCED TO',I6)
        NLEVEL=NSTATE
        EXIT

 9100   WRITE(6,608) IASYMU,NSTATE
  608   FORMAT(/'  * * * ERROR.   E.O.F. ON UNIT',I4,
     &         '  BEFORE ATAU CARDS FOR NSTATE =',I5)
        WRITE(6,699)
  699   FORMAT(/'  * * * TERMINAL ERROR.')
        STOP
      ENDDO
!  End of long DO loop #1
C  THIS COMPLETES READ(IASYMU) LOOP
C
c only IN SET6
      IF (IASYMU.NE.5) CLOSE(IASYMU)
      IF (NLEVEL.EQ.0) NLEVEL=NSTATE
 2400 IF (LSET6 .AND. LIN) THEN
C  WE FILTERED ON JLEVEL(), MAKE SURE WE HAVE THEM ALL
        IF (NSTATE.NE.NLEVEL) THEN
          WRITE(6,*) ' ALL LEVELS SPECIFIED BY JLEVEL() WERE NOT FOUND'
          WRITE(6,*) ' *** TERMINAL ERROR.'
          STOP
        ENDIF
      ENDIF
C  MAKE SURE EACH VALUE IS THERE AND REORDER IF NECESSARY
C  SO THAT JSTATE(I,6)=I (EXPECTED BY PRBR, EG)
!  Start of long IF block #3
      IF (LIN) THEN
        IF (.NOT.LSET6) THEN
          JMIN=JSTATE(1)
          JMAX=JMIN
          DO I=1,NSTATE
            JMIN=MIN(JMIN,JSTATE(I))
            JMAX=MAX(JMAX,JSTATE(I))
          ENDDO
        ELSE
        LOOP_I: DO I=1,NLEVEL
          DO IX=1,NSTATE
            IF (I.NE.JSTATE(5*NSTATE+IX)) CYCLE
            IF (I.EQ.IX) CYCLE LOOP_I

            IF (LSET6) THEN
              DO IC=1,6
                ITMP=JSTATE((IC-1)*NSTATE+I)
                JSTATE((IC-1)*NSTATE+I)=JSTATE((IC-1)*NSTATE+IX)
                JSTATE((IC-1)*NSTATE+IX)=ITMP
              ENDDO
            ENDIF
            CYCLE LOOP_I

          ENDDO
          WRITE(6,684) I,JLEVEL(NJTERM*(I-1)+1),JLEVEL(NJTERM*(I-1)+2)
  684     FORMAT('  INPUT SET',I4,'  J, TAU =',2I5,
     &           '  NOT FOUND ON IASYMU'/'  *** TERMINAL ERROR')
          STOP
        ENDDO LOOP_I
        ENDIF
      ELSE
C  SET J,TAU INTO JLEVEL; GET JMIN,JMAX
        NLEVEL=NSTATE
C  SET JLEVEL(), JMIN, AND JMAX.
        IF (LSET6) THEN
          JMIN=JSTATE(1)
          JMAX=JMIN
        ENDIF
        DO I=1,NSTATE
          IF (LSET6) THEN
            JI=JSTATE(I)
            JMIN=MIN(JMIN,JI)
            JMAX=MAX(JMAX,JI)
          ENDIF
          JLEVEL(2*I-1)=JI
          JLEVEL(2*I)=JSTATE(I+NSTATE)
        ENDDO
      ENDIF
!  End of long IF block #3
C  CHECK THAT FUNCTIONS ARE ORTHOGONAL
      CALL CHECK6(NSTATE,JSTATE,ATAU,.FALSE.)

      RETURN
      END
