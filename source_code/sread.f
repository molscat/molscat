      SUBROUTINE SREAD(IU,N,S,IEND)
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  ON ENTRY: IU IS NUMBER LABEL OF (UNFORMATTED) INPUT STREAM;
C            S IS ARRAY OF SIZE (N,N).
C  ON EXIT:  IEND=0 INDICATES SUCCESSFUL READ;
C                 1 INDICATES FAILURE.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(N,N)
C
      IEND=0
      READ(IU,END=9999) ((S(I,J),J=1,I),I=1,N)
      DO 1000 I=1,N
      DO 1000 J=1,I-1
 1000 S(J,I)=S(I,J)
      RETURN
 9999 IEND=1
      RETURN
      END
C=======================================================================
      SUBROUTINE SKREAD
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  WRITTEN BY CR Le Sueur Aug 2018
      IMPLICIT NONE
      SAVE
C  THE FOLLOWING VARIABLES ARE IN THE MODULE EFVS.
      INTEGER, PARAMETER            :: MXEFV=10,LEFVN=20,LEFVU=6
      INTEGER, INTENT(OUT)          :: NEFV, ISVEFV, ITPSUB

      CHARACTER(LEFVN), INTENT(OUT) :: EFVNAM(0:MXEFV)
      CHARACTER(LEFVU), INTENT(OUT) :: EFVUNT(0:MXEFV)
C
C  PASSED VARIABLES FOR HEADER
      INTEGER, INTENT(IN)           :: IRDUNT,NSTATE,NQN,IPROGM
      INTEGER, INTENT(OUT)          :: JSTATE(NSTATE*NQN),NLVL,NNRG,
     1                                 NFIELD,NDGVL,NCONST,NRSQ,IBOUND

      DOUBLE PRECISION, INTENT(OUT) :: ELEVEL(1),ENERGY(1)

      LOGICAL, INTENT(IN)           :: LFIN
C  PASSED VARIABLES FOR LOOPS
      INTEGER, INTENT(IN)           :: IREAD
      INTEGER, INTENT(OUT)          :: JTOT,INRG,IEXCH,M,NOPEN,IFIELD,
     1                                 L(1),INDLEV(1)

      DOUBLE PRECISION, INTENT(OUT) :: ENERGN,WT,WVEC(1),EREF,CENT(1),
     1                                 SREAL(1),SIMAG(1),AKMAT(1)
      DOUBLE PRECISION, INTENT(OUT) :: EFV(0:MXEFV)

      INTEGER, INTENT(OUT)          :: IFAIL
C
C  LOCAL VARIABLES
      INTEGER       :: NEFVP, I, IEFV, IEND
      CHARACTER(80) :: LABEL
      LOGICAL       :: KMAT, LFMT
C=======================================================================
C  THIS LINE MUST BE READ BEFORE CALLING THIS SUBROUTINE.
C  THIS READ STATEMENT MUST NOT CHANGE IN THE POST-PROCESSORS!
C
C     READ(IRDUNT) LABEL,ITYPE,NSTATE,NQN,URED,IPROGM
C=======================================================================
C  ENTRY POINT FOR READING HEADER FOR S MATRIX FILE (SBE, DCS, RESTART)
C                                  OR K MATRIX FILE (RESFIT/SAVER)
      ENTRY HDREAD(IRDUNT,LFIN,NSTATE,NQN,IPROGM,
     1            JSTATE,NLVL,ELEVEL,NDGVL,NCONST,NRSQ,IBOUND,
     2            ITPSUB,NEFV,ISVEFV,EFVNAM,EFVUNT,NFIELD,NNRG,ENERGY)
C=======================================================================
      LFMT=LFIN
      IF (IPROGM.GE.17) LFMT=.FALSE.

      IF (LFMT) THEN
        WRITE(6,*) ' *** ERROR: ATTEMPT TO READ FORMATTED '//
     1             'S-MATRIX FILE'
        WRITE(6,*) '     CODE HAS BEEN COMMENTED OUT'
        STOP
C       READ(IRDUNT,101) (JSTATE(I),I=1,NSTATE*NQN)
C       IF (IPROGM.GE.3) READ(IRDUNT,102) NLVL,(ELEVEL(I),I=1,NLVL)
      ELSE

        READ(IRDUNT) (JSTATE(I),I=1,NSTATE*NQN)
        IF (IPROGM.GE.3) READ(IRDUNT) NLVL,(ELEVEL(I),I=1,NLVL)
        IF (IPROGM.GE.17) THEN
          READ(IRDUNT) NDGVL,NCONST,NRSQ,IBOUND,ITPSUB
          READ(IRDUNT) NEFV,ISVEFV,(EFVNAM(IEFV),EFVUNT(IEFV),
     2                              IEFV=1,NEFV)
          NEFVP=0
          IF (NEFV.GT.0) NEFVP=MAX(NEFV,ISVEFV)
        ENDIF
      ENDIF

      IF (LFMT) THEN
C       READ(IRDUNT,102) NNRG,(ENERGY(I),I=1,NNRG)
      ELSE
        IF (IPROGM.LT.17) THEN
          NFIELD=1
          READ(IRDUNT) NNRG,(ENERGY(I),I=1,NNRG)
        ELSE
          READ(IRDUNT) NFIELD,NNRG,(ENERGY(I),I=1,NNRG)
        ENDIF
      ENDIF

      RETURN

C  FROM SBE AND DCS
  101 FORMAT(20I4)
  102 FORMAT(I4/(5E16.8))
C=======================================================================
C  ENTRY POINT FOR READING K MATRIX FILE (RESFIT/SAVER)
      ENTRY KLPRD(IRDUNT,JTOT,INRG,ENERGN,EFV,EREF,IEXCH,WT,M,NOPEN,
     1            IFIELD,INDLEV,L,CENT,WVEC,SREAL,SIMAG,AKMAT,
     2            IPROGM,IBOUND,ITPSUB,IREAD,IFAIL)
      KMAT=.TRUE.
      GOTO 20
C=======================================================================
C  ENTRY POINT FOR READING S MATRIX FILE (SBE, DCS, RESTART)
      ENTRY SLPRD(IRDUNT,JTOT,INRG,ENERGN,EFV,EREF,IEXCH,WT,M,NOPEN,
     1            IFIELD,INDLEV,L,CENT,WVEC,SREAL,SIMAG,
     2            IPROGM,IBOUND,ITPSUB,IREAD,IFAIL)
      KMAT=.FALSE.
C=======================================================================

C  THIS READS THE LOOP INFORMATION REQUIRED BY THE POST-PROCESSING
C  PROGRAMS.  IT IS FULL OF READ STATEMENTS FOR LEGACY OUTPUTS.
   20 IFAIL=0
C  FOR DCS AND SBE CODES
      IF (LFMT) THEN
C       READ(IRDUNT,103,END=1000) JTOT,INRG,ENERGN,IEXCH,WT,M
C       READ(IRDUNT,104,END=2000) NOPEN,(JSINDX(I),L(I),WVEC(I),
C    1                                   I=1,NOPEN)
      ELSE
        IF (IPROGM.LT.14) THEN
          READ(IRDUNT,END=1000) JTOT,INRG,ENERGN,IEXCH,WT,M
          READ(IRDUNT,END=2000) NOPEN,(INDLEV(I),L(I),WVEC(I),
     1                                 I=1,NOPEN)
        ELSEIF (IPROGM.LT.17) THEN
          READ(IRDUNT,END=1000) JTOT,INRG,ENERGN,EFV(1),
     1                          EFV(2),IEXCH,WT,M,NOPEN
          READ(IRDUNT,END=2000) (INDLEV(I),L(I),WVEC(I),I=1,NOPEN)
        ELSE
C  BELOW IS THE STRUCTURE FOR OUTPUT GENERATED FROM THE CURRENT VERSION OF
C  MOLSCAT (IPROGM=17).
C
C  IREAD CAN BE USED TO READ JUST ONE LINE AT A TIME, OR, IF SET TO 0
C  THIS SUBROUTINE WILL READ ALL THE LINES FOR THE CURRENT CYCLE.
          IF (IREAD.EQ.1 .OR. IREAD.EQ.0) THEN
            READ(IRDUNT,END=1000) JTOT,INRG,M,IFIELD,ENERGN,
     1                            (EFV(IEFV),IEFV=0,NEFVP),
     2                            EREF,IEXCH,WT,NOPEN
          ENDIF
C  NOTE THAT FOR DIAGONAL BASES, THE QUANTITY REFERRED TO AS INDLEV HERE
C  IS ACTUALLY JSINDX
          IF (IREAD.EQ.2 .OR. IREAD.EQ.0) THEN
            IF (IBOUND.EQ.0) THEN
              READ(IRDUNT,END=2000) (INDLEV(I),L(I),WVEC(I),I=1,NOPEN)
            ELSE
              READ(IRDUNT,END=2000) (INDLEV(I),CENT(I),WVEC(I),I=1,
     1                                                          NOPEN)
            ENDIF
          ENDIF
          IEND=0
          IF (IREAD.EQ.3 .OR. IREAD.EQ.0) THEN
            IF (KMAT) THEN
              READ(IRDUNT,END=3000) (SREAL(I),SIMAG(I),I=1,NOPEN)
            ELSE
              CALL SREAD(IRDUNT,NOPEN,SREAL,IEND)
              IF (IEND.GT.0) GOTO 3000
            ENDIF
          ENDIF
          IF (IREAD.EQ.4 .OR. IREAD.EQ.0) THEN
            IF (KMAT) THEN
              CALL SREAD(IRDUNT,NOPEN,AKMAT,IEND)
              IF (IEND.GT.0) GOTO 4000
            ELSE
              CALL SREAD(IRDUNT,NOPEN,SIMAG,IEND)
              IF (IEND.GT.0) GOTO 4000
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      RETURN

 1000 IFAIL=1
      RETURN

 2000 IFAIL=2
      RETURN

 3000 IFAIL=3
      RETURN

 4000 IFAIL=4
      RETURN
C  FROM SBE AND DCS
  103 FORMAT(2I4,E16.8,I4,E16.8,I4)
  104 FORMAT(I4/(2I4,E16.8))

      END
C=======================================================================
