      SUBROUTINE SET6(LEVIN,EIN,NSTATE,JSTATE,ATAU,EFACT,IUNIT,IPRINT)
C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE basis_data, ONLY: ELEVEL, EMAX, JLEVEL, JMAX, JMIN,
     b                      JSTEP, NLEVEL, MXELVL, ROTI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  REVISED FOR VERSION 14:
C  THREE POSSIBLE METHODS OF SPECIFYING ASYMMETRIC TOP LEVELS
C    1. A,B,C .GT.0 IMPLIES GENERATE VIA SET6C
C    2. NLEVEL.GE.0 IMPLIES READ FROM IASYMU (FILTER ON JMIN,JMAX,
C       JSTEP,EMAX); IF (NLEVEL.EQ.0) READ TO END-OF-FILE
C    3. NLEVEL.LT.0 IMPLIES READ FROM IASYMU BUT ACCEPT ONLY THOSE
C       J,ITAU CORRESPONDING TO JLEVEL(2*I-1),JLEVEL(2*I),
C       I=1,ABS(NLEVEL)
C
C  BELOW REPLACES GENERIC SAVE IN V11, WHICH APPEARED UNNECESSARY.
      SAVE IFIRST,NOMEM,NL12,IXMX,ISTART
C
C  THIS ROUTINE HANDLES INPUT, ALSO MATRIX ELEMENTS FOR ITYPE=6.
C  LATTER ARE OBTAINED VIA ENTRIES CPL16, CPL6, CPL26.
C  FIRST VERSION WRITTEN AT MPI, MUNCHEN, JULY 1976.
C  CURRENT VERSION 11 MAR 93 SAVES COUPLING ELEMENTS IN X ARRAY
C    CPL16 (EFF. POTL) COULD BE CHANGED, BUT PROBABLY NO LONGER USED
C
C  N.B. NKVAL HERE COULD BE OBTAINED AS NEEDED - NK=2*J+1.
C       THIS CODE IS MORE FLEXIBLE AS NOT ALL NEED BE STORED BUT
C       K-VALUE COULD BE OBTAINED VIA ADDITIONAL VECTOR KVAL(IST+1).
C
      LOGICAL LEVIN,EIN
      LOGICAL LIN
      LOGICAL NOMEM,LODD
      INTEGER JSTATE(2),IV(*)
      DIMENSION ATAU(2)
C  N.B. JSTATE AND ATAU OCCUPY SAME STORAGE PASSED FROM DRIVER/BASIS.
C       IXNEXT MUST BE INCREMENTED TO REFLECT *ATAU* STORAGE USED
C       (BUT NOT JSTATE, BECAUSE BASIN INCREMENTS IXNEXT BY NQN*NSTATE)
C
C  SPECIFICATIONS FOR CPL16, CPL6, CPL26 ENTRIES.
      INTEGER JSINDX(N),L(N),LAM(2)
      DIMENSION VL(2)
      INTEGER IPRINT
      LOGICAL LFIRST,L20
      CHARACTER(1) PLUR(2)
C
C  DYNAMIC STORAGE COMMON BLOCK ...
      COMMON /MEMORY/ MX,IXNEXT,NIPR,IDUMMY,X(1)
      COMMON /VLSAVE/ IVLU
C
C  DEFAULT INPUT UNIT IS STANDARD INPUT ...
      DATA IDU /5/
      DATA PI /3.14159 26535 89793 D0/
      DATA EPS /1.D-9/, Z0 /0.D0/
      DATA PLUR /' ','S'/
C
C  STATEMENT FUNCTIONS
      F(NN)=DBLE(NN+NN+1)
      LODD(I)=I-2*(I/2).NE.0

      A=ROTI(1)
      B=ROTI(3)
      C=ROTI(5)
C
C  IF ROTATION CONSTANTS ARE INPUT, GENERATE BASIS VIA SET6C
      IF (A.GT.0.D0 .AND. B.GT.0.D0 .AND. C.GT.0.D0) THEN
        CALL SET6C(JSTATE,ATAU,NSTATE,EIN,IPRINT)
C  OPTION ADDED (AUG 94) TO OUTPUT ROTOR WFNS TO IASYMU
        IF (IUNIT.LE.0 .OR. IUNIT.GE.100 .OR. IUNIT.EQ.IDU) RETURN

        IF (IPRINT.GE.1) WRITE(6,1098) IUNIT
 1098   FORMAT(/'  *** SET6 WILL OUTPUT ROTOR WAVEFUNCTIONS TO UNIT',
     1         I4,' FOR FUTURE INPUT'/
     2         '  ENERGIES ARE WRITTEN IN UNITS SPECIFIED BY EUNITS')
        DO 1011 I=1,NSTATE
          JI  =JSTATE(       I)
          ITAU=JSTATE(  NSTATE+I)
          ISTA=JSTATE(3*NSTATE+I)
          NK  =JSTATE(4*NSTATE+I)
          INDX=JSTATE(5*NSTATE+I)
          WRITE(IUNIT,500,ERR=1099) JI,ITAU,ELEVEL(INDX)/EFACT
C         WRITE(IUNIT,*,ERR=1099) JI,ITAU,ELEVEL(INDX)/EFACT
 1011     WRITE(IUNIT,501,ERR=1099) (ATAU(ISTA+II),II=1,NK)
        CLOSE(IUNIT)
        RETURN

 1099   WRITE(6,*) '  *** SET6. ERROR WRITING TO IASYMU; WFNS NOT SAVED'
        RETURN

      ENDIF
C
C  OTHERWISE, INPUT FROM UNIT IASYMU
      IF (IUNIT.GT.0 .AND. IUNIT.LT.100) GOTO 1000

      WRITE(6,601) IUNIT,IDU
  601 FORMAT(/'  ILLEGAL UNIT =',I12,'  SPECIFIED FOR IASYMU, ',
     1       'DEFAULTED TO ',I4)
      IUNIT=IDU
C
 1000 IF (IPRINT.GE.1) WRITE(6,602) IUNIT
  602 FORMAT(/'  ASYMMETRIC TOP BASIS WILL BE INPUT FROM UNIT IASYMU =',
     1       I4//' NB: ENERGIES MUST BE IN UNITS SPECIFIED BY EUNITS'/)
      LIN=.FALSE.

      IF (IUNIT.NE.IDU) REWIND(IUNIT)

      IF (LEVIN) THEN
        NREAD=NLEVEL
        IF (IPRINT.GE.1) WRITE(6,603) NLEVEL,PLUR(MIN(NLEVEL,2))
  603   FORMAT(' ',10X,I6,'  INPUT LEVEL',A1,' SPECIFIED BY NLEVEL.')
      ELSE
        IF (IUNIT.EQ.IDU) THEN
          WRITE(6,*) '  *** SET6. CANNOT READ FROM STD INPUT FOR',
     1               ' NLEVEL.LE.0'
          STOP
        ENDIF
        NREAD=1000000
        IF (NLEVEL.LT.0) THEN
          LIN=.TRUE.
          IF (IPRINT.GE.1) THEN
            WRITE(6,613) NLEVEL
  613       FORMAT(5X,'NEGATIVE NLEVEL =',I5,
     1             ' WILL SCREEN INPUT ON &BASIS JLEVEL()')
            IF (EIN) THEN
              WRITE(6,*) '      ENERGIES TAKEN FROM &BASIS ELEVEL'
            ELSE
              WRITE(6,*) '      ENERGIES TAKEN FROM IASYMU'
            ENDIF
          ENDIF
          NLEVEL=-NLEVEL
        ENDIF
      ENDIF
C
      NSTATE=0
      IOFF=0
      NKVAL=0
      DO 2000 III=1,NREAD
        READ(IUNIT,*,END=9000) JI,ITAU,EINP
c       READ(IUNIT,500,END=9000) JI,ITAU,EINP
  500   FORMAT(2I5,G18.10)
        NSTATE=NSTATE+1
        IF (NSTATE.GT.MXELVL) THEN
          WRITE(6,*) '  *** SET6. DIMENSION OF ELEVEL EXCEEDED',NSTATE
          STOP
        ENDIF
        JI=ABS(JI)
        NK=2*JI+1
        IF (LIN) THEN
C  CODE BELOW FILTERS IASYMU INPUT ON JLEVEL
          DO 2099 IND=1,NLEVEL
            IF (JLEVEL(2*IND-1).NE.JI .OR. JLEVEL(2*IND).NE.ITAU)
     &        GOTO 2099
            INDX=IND
            IF (.NOT.EIN) ELEVEL(INDX)=EINP*EFACT
            GOTO 2090

 2099     CONTINUE
          IF (IPRINT.GE.1) WRITE(6,683) JI,ITAU,EINP
  683     FORMAT(/'  INPUT LEVEL J, TAU, E =',2I5,F13.3,' SKIPPED.',
     1           ' NOT IN JLEVEL LIST')
          GOTO  2070

        ELSE
          ELEVEL(NSTATE)=EINP*EFACT
C  SEE IF WE SHOULD SKIP ON JMIN,JMAX,JSTEP OR EMAX
          IF (JMAX.LE.0) GOTO 2080

          JDIF=JI-JMIN
          JDIF=JDIF-JSTEP*(JDIF/JSTEP)
          IF (JDIF.EQ.0 .AND. JI.GE.JMIN .AND. JI.LE.JMAX) GOTO 2080

          IF (IPRINT.GE.1)
     1      WRITE(6,681) JI,ITAU,ELEVEL(NSTATE),JMIN,JSTEP,JMAX
  681     FORMAT(/'  INPUT LEVEL J, TAU, E =',2I5,F13.3,' SKIPPED. ',
     1           'J NOT IN RANGE',I4,' (',I4,')',I4)
          GOTO 2070

 2080     IF (EMAX.LE.0.D0) GOTO 2090
          IF (ELEVEL(NSTATE).LE.EMAX) GOTO 2090

          IF (IPRINT.GE.1) WRITE(6,680) JI,ITAU,ELEVEL(NSTATE),EMAX
  680     FORMAT(/'  INPUT LEVEL J, TAU, E =',2I5,F13.3,' SKIPPED ',
     1           'DUE TO EMAX =',F11.3)
        ENDIF
C
C  REACH BELOW IF WE ARE SKIPPING THIS SET
 2070   NSTATE=NSTATE-1
        READ( IUNIT,501,END=9100) (ATAUX,I=1,NK)
        GOTO 2000
C
C  REACH BELOW IF WE ARE INCLUDING THIS SET
 2090   CONTINUE
C  SHIFT ATAU BY 6 WORDS TO MAKE ROOM FOR INCOMING JSTATE.
        IOFF=IOFF+6
        DO 2020 I=1,NKVAL
 2020     ATAU(IOFF+NKVAL+1-I)=ATAU(IOFF+NKVAL-5-I)
        INST=IOFF+NKVAL
        READ(IUNIT,501,END=9100) (ATAU(INST+I),I=1,NK)
  501   FORMAT(6F12.8)
C  OUTPUT INFORMATION READ.
        IF (IPRINT.GE.1) WRITE(6,604) NSTATE,JI,ITAU,EINP,ELEVEL(NSTATE)
  604   FORMAT(/'  INPUT LEVEL',I4,'  J, TAU =',2I4,'  ENERGY =',F15.5,
     1         '  =  ',F15.5,' CM-1')
        MJI=-JI
        IF (IPRINT.GE.1) WRITE(6,605) (ATAU(INST+1+JI+I),I, I=MJI,JI)
  605   FORMAT(10X,'INPUT COEFFICIENTS ARE'/(10X,6(F12.6,'(',I3,')')))
C
C  GET PARITY CODE FROM ATAU SYMMETRIES. . .
        IPAR=IPASYM(JI,NK,ATAU(INST+1))
C  IPAR=-1 IS ERROR RETURN FROM IPASYM.
        IF (IPAR.NE.-1) GOTO 2001

        WRITE(6,699)
        STOP

C  REORDER JSTATE TO RECEIVE NEW ROW.
 2001   NRM1=NSTATE-1
        IF (NRM1.LE.0) GOTO 2100

        IOLD=6*NRM1
        IX=6*NSTATE
        DO 2110 II=1,6
          IX=IX-1
          DO 2120  I=1,NRM1
            JSTATE(IX)=JSTATE(IOLD)
            IX=IX-1
 2120       IOLD=IOLD-1
 2110   CONTINUE
C
 2100   JSTATE(  NSTATE)=JI
        JSTATE(2*NSTATE)=ITAU
        JSTATE(3*NSTATE)=IPAR
        JSTATE(4*NSTATE)=NKVAL
        JSTATE(5*NSTATE)=NK
        IF (LIN) THEN
          JSTATE(6*NSTATE)=INDX
        ELSE
          JSTATE(6*NSTATE)=NSTATE
        ENDIF
        NKVAL=NKVAL+NK
        GOTO 2000
C
C  * * * END OF FILE CONDITIONS * * *
 9000   IF (LEVIN) GOTO 2200

        WRITE(6,606) IUNIT,NSTATE
  606   FORMAT(/'  END OF FILE ENCOUNTERED ON UNIT',I4,'   AFTER',I5,
     &         '  FUNCTIONS.')
        GOTO 2400

 2200   WRITE(6,607) IUNIT,NSTATE
  607   FORMAT(/'  PREMATURE E.O.F. ON UNIT',I4,
     &         '.  NLEVEL REDUCED TO',I6)
        GOTO 2400

 9100   WRITE(6,608) IUNIT,NSTATE
  608   FORMAT(/'  * * * ERROR.   E.O.F. ON UNIT',I4,
     &         '  BEFORE ATAU CARDS FOR NSTATE =',I5)
        WRITE(6,699)
  699   FORMAT(/'  * * * TERMINAL ERROR.')
        STOP
 2000 CONTINUE
C  THIS COMPLETES READ(IASYMU) LOOP
C
 2400 IF (LIN) THEN
C  WE FILTERED ON JLEVEL(), MAKE SURE WE HAVE THEM ALL
        IF (NSTATE.NE.NLEVEL) THEN
          WRITE(6,*) '  ALL LEVELS SPECIFIED BY JLEVEL() WERE NOT FOUND'
          WRITE(6,*) '  *** TERMINAL ERROR.'
          STOP
        ENDIF
C  MAKE SURE EACH VALUE IS THERE AND REORDER IF NECESSARY
C  SO THAT JSTATE(I,6)=I (EXPECTED BY PRBR, EG)
        DO 2409 I=1,NLEVEL
          DO 2408 IX=1,NSTATE
            IF (I.NE.JSTATE(5*NSTATE+IX)) GOTO 2408
            IF (I.EQ.IX) GOTO 2409

            DO 2407 IC=1,6
              ITMP=JSTATE((IC-1)*NSTATE+I)
              JSTATE((IC-1)*NSTATE+I)=JSTATE((IC-1)*NSTATE+IX)
 2407         JSTATE((IC-1)*NSTATE+IX)=ITMP
            GOTO 2409

 2408     CONTINUE
          WRITE(6,684) I,JLEVEL(2*I-1),JLEVEL(2*I)
  684     FORMAT('  INPUT SET',I4,'  J, TAU =',2I5,
     &           '  NOT FOUND ON IASYMU'/'  *** TERMINAL ERROR')
          STOP
 2409   CONTINUE
      ELSE
C  SET J,TAU INTO JLEVEL; GET JMIN,JMAX
        NLEVEL=NSTATE
C  SET JLEVEL(), JMIN, AND JMAX.
        JMIN=JSTATE(1)
        JMAX=JMIN
        DO 2401 I=1,NSTATE
          JI=JSTATE(I)
          JLEVEL(2*I-1)=JI
          JLEVEL(2*I)=JSTATE(I+NSTATE)
          JMIN=MIN(JMIN,JI)
 2401     JMAX=MAX(JMAX,JI)
      ENDIF
C  CORRECT JSTATE(LEV,4) FOR SPACE TAKEN BY JSTATE. . .
      IF (IOFF.NE.6*NSTATE) THEN
        WRITE(6,698) IOFF, NSTATE
  698   FORMAT('  SET6. INDEXING ERROR.  IOFF,NSTATE =',2I6)
        STOP
      ENDIF
      IX=3*NSTATE+1
      IXTOP=4*NSTATE
      DO 2410 I=IX,IXTOP
 2410   JSTATE(I)=JSTATE(I)+IOFF
C  CHECK THAT FUNCTIONS ARE ORTHOGONAL
      CALL CHECK6(NSTATE,JSTATE,ATAU)
C  CHECK THAT ENERGIES ARE NOT ALL IDENTICALLY ZERO.
      DO 2500 I=1,NSTATE
        IF (ELEVEL(I).NE.0.D0) GOTO 2510

 2500 CONTINUE
      IF (NLEVEL.GT.1) THEN
        IF (IPRINT.GE.1) WRITE(6,609)
  609   FORMAT('  *** WARNING.  SET6. ENERGIES ARE ALL ZERO')
      ENDIF
 2510 IXNEXT=IXNEXT+NKVAL

      RETURN
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * END OF SET6 * *
C
C  THESE ENTRY POINTS COMPUTE COUPLING MATRIX ELEMENTS . . .
C
      ENTRY CPL16(N,MXLAM,LAM,NSTATE,JSTATE,JSINDX,L,VL,IV,ATAU)
      IGO1=3003
      IGO2=3033
      L20=.FALSE.
      GOTO 3000
C
C ------------------------------------------------------ END OF CPL16
      ENTRY CPL6(N,MXLAM,LAM,NSTATE,JSTATE,JSINDX,L,JTOT,ATAU,
     1           VL,IPRINT,LFIRST)
      IF (LFIRST) THEN
        IFIRST=-1
        LFIRST=.FALSE.
        NOMEM=.FALSE.
      ENDIF
      L20=.FALSE.
      IF (IFIRST.GT.-1) GOTO 5500

      IF (NOMEM) GOTO 5900

      NL12=NSTATE*(NSTATE+1)/2
      IXMX=NL12*MXLAM
      ISTART=MX+1
      NAVAIL=ISTART-IXNEXT
      IF (IXMX.LE.NAVAIL) GOTO 5100

      IF (IPRINT.GE.3) WRITE(6,694) IXMX,NAVAIL
  694 FORMAT(/'  CPL6 (MAR 93).   UNABLE TO STORE JTOT-INDEPENDENT PART'
     1       /'                   REQUIRED AND AVAILABLE STORAGE =',2I9)
      NOMEM=.TRUE.
      GOTO 5900

 5100 IX=0
      DO 5200 LL=1,MXLAM
        LM=LAM(2*LL-1)
        XLM=LM
        MU=LAM(2*LL)
        XMU=MU
        DO 5201 IC=1,NSTATE
          JC  =JSTATE(IC)
          XJC=JC
          ISTC=JSTATE(IC+3*NSTATE)
          NKC =JSTATE(IC+4*NSTATE)
        DO 5201 IR=1,IC
          IX=IX+1
          JR  =JSTATE(IR)
          XJR=JR
          ISTR=JSTATE(IR+3*NSTATE)
          NKR =JSTATE(IR+4*NSTATE)
          XCPL=Z0
          KKC=-JC
          DO 5300 KC=1,NKC
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
            IF (ABS(ATAU(ISTC+KC)).LE.EPS) GOTO 5300

            XKC=KKC
            KKR=-JR
            DO 5400 KR=1,NKR
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
              IF (ABS(ATAU(ISTR+KR)).LE.EPS) GOTO 5400

              XKR=KKR
              AF=ATAU(ISTR+KR)*ATAU(ISTC+KC)
              IF (LODD(KKR)) AF=-AF
              IF (KKR-KKC.NE.MU) GOTO 5401

              XCPL=XCPL+AF*THRJ(XJC,XJR,XLM,XKC,-XKR,XMU)
              IF (MU.EQ.0) GOTO 5400
 5401         IF (KKC-KKR.NE.MU) GOTO 5400

C  ADJUST FOR (-1)**MU IN POTENTIAL. . .
              IF (LODD(MU)) AF=-AF
              XCPL=XCPL+AF*THRJ(XJC,XJR,XLM,XKC,-XKR,-XMU)
 5400         KKR=KKR+1
 5300       KKC=KKC+1
C  NOW GET 'CONSTANT FACTORS'
          XFCT=PARSGN(JC+JR)*SQRT((F(JC)*F(JR)*F(LM))/(4.D0*PI))
 5201     X(ISTART-IX)=XCPL*XFCT
 5200 CONTINUE

      IF (IPRINT.GE.4) WRITE(6,695) IXMX
  695 FORMAT(/'  CPL6 (MAR 93).   JTOT-INDEPENDENT PARTS OF COUPLING',
     1        ' MATRIX STORED.'/
     2        '                   REQUIRED STORAGE =',I8)
C  RESET MX, IFIRST TO REFLECT STORED VALUES
      MX=MX-IXMX
      IFIRST=0
C
C  NOW GET COUPLING MATRIX ELEMENTS FROM STORED PARTS
 5500 PJT=PARSGN(JTOT)
      IF (IVLU.GT.0) REWIND IVLU
      DO 5600 LL=1,MXLAM
        LM=LAM(2*LL-1)
        MU=LAM(2*LL)
C
C  STORAGE FOR 3J AND 6J SYMBOLS
C
        ITL=IXNEXT
        IT6=ITL+2*LM+1
        IXNEXT=IT6+2*LM+1
        NUSED=0
        CALL CHKSTR(NUSED)
C
        IX1=(LL-1)*NL12
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
C
        LSAV=-1
        DO 5700 IC=1,N
          LEVC=JSINDX(IC)
          JC=JSTATE(LEVC)
          LC=L(IC)
          IF (LC.NE.LSAV) THEN
            CALL J3J000(DBLE(LC),DBLE(LM),IVALL,X(ITL),XLMIN)
            LMIN=ABS(LC-LM)
            LMAX=LC+LM
            LSAV=LC
          ENDIF
C
          LSAV6=-1
        DO 5700 IR=1,IC
          LEVR=JSINDX(IR)
          JR=JSTATE(LEVR)
          LR=L(IR)
C
          IF (LEVR.GE.LEVC) THEN
            IX2=LEVR*(LEVR-1)/2+LEVC
          ELSE
            IX2=LEVC*(LEVC-1)/2+LEVR
          ENDIF
          INDX=IX1+IX2
C
          IF (X(ISTART-INDX).EQ.0.D0 .OR. LR.LT.LMIN .OR. LR.GT.LMAX
     2                               .OR. LODD(LR+LMAX)) THEN
            VL(IX)=0.D0
          ELSE
            IF (LR.NE.LSAV6) THEN
              IVAL6=MX-IT6+1
              CALL J6J(DBLE(LR),DBLE(JTOT),DBLE(LC),DBLE(JC),DBLE(LM),
     1                 IVAL6,XJMIN6,X(IT6))
              JMIN6=INT(XJMIN6)
              LSAV6=LR
            ENDIF
            IF (JR.LT.JMIN6 .OR. JR.GE.JMIN6+IVAL6) THEN
              VL(IX)=0.D0
            ELSE
              INDL=ITL+(LR-LMIN)/2
              IND6=IT6+JR-JMIN6
              VL(IX)=PJT*SQRT(F(LC)*F(LR))*X(ISTART-INDX)*
     &                   X(INDL)*X(IND6)
            ENDIF
          ENDIF
          IF (VL(IX).NE.0.D0) NNZ=NNZ+1
          IF (IVLU.EQ.0) THEN
            IX=IX+MXLAM
          ELSE
            IX=IX+1
          ENDIF
 5700   CONTINUE
        IF (NNZ.EQ.0) WRITE(6,697) LM,MU
        IF (IVLU.GT.0) WRITE(IVLU) (VL(I),I=1,N*(N+1)/2)
        IXNEXT=ITL
 5600 CONTINUE
      RETURN
C
C  IF WE CANNOT STORE PARTIAL COUPLING MATRIX, RECALCULATE.
 5900 IGO1=3001
      IGO2=3011
      GOTO 3000
C--------------------------------------------------------- END OF CPL6
C
      ENTRY CPL26(N,MXLAM,LAM,NSTATE,JSTATE,JSINDX,MVAL,ATAU,
     1            VL,IPRINT,LFIRST)
C
C  IF LFIRST IS TRUE (FIRST CALL), DO SOME INITIALIZATION
      IF (LFIRST) THEN
        IFIRST=-1
        LFIRST=.FALSE.
        NOMEM=.FALSE.
      ENDIF
      L20=.TRUE.
C
      IF (IFIRST.GT.-1) GOTO 4500

C  FIRST TIME THROUGH SET UP SOME STORAGE POINTERS
      NL12=NSTATE*(NSTATE+1)/2
      IXMX=NL12*MXLAM
      ISTART=MX+1
C
 4500 MVABS=ABS(MVAL)
C  SEE IF VALUES ARE STORED FOR THIS HIGH AN MVALUE
C  IF NOT, TRY TO STORE THEM IN XCPL().
      IF (MVABS.LE.IFIRST .OR. NOMEM) GOTO 4900

      MV=IFIRST+1
C  FIRST CHECK THAT WE STILL HAVE A CONTINUOUS BLOCK OF HI MEMORY.
 4600 IF (MX.EQ.ISTART-(IFIRST+1)*IXMX-1) GOTO 4610

      IF (IPRINT.GE.1) WRITE(6,642) MV,ISTART-1,MX,IXMX*(IFIRST+1)
  642 FORMAT(/'  CPL26 (FEB 93).  HIGH MEMORY FRAGMENTED.  CANNOT',
     1       ' STORE COUPLING COEFFS FOR MVAL =',I3/ 19X,'ORIGINAL ',
     2       'MINUS CURRENT MEMORY LIMITS .NE. NO. USED =',3I12)
      NOMEM=.TRUE.
      GOTO 4900

C  TEST FOR AVAILABLE STORAGE; NEED IXMX FOR THIS MVAL
 4610 NAVAIL=MX-IXNEXT+1
      IF (IXMX.LE.NAVAIL) GOTO 4601

      IF (IPRINT.GE.3) WRITE(6,692) MV,IXMX,NAVAIL
  692 FORMAT(/'  CPL26 (FEB 93).   UNABLE TO STORE 3-J VALUES FOR ',
     1        'MVAL =',I3/
     2        '                    REQUIRED AND AVAILABLE STORAGE =',
     3        2I9)
C  SET NOMEM TO REFLECT INABILITY TO ADD MORE M-VALUES
      NOMEM=.TRUE.
      GOTO 4900
C
C  REDUCE 'TOP OF MEMORY' AND STORE COUPLING VALUES FOR THIS MVAL
 4601 MX=MX-IXMX
C  START INDEX AFTER M-BLOCKS ALREADY STORED (STARTING WITH MV=0)
      IX=MV*IXMX
      NZERO=0
      DO 4200 LL=1,MXLAM
        LM=LAM(2*LL-1)
        MU=LAM(2*LL)
        DO 4201 IC=1,NSTATE
          JC  =JSTATE(IC)
          ISTC=JSTATE(IC+3*NSTATE)
          NKC =JSTATE(IC+4*NSTATE)
        DO 4201 IR=1,IC
          JR  =JSTATE(IR)
          ISTR=JSTATE(IR+3*NSTATE)
          NKR =JSTATE(IR+4*NSTATE)
          IX=IX+1
          XCPL=Z0
          KKC=-JC
          DO 4300 KC=1,NKC
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
            IF (ABS(ATAU(ISTC+KC)).LE.EPS) GOTO 4300

            KKR=-JR
            DO 4400 KR=1,NKR
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
              IF (ABS(ATAU(ISTR+KR)).LE.EPS) GOTO 4400

              AF=ATAU(ISTR+KR)*ATAU(ISTC+KC)
              IF (KKR-KKC.NE.MU) GOTO 4401

              XCPL=XCPL+AF*GSYMTP(JC,KKC,JR,KKR,MV  ,LM,MU)
              IF (MU.EQ.0) GOTO 4400
 4401         IF (KKC-KKR.NE.MU) GOTO 4400

C  ADJUST FOR (-1)**MU IN POTENTIAL. . .
              IF (LODD(MU)) AF=-AF
              XCPL=XCPL+AF*GSYMTP(JC,KKC,JR,KKR,MV,LM,-MU)
 4400         KKR=KKR+1
 4300       KKC=KKC+1
 4201     X(ISTART-IX)=XCPL
 4200 CONTINUE
      IF (IPRINT.GE.4) WRITE(6,693) MV,IXMX,NAVAIL
  693 FORMAT(/'  CPL26 (FEB 93).   3-J VALUES STORED FOR MVAL =',I3
     1       /'                    REQUIRED AND AVAILABLE STORAGE =',
     2       2I9)
C  RESET IFIRST TO REFLECT HIGHEST M-VALUE STORED.
      IFIRST=MV
C  SEE IF CURRENT MVALUE REQUIRES MORE STORED M-VALUES.
      MV=MV+1
      IF (MV.LE.MVABS) GOTO 4600
C
 4900 IF (MVABS.GT.IFIRST) GOTO 4800

C  MVABS.LE.IFIRST.  COEFFS STORED.  FILL VL() FROM XCPL
      IXM=MVABS*IXMX
      IF (IVLU.GT.0) REWIND IVLU
      DO 4513 LL=1,MXLAM
        LM=LAM(2*LL-1)
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
        DO 4503 ICOL=1,N
          I1=JSINDX(ICOL)
          J1=JSTATE(I1)
        DO 4503 IROW=1,ICOL
          I2=JSINDX(IROW)
          J2=JSTATE(I2)
          IF (I1.GT.I2) THEN
            IX12=I1*(I1-1)/2+I2
          ELSE
            IX12=I2*(I2-1)/2+I1
          ENDIF
          IXX=IXM+(LL-1)*NL12+IX12
          VL(IX)=X(ISTART-IXX)
C  WE HAVE STORED COUPLING FOR POSITIVE MVALUES; CORRECT IF NECESSARY
C  FOR PARITY OF THRJ(J1,LM,J2,-MVAL,0,MVAL)
          IF (MVAL.LT.0 .AND. LODD(J1+J2+LM)) VL(IX)=-VL(IX)
          IF (VL(IX).NE.Z0) NNZ=NNZ+1
          IF (IVLU.EQ.0) THEN
            IX=IX+MXLAM
          ELSE
            IX=IX+1
          ENDIF
 4503   CONTINUE
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) MVAL,LL
        ENDIF
  612   FORMAT('  * * * NOTE.  FOR MVALUE =',I4,',  ALL COUPLING '
     1         'COEFFICIENTS ARE 0.0 FOR EXPANSION TERM',I4)
        IF (IVLU.GT.0) WRITE(IVLU) (VL(I),I=1,N*(N+1)/2)
 4513 CONTINUE

      RETURN
C
C  MV.GT.IFIRST ==> VALUES NOT STORED.  CALCULATE THEM VIA OLD CODE
 4800 IGO1=3002
      IGO2=3022
      GOTO 3000
C
C  -------------------- OLD CODE REJOINS HERE ---------------------
C
 3000 IF (IVLU.GT.0) REWIND IVLU
      NZERO=0
      DO 3100 LL=1,MXLAM
        LM=LAM(2*LL-1)
        MU=LAM(2*LL)
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
C
        DO 3200 IC=1,N
          JC  =JSTATE(JSINDX(IC))
          ISTC=JSTATE(JSINDX(IC)+3*NSTATE)
          NKC =JSTATE(JSINDX(IC)+4*NSTATE)
        DO 3200 IR=1,IC
          JR  =JSTATE(JSINDX(IR))
          ISTR=JSTATE(JSINDX(IR)+3*NSTATE)
          NKR =JSTATE(JSINDX(IR)+4*NSTATE)
C
          VL(IX)=0.D0
          KKC=-JC
          DO 3300 KC=1,NKC
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
            IF (ABS(ATAU(ISTC+KC)).LE.EPS) GOTO 3300

            KKR=-JR
            DO 3400 KR=1,NKR
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
              IF (ABS(ATAU(ISTR+KR)).LE.EPS) GOTO 3400

              AF=ATAU(ISTR+KR)*ATAU(ISTC+KC)
              IF (KKR-KKC.NE.MU) GOTO 3500

              IF (IGO1.EQ.3001) THEN
                VL(IX)=VL(IX)+
     1                   AF*FSYMTP(JC,KKC,L(IC),JR,KKR,L(IR),JTOT,LM,MU)

              ELSEIF (IGO1.EQ.3002) THEN
                VL(IX)=VL(IX)+AF*GSYMTP(JC,KKC,JR,KKR,MVAL,LM,MU)

              ELSEIF (IGO1.EQ.3003) THEN
                VL(IX)=VL(IX)+AF*ESYMTP(JC,KKC,JR,KKR,LM,MU)
              ENDIF
C
              IF (MU.EQ.0) GOTO 3400
 3500         IF (KKC-KKR.NE.MU) GOTO 3400

C  ADJUST FOR (-1)**MU IN POTENTIAL. . .
              AF=AF*PARSGN(MU)
              IF (IGO2.EQ.3011) THEN
                VL(IX)=VL(IX)+
     &                  AF*FSYMTP(JC,KKC,L(IC),JR,KKR,L(IR),JTOT,LM,-MU)

              ELSEIF (IGO2.EQ.3022) THEN
                VL(IX)=VL(IX)+AF*GSYMTP(JC,KKC,JR,KKR,MVAL,LM,-MU)

              ELSEIF (IGO2.EQ.3033) THEN
                VL(IX)=VL(IX)+AF*ESYMTP(JC,KKC,JR,KKR,LM,-MU)
C
              ENDIF
 3400         KKR=KKR+1
 3300       KKC=KKC+1
          IF (VL(IX).NE.0.D0) NNZ=NNZ+1
          IF (IVLU.EQ.0) THEN
            IX=IX+MXLAM
          ELSE
            IX=IX+1
          ENDIF
 3200   CONTINUE
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,697) LM,MU
        ENDIF
  697   FORMAT('  * * * NOTE.  ALL COUPLING COEFFICIENTS ARE ZERO '
     1         ' FOR LM, MU = ',2I4)
        IF (IVLU.GT.0) WRITE(IVLU) (VL(I),I=1,N*(N+1)/2)
 3100 CONTINUE

      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14) THEN
        IF (L20) WRITE(6,620) 'MVAL',MVAL,NZERO
        IF (.NOT.L20) WRITE(6,620) 'JTOT',JTOT,NZERO
      ENDIF
  620 FORMAT('  * * * NOTE.  FOR ',A,' =',I4,',  ALL COUPLING ',
     1       'COEFFICIENTS ARE 0.0 FOR',I5,' POTENTIAL EXPANSION TERMS')

      RETURN
      END
C -------------------------------------------------------- END OF CPL26
