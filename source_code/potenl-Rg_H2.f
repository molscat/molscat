      SUBROUTINE POTENL(IC, MXLMB, LMB, RR, P, ITYP, IPRINT)
C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE potential, ONLY: LAMBDA, MXLMDA, RMNAME, EPNAME
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  ROUTINE TO HANDLE THE INTERMOLECULAR POTENTIAL
C  FOR AN ATOM + VIBRATING DIATOM SYSTEM (ITYP=10*N+7)
C  THIS VERSION IS FOR THE RG-H2 POTENTIALS OF
C  LE ROY & CARLEY (ADV CHEM PHYS 42, 353 (1980) AND
C  LE ROY & HUTSON (J CHEM PHYS 86, 837 (1987)).
C  IT IS RESTRICTED TO V=0 AND 1 BUT EASILY GENERALISED.
C  FULL SETS OF EXPECTATION VALUES OF THE INTERMOLECULAR POTENTIAL
C  BETWEEN (V,J) AND (V',J') DIATOM INTERNAL STATES ARE RETURNED.
C
      CHARACTER(80) FILNAM
      DIMENSION LMB(5,45), P(MXLMB)
      DIMENSION VLAM(15), VT(4), VPRM(4,2,4), VPRMR(4,2,4), BETA(2),
     1          AA(2,210), C8(2,210), C6(2,210), XI(9,210), IV(210),
     2          IVD(210), IJ(210), DEL(210), XG(20), WG(20), SM(15),
     3          FSM(15,15), PWT(16,16), MP(2), LAMB(2), D(4)
      DATA Z0/0.D0/,ZH/0.5D0/,Z1/1.D0/,Z2/2.D0/,Z4/4.D0/,Z6/6.D0/,
     1     Z8/8.D0/,Z12/12.D0/
      DATA MXJ/210/
      DATA LAMB/0,2/,X0/0.7666348D0/,DREL/-0.1665563536D0/
      DATA FILNAM /'data/h2even.dat'/
      NAMELIST /POTL/ NLEG,NSTR,IHET,JMAXV0,JMAXV1,NGP,SCAL,SCALMX,
     1                FILNAM
      P2(X) = 1.5D0*X*X-ZH
C
C  IC IS A CONTROL PARAMETER:
C  IC=-1  INITIALISATION OF POTENTIAL ETC.
C  IC=0   CALCULATION OF POTENTIAL
C
      IF (IC.EQ.0) GOTO 530
      IF (IC.EQ.-1) GOTO 120
      WRITE(6,110) IC
  110 FORMAT ('  ***** ERROR IN POTENL - IC =',I6)
      STOP
C-----------------------------------------------------------------------
C
C  INITIALIZE. PARAMETERS HAVE THE FOLLOWING MEANINGS ON ENTRY:
C  LMB     AN EMPTY INTEGER ARRAY
C  MXLMB   THE TOTAL NUMBER OF 4-BYTE WORDS IN LMB
C  RR      NOT USED
C  P       NOT USED
C  ITYP    THE COLLISION TYPE (MUST BE 10*N+7 FOR THIS ROUTINE)
C
C  FIRST CHECK IF ITYP VALUE IS ACCEPTABLE
C
  120 ITYPE = ITYP-10*(ITYP/10)
      IF (ITYPE.EQ.7) GOTO 140
      WRITE(6,130) ITYP
  130 FORMAT ('  ***** ERROR IN POTENL - CANNOT HANDLE ITYPE =',I6)
      STOP
C
C  READ &POTL DATA. THE PARAMETERS HAVE THE FOLLOWING MEANINGS:
C  NLEG   THE TOTAL NUMBER OF LEGENDRE POLYNOMIALS IN THE ANGULAR
C         EXPANSION OF THE POTENTIAL. THE POLYNOMIAL INDICES RUN
C         OVER 0,1,2, .... NLEG-1  FOR A HETERONUCLEAR DIATOM
C         AND OVER 0,2,4, .... 2*(NLEG-1) FOR A HOMONUCLEAR DIATOM.
C  NSTR   THE NUMBER OF TERMS IN THE DIATOM STRETCH EXPANSION
C         OF EACH LEGENDRE COMPONENT THAT ARE INCLUDED IN THE
C         QUADRATURE TO AVERAGE OVER V AND V' WHEN PERFORMING A
C         CENTRE-OF-MASS SHIFT (FOR HD).
C  IHET   IHET=0 FOR A HOMONUCLEAR DIATOM, IHET=1 FOR A
C         HETERONUCLEAR DIATOM
C  JMAXV0  THE POTENTIAL ARRAY RETURNED BY POTENL WILL CONTAIN
C          TERMS CORRESPONDING TO V=0 DIATOM STATES UP TO JMAXV0.
C          SUBJECT TO THIS LIMIT, THE ROUTINE WILL INCLUDE ALL
C          V=0 J'S FOR WHICH DIATOM MATRIX ELEMENTS ARE PROVIDED.
C          SET JMAXV0 NEGATIVE IF NO V=0 STATES ARE REQUIRED.
C  JMAXV1  AS FOR JMAXV0, BUT FOR V=1 DIATOM STATES.
C  NGP     NO. OF GAUSSIAN QUADRATURE POINTS USED FOR INTEGRATING
C          OVER COS THETA AND XI FOR CENTRE-OF-MASS SHIFT
C  SCAL    IF SCAL > 1.0, MATRIX ELEMENTS OFF-DIAGONAL IN V WILL
C          BE INCREASED BY A FACTOR OF SCAL. THIS IS USEFUL WHEN
C          SEARCHING FOR VERY NARROW RESONANCES.
C  SCALMX  LARGEST PERMITTED VALUE OF OFF-DIAGONAL MATRIX ELEMENTS
C          WHEN SCAL OPTION IS BEING USED. HELPS TO PREVENT
C          UNPHYSICAL BEHAVIOUR OF THE POTENTIAL.
C  FILNAM  NAME OF FILE CONTAINING MATRIX ELEMENTS OF POWERS OF XI
C
  140 IF (IPRINT.GE.1)
     1  WRITE(6,'(2X,A,A)') 'POTENL ROUTINE FOR RARE GAS ATOM + H2 ',
     2                      '(UPDATED APR 2019)'
C  INITIALISE NAMELIST VARIABLES BEFORE READ
      IHET = 0
      NLEG = 1
      NSTR = 1
      JMAXV0 = 9999
      JMAXV1 = 9999
      NGP = 10
      SCAL = 0.D0
      SCALMX = 1.D5
C-----------------------------------------------------------------------
      READ(5,POTL)
C-----------------------------------------------------------------------
      IHET = MIN(1,MAX(0,IHET))
C  NOFK  IS THE NUMBER OF STRETCHING-DEPENDENT TERMS IN FINAL POTENTIAL
C  FOLLOWING TRANSFORMATION (FOR (IHET.GT.0)).  ***
      NOFK = MIN(9,NSTR)
      IF (IHET.EQ.0) NOFK = 4
      IF (IPRINT.GE.1)
     1  WRITE(6,150) NLEG,NSTR,IHET,JMAXV0,JMAXV1,NGP,SCAL,SCALMX
  150 FORMAT (/'  /POTL/ DATA IS   NLEG =',I3,'   NSTR =',I3,
     1        '   IHET =',I2,'   JMAXV0 =',I3,'    JMAXV1 =',I3/
     2        19X,'NGP  =',I3,'    SCAL =',F8.1,'  SCALMX =',F8.1)
C
C  INITIALIZE POTENTIAL AND SET EXIT VALUES OF ALL PARAMETERS. EXIT
C  VALUES ARE DEFINED AS FOLLOWS (WHERE NINPOT IS THE TOTAL NUMBER
C  OF DIATOM VIB-ROTATIONAL LEVELS ACTUALLY INCLUDED IN THE POTENTIAL
C  ARRAY RETURNED BY THIS ROUTINE):
C
C  LMB     ON EXIT, LMB MUST CONTAIN SETS OF INDICES SPECIFYING
C          EACH SYMMETRY TERM IN THE POTENTIAL. IN EFFECT, LMB IS
C          A 2-DIMENSIONAL ARRAY LMB(5 , NLEG*NINPOT*(NINPOT+1)/2).
C          LMB(1,I) GIVES THE LEGENDRE POLYNOMIAL INDEX, AND
C          LMB(2,I),LMB(3,I),LMB(4,I) AND LMB(5,I) GIVE THE DIATOM
C          V,J AND V',J' QUANTUM NUMBERS RESPECTIVELY. EACH OF
C          THESE SETS OF FIVE INDICES DEFINES A SYMMETRY TERM.
C          MOLSCAT DOES NOT REQUIRE THE SYMMETRY TERMS TO BE GIVEN
C          IN ANY PARTICULAR ORDER, BUT JUST THAT THE I-TH ELEMENT
C          OF THE POTENTIAL ARRAY RETURNED BY THIS ROUTINE SHOULD
C          CORRESPOND TO THE I-TH SET OF INDICES IN LMB.
C          WE CHOOSE HERE TO ARRANGE LMB WITH THE LEGENDRE
C          POLYNOMIAL INDEX VARYING FASTEST, THEN THE (V,J) INDICES
C          AND THEN THE (V',J') INDICES. THE ORDERING OF (V,J)
C          AND (V',J') INDICES FOLLOWS THE ORDER IN WHICH DIATOM
C          MATRIX ELEMENTS ARE READ IN (SEE BELOW).
C  MXLMB   EQUAL TO NLEG*NINPOT*(NINPOT+1)/2, I.E. THE TOTAL NUMBER
C          OF SYMMETRY TERMS IN THE POTENTIAL ARRAY.
C  RR      GIVES RM, THE UNIT OF LENGTH (IN ANGSTROMS)
C  P(1)    GIVES EPSIL, THE UNIT OF ENERGY (IN WAVENUMBERS)
C
C  SET RM AND EPSIL
C
      RR = 1.D0
      P(1) = 1.D0
      IF (IPRINT.GE.1) WRITE(6,160)
  160 FORMAT ('  LENGTH UNITS ARE ANGSTROMS.    ENERGY UNITS ARE ',
     1        'WAVENUMBERS (1/CM)')
C***********************************************************************
C**                                                                  ***
C** R J LE ROY POTENTIAL INITIALIZATION CODE                         ***
C**                                                                  ***
C***********************************************************************
C** GENERATE LEGENDRE EXPANSION COMPONENTS OF POTENTIAL  VLAM(L)  (CM-1)
C** ON THIS ENTRY, READ IN PARAMETERS OF  GBC3(6,8)-TYPE POTENTIAL FOR
C                  HOMONUCLEAR DIATOM & TRANSFORM TO SHIFTED CENTRE OF
C                  MASS  IF (IHET.GT.0)
C-----------------------------------------------------------------------
C** NOFPOT  IS THE NUMBER OF VIB/ROT THRESHOLDS FOR WHICH DIATOM
C           MATRIX ELEMENTS ARE SUPPLIED.
C** IDAMP=1 FOR THE BC DAMPING FUNCTION (LE ROY/CARLEY 1980)
C   IDAMP=2 FOR THE TT DAMPING FUNCTION (LE ROY/HUTSON 1987)
      READ (5,*) NOFPOT,IDAMP
C** READ POWER  MP(LAM)  & EXPONENTIAL CONSTANT  BETA(LAM) FOR POTENTIAL
C   REPULSION TERM  R**(-MP(LAM))*EXP(-BETA(LAM)*R) .
      READ (5,*) (MP(L),BETA(L),L=1,2),R0
C-----------------------------------------------------------------------
      IF (BETA(1).EQ.BETA(2)) GOTO 5
      WRITE(6,601)
  601 FORMAT('  *** ERROR *** BETA(1) .NE. BETA(2) NOT IMPLEMENTED',
     1       ' IN THIS VERSION OF POTENL')
      STOP
C
    5 IF (NOFPOT*(NOFPOT+1)/2.LE.210) GOTO 165

      WRITE(6,163) NOFPOT
  163 FORMAT('   POTENL: NOFPOT =',I4,' TOO GREAT FOR ARRAY DIMENSIONS')
      STOP

  165 IF (IPRINT.GE.1) WRITE(6,430)
      IF (IDAMP.EQ.1 .AND. IPRINT.GE.1) WRITE(6,431) R0
      IF (IDAMP.EQ.2 .AND. IPRINT.GE.1) WRITE(6,432)
      IF (IPRINT.GE.1) WRITE(6,433)
      DO 210 LAM=1,2
        MM = MP(LAM)
        DO 190 K=1,4
          IK = K-1
C**  VT(I)  ARE  EPSILON, REQ & C6, FOR I=1-3, FOR THIS (LAM,K)
C-----------------------------------------------------------------------
          READ (5,*) (VT(I),I=1,3)
C-----------------------------------------------------------------------
          IF (VT(2).LE.Z0) GOTO 180

          EPS = VT(1)
          RE = VT(2)
          BETRE = BETA(LAM)*RE
          C6RE = VT(3)/RE**6
          P3 = EXP(BETRE)
          CALL DAMP(IDAMP,RE,R0,BETA(LAM),D)
          P4 = Z8-BETRE-RE*D(3)
          VPRM(1,LAM,K) = EPS
          VPRMR(1,LAM,K) = ((Z2+(D(1)-D(3))*RE)*C6RE*D(2)
     1                     + (RE*D(3)-Z8)*EPS)*P3/P4
          VPRM(2,LAM,K) = RE
          VPRMR(2,LAM,K) = (VPRMR(1,LAM,K)/P3+EPS-C6RE*D(2))*RE**8/D(4)

  180     VPRM(3,LAM,K) = VT(3)
          VPRMR(3,LAM,K) = VT(3)
          IF (VT(2).LE.Z0) GOTO 200

  190     IF (IPRINT.GE.1) WRITE(6,440) LAMB(LAM),IK,MP(LAM),BETA(LAM),
     1                                  (VPRM(I,LAM,K),I=1,3),
     2                                  VPRMR(2,LAM,K),VPRMR(1,LAM,K)
        GOTO 210

  200   VPRMR(1,LAM,K) = Z0
        VPRMR(2,LAM,K) = Z0
        IF (IPRINT.GE.1) WRITE(6,450) LAMB(LAM),VPRM(3,LAM,4)
  210 CONTINUE
      IF (IPRINT.GE.1) WRITE(6,470)
C** READ IN THE VIBRATIONAL & ROTATIONAL QUANTUM NUMBERS  IV & IJ  AND
C     EXPECTATION VALUES AND (FOR NOFPOT>1) MATRIX ELEMENTS OF POWERS  K
C     TO  8  OF THE STRETCHING COORDINATE  XI .
C** IN ORDER: DIAGONAL TERMS  <V,J:XI**K:V,J>, K=0 TO 8  IN ORDER OF
C     INCREASING THRESHOLD ENERGY, EACH IMMEDIATELY FOLLOWED BY MATRIX
C     ELEMENTS  <V,J:XI**K:V',J'> ,K=0 TO 8  COUPLING IT TO EACH HIGHER
C     THRESHOLD  (V',J') .
C** WHILE READING THESE QUANTITIES, IGNORE ANY WHICH ARE EXCLUDED BY
C     JMAXV0 AND JMAXV1. ALSO, COMPUTE NINPOT, THE NUMBER OF DIATOM STAT
C     WHICH WILL ACTUALLY BE INCLUDED IN THE POTENTIAL ARRAY. ALSO, SET
C     UP THE LMB ARRAY.
C** IF ((IHET.LE.0) .OR. (ISTR.LE.1)) PREPARE VIBRATIONALLY AVERAGED
C     STRENGTH COEFFICIENTS FOR POTENTIAL MATRIX ELEMENTS
      JJMX = NOFPOT
      NINPOT = 0
      KLM = 0
      IFAC = 2-IHET
      J = 1
      OPEN(1,FILE=FILNAM,STATUS='OLD',ERR=999)
  220 DO 320 JJ=1,JJMX
C-----------------------------------------------------------------------
        READ (1,*) IV(J),IJ(J),IVD(J),IJD,(XI(I,J),I=1,9)
C-----------------------------------------------------------------------
C** CHECK FOR EXCLUSION OF THIS (V,J) AND SET UP LMB
        IF (IVD(J).EQ.0.AND.IJD.GT.JMAXV0) GOTO 320
        IF (IVD(J).EQ.1.AND.IJD.GT.JMAXV1) GOTO 320
        IF (IV(J).EQ.0.AND.IJ(J).GT.JMAXV0) GOTO 320
        IF (IV(J).EQ.1.AND.IJ(J).GT.JMAXV1) GOTO 320

        IF (J.LE.MXJ) GOTO 235
        WRITE(6,231) MXJ
  231   FORMAT ('  ****** ERROR IN POTENL - ARRAYS AA, C6, ETC. (',I4,
     1          ') ARE NOT BIG ENOUGH')
        STOP

  235   IF (5*J*NLEG.LE.MXLMB) GOTO 240
        WRITE(6,236) MXLMB
  236   FORMAT ('  ****** ERROR IN POTENL - ARRAY LMB(',I7,
     1          ') IS NOT BIG ENOUGH')
        STOP

  240   IF (JJ.EQ.1) NINPOT = NINPOT+1
        KLEG = -IFAC
        DO 250 ILEG=1,NLEG
          KLEG = KLEG+IFAC
          KLM = KLM+1
          LMB(1,KLM) = KLEG
          LMB(2,KLM) = IV(J)
          LMB(3,KLM) = IJ(J)
          LMB(4,KLM) = IVD(J)
          LMB(5,KLM) = IJD
  250   CONTINUE
        DEL(J) = Z0
        IF (NSTR.GT.1) GOTO 280

        DEL(J) = (XI(2,J)+Z1)*X0*DREL
        DO 270 LAM=1,2
          AAT = Z0
          C8T = Z0
          C6T = Z0
          DO 260 K=1,4
            XIPW = XI(K,J)
            AAT = AAT+VPRMR(1,LAM,K)*XIPW
            C8T = C8T+VPRMR(2,LAM,K)*XIPW
  260       C6T = C6T+VPRMR(3,LAM,K)*XIPW
          AA(LAM,J) = AAT
          C8(LAM,J) = C8T
          C6(LAM,J) = C6T
  270   CONTINUE
        GOTO 290

C       CONVERT POWERS OF XI TO P_N(XI) FOR GAUSS-LEGENDRE QUADRATURE.
C       THIS VERSION OF THE ARRAYS IS USED ONLY FOR CENTRE-OF-MASS SHIFT
C       AND ONLY WHEN TRANSFORM IS DONE BEFORE INTEGRATION OVER XI.
C
C       THE TRANSFORMATION IS BASED ON LIU ET AL., JCP 68, 5028 (1978)
C       BUT NOTE THAT THEIR XI IS (R-X0)/(R+X0) NOT (R-X0)/X0 AS HERE.
C       USING GAUSS-LEGENDRE QUADRATURE FOR XI IMPLIES A RANGE -1<XI<1
C       WHICH CORRESPONDS TO 0<R<INFTY IN LIU ET AL. BUT 0<R<2*X0 HERE.
C       THIS APPROACH IS THEREFORE INACCURATE FOR VIBRATIONAL STATES
C       WHOSE WAVEFUNCTION EXTENDS OUTSIDE 2*X0.
C
  280   XI(9,J) = (6435*XI(9,J)-12012*XI(7,J)+6930*XI(5,J)-1260*XI(3,J)+
     1            35*XI(1,J))/128.D0
        XI(8,J) = (429*XI(8,J)-693*XI(6,J)+315*XI(4,J)-35*XI(2,J))/16.D0
        XI(7,J) = (231*XI(7,J)-315*XI(5,J)+105*XI(3,J)-5*XI(1,J))/16.D0
        XI(6,J) = (63*XI(6,J)-70*XI(4,J)+15*XI(2,J))/8.D0
        XI(5,J) = (35*XI(5,J)-30*XI(3,J)+3*XI(1,J))/8.D0
        XI(4,J) = (5*XI(4,J)-3*XI(2,J))/2.D0
        XI(3,J) = (3*XI(3,J)-XI(1,J))/2.D0
  290   IF (JJ.GE.2) GOTO 300

        IF (NSTR.LE.1 .AND. IPRINT.GE.1)
     1    WRITE(6,460) IV(J),IJ(J),(LAMB(L),C6(L,J),
     2                              C8(L,J),AA(L,J),L=1,2)
        IF (IPRINT.GE.1) WRITE(6,500) IV(J),IJ(J),(XI(K,J),K=1,NOFK)
        GOTO 310

  300   IF (NSTR.LE.1 .AND. IPRINT.GE.1)
     1    WRITE(6,480) IVD(J),IJD,IV(J),IJ(J),(LAMB(L),C6(L,J),
     2                                         C8(L,J),AA(L,J),L=1,2)
        IF (IPRINT.GE.1) WRITE(6,510) IVD(J),IJD,IV(J),IJ(J),
     1                                (XI(K,J),K=1,NOFK)
  310   J = J+1
  320 CONTINUE

      JJMX = JJMX-1
      IF (JJMX.GT.0) GOTO 220

C     FINISHED WITH CHANNEL 1. CLOSE FILE ON VAX.
      CLOSE(1)
      IF (NINPOT.GT.0) GOTO 340

      WRITE(6,330)
  330 FORMAT ('  ****** ERROR IN POTENL - FINAL LIST OF DIATOM ',
     1        'STATES IS EMPTY - CHECK JMAXV0,JMAXV1,MATRIX ELEMENTS')
      STOP

  340 NOFJ = NINPOT*(NINPOT+1)/2
      MXLMB = NLEG*NOFJ
      IF (NOFJ.EQ.J-1) GOTO 360

      WRITE(6,350) NOFJ,J
  350 FORMAT ('  ****** ERROR IN POTENL - INCORRECT QUANTUM NUMBER ',
     1        'INDEXING',2I7)
      STOP

  360 IF (NSTR.GT.1.OR.NINPOT.EQ.1) GOTO 380

      INOW = NOFJ+1
      DO 370 I=2,NINPOT
        INOW = INOW-1
        JTOP = I-1
      DO 370 J=1,JTOP
        INOW = INOW-1
        ID = NOFJ+1-I*(I+1)/2
        JD = NOFJ+1-J*(J+1)/2
        DEL(INOW) = ZH*(DEL(ID)+DEL(JD))
  370 CONTINUE

  380 CONTINUE
      IF (IHET.LE.0) GOTO 420

      IF (IPRINT.GE.1) WRITE(6,490) DREL
      IF (NSTR.LE.1 .AND. IPRINT.GE.1) WRITE(6,520) (DEL(J),J=1,NOFJ)
C** PREPARE GAUSSIAN WEIGHTS AND POINTS FOR QUADRATURE(S)
C   IN CENTRE-OF-MASS TRANSFORMATION
C   AND ALSO FORM PWT ARRAY CONTAINING PRODUCTS OF LEGENDRE
C   POLYNOMIALS AND GAUSSIAN WEIGHTS.
      IF (NGP.LE.16) GOTO 400

      IF (IPRINT.GE.1) WRITE(6,390)
  390 FORMAT ('  *** TOO MANY GAUSSIAN POINTS REQUESTED FOR ARRAY ',
     1        'DIMENSIONS. NGP RESET TO 16.'/)
      NGP = 16

  400 CONTINUE
C     GET GAUSS-LEGENDRE QUADRATURE POINTS AND WEIGHTS
      CALL GAUSSP(-1.D0,1.D0,NGP,XG,WG)
C     WRITE(6,609) NGP, (XG(I),WG(I),I=1,NGP)
C     609 FORMAT(/2X,'POINTS AND WEIGHTS FOR SIMPLE GAUSSIAN',I3,
C     1          '-PT QUADRATURE'//(2(F25.20,F23.20)))
      DO 410 I=1,NGP
        XX = XG(I)
        PLEG = Z0
        PLEGP = Z1
        A3N = -Z1
        A4N = -Z1
      DO 410 ILEG=1,NGP
        PWT(ILEG,I) = PLEGP*WG(I)
        A3N = A3N+Z2
        A4N = A4N+Z1
        PLEGB = PLEG
        PLEG = PLEGP
  410   PLEGP = (A3N*XX*PLEG-A4N*PLEGB)/DBLE(ILEG)

  420 RETURN
  999 WRITE(6,*) ' *** ERROR: FILE '//TRIM(FILNAM)//' NOT FOUND'
      STOP

  430 FORMAT ('  POTENTIAL REPULSION:  R**(-MM)*EXP(-BETA*R).')
  431 FORMAT ('  DAMPING FUNCTION IS:  D(R)=EXP(-4*(R0/R)**3)',
     1        '  TURNS ON AT   R0 =',F9.6)
  432 FORMAT ('  DAMPING FUNCTION IS:  1 - ] SUM(K=0,N)',
     1        ' (BETA*R)**K/K! ( EXP(-BETA*R)')
  433 FORMAT (/16X,'LAMBDA   K   MM   BETA/A-1   EPS/CM-1   REQ/A',
     1        5X,'C6(L,K)      C8(L,K)     A(L,K)'/16X,84('-'))
  440 FORMAT (14X,2I6,I5,F10.4,F11.4,F10.6,2F12.2,F13.2)
  450 FORMAT (14X,I6,'     3 ',35X,F12.2)
  460 FORMAT (/3X,'V =',I2,', J =',I2,'   VIBRATIONALLY AVERAGED ',
     1        'POTENTIAL:  LAMBDA =',I2,2F12.2,F13.2/61X,I2,2F12.2,
     2        F13.2)
  470 FORMAT (16X,84('-'))
  480 FORMAT (/9X,'(V =',I2,', J =',I2,': V'' =',I2,', J'' =',I2,')',
     1        ' POTENTIAL MATRIX:  LAMBDA =',I2,2F12.2,F13.2/
     2        67X,I2,2F12.2,F13.2)
  490 FORMAT ('  LEGENDRE EXPANSION TRANSFORMATION BASED ON CENTRE-',
     1        'OF-MASS SHIFT OF',F14.10,' TIMES THE BOND LENGTH')
  500 FORMAT (3X,'V =',I2,', J =',I2,'   DIAGONAL EXPECTATION VALUES',
     1        6X,5F12.8/(51X,5F12.8))
  510 FORMAT (9X,'(V =',I2,', J =',I2,': V'' =',I2,', J'' =',I2,')',
     1        ' MATRIX ELEMENTS  ',5F12.8/(57X,5F12.8))
  520 FORMAT (/10X,'YIELDS   DEL(J) =',5F13.8)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  POTENTIAL EVALUATION CALL. ON ENTRY, MXLMB HAS
C  THE VALUE ASSIGNED ON INITIALIZATION (SEE ABOVE).
C  R GIVES THE INTERMOLECULAR DISTANCE (IN ANGSTROMS) AT WHICH EACH
C  SYMMETRY TERM IS TO BE EVALUATED. LMB IS A DUMMY, AND DOES *NOT*
C  CONTAIN THE INDICES WHICH WERE PUT INTO IT IN THE INITIALIZATION
C  CALL. ITYP IS ALSO A DUMMY. THE ARRAY P IS NOT IN USE ON ENTRY.
C
C  ON EXIT, THE ARRAY P(MXLMB) MUST CONTAIN THE REQUIRED POTENTIALS
C  (IN WAVENUMBERS). THE ORDER OF STORAGE IN P MUST BE THE SAME
C  AS THAT IN THE ORIGINAL LMB ARRAY (SEE ABOVE).
C
  530 IF (IHET.GT.0) GOTO 560
C*********************************************************************
C***                                                               ***
C*** BASED ON R J LE ROY'S POTENTIAL ROUTINE                       ***
C***                                                               ***
C*********************************************************************
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  ARRIVE HERE TO PREPARE POTENTIAL ARRAYS FOR AN UNSHIFTED
C  HOMONUCLEAR DIATOM, USING PREAVERAGED AA, C8, C6, ETC.
      IND = 0
      R6=RR**(-6)
      R8=R6/(RR*RR)
      EXPBR=EXP(-BETA(1)*RR)
      DO 550 J=1,NOFJ
        CALL VHOM (RR,R0,J,AA,C8,C6,BETA,R6,R8,EXPBR,VLAM,IDAMP)
        DO 540 ILEG=1,NLEG
          IND = IND+1
  540     P(IND) = VLAM(ILEG)
  550 CONTINUE
      RETURN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  ARRIVE HERE TO PREPARE POTENTIAL ARRAYS FOR A HETERONUCLEAR DIATOM BY
C  SHIFTING CENTRE OF MASS BY (DEL(J),J=1,NOFJ) ON HOMONUCLEAR POTENTIAL
C  THAT HAS ALREADY BEEN VIBRATIONALLY AVERAGED. THIS IS NOT THE BEST
C  PROCEDURE IF THE POTENTIAL IS AVAILABLE AS AN EXPLICIT FUNCTION OF
C  DIATOM BOND LENGTH.
  560 IF (NSTR.GT.1) GOTO 620
      IND = 0
      DO 610 J=1,NOFJ
        DO 570 ILEG=1,NLEG
  570     SM(ILEG) = Z0
        TT = DEL(J)/RR
        DO 590 IGP=1,NGP
          XX = XG(IGP)
          FXG = SQRT(Z1+TT*(TT+Z2*XX))
          XGP = (XX+TT)/FXG
          RRP = RR*FXG
          R6=RRP**(-6)
          R8=R6/(RRP*RRP)
          EXPBR=EXP(-BETA(1)*RRP)
          CALL VHOM (RRP,R0,J,AA,C8,C6,BETA,R6,R8,EXPBR,VLAM,IDAMP)
          YP = VLAM(1)+P2(XGP)*VLAM(2)
          DO 580 ILEG=1,NLEG
  580       SM(ILEG) = SM(ILEG)+YP*PWT(ILEG,IGP)
  590   CONTINUE
        DO 600 ILEG=1,NLEG
          IND = IND+1
          P(IND) = (ILEG-ZH)*SM(ILEG)
  600   CONTINUE
  610 CONTINUE
      RETURN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  ARRIVE HERE TO PREPARE POTENTIAL ARRAYS FOR A HETERONUCLEAR DIATOM BY
C  PERFORMING CENTRE-OF-MASS SHIFT PRIOR TO INTEGRATION OVER XI.
C  THIS IS THE PROPER PROCEDURE FOR THE CENTRE-OF-MASS SHIFT WHEN THE
C  POTENTIAL IS AVAILABLE AS AN EXPLICIT FUNCTION OF DIATOM BOND LENGTH.
C
C  THE TRANSFORMATION IS BASED ON LIU ET AL., JCP 68, 5028 (1978)
C  BUT NOTE COMMENT ABOVE ABOUT USE OF XI = (R-X0)/X0 HERE SO THAT
C  THIS IS INACCURATE FOR VIBRATIONAL STATES THAT EXTEND OUTSIDE 2*X0.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  620 DO 630 ILEG=1,NLEG
      DO 630 ISTR=1,NOFK
  630   FSM(ISTR,ILEG) = Z0
      DO 650 IGSTR=1,NGP
        XID = XG(IGSTR)
        TT = DREL*X0*(Z1+XID)/RR
      DO 650 IANG=1,NGP
        XX = XG(IANG)
        FXG = SQRT(Z1+TT*(TT+Z2*XX))
        XGP = (XX+TT)/FXG
        RRP = RR*FXG
        CALL FVHOM (RRP,XGP,XID,R0,VPRMR,BETA,VPOT,IDAMP)
        DO 640 ILEG=1,NLEG
          YP = VPOT*PWT(ILEG,IANG)
        DO 640 ISTR=1,NOFK
  640     FSM(ISTR,ILEG) = FSM(ISTR,ILEG)+YP*PWT(ISTR,IGSTR)
  650 CONTINUE
C
      DO 660 ILEG=1,NLEG
      DO 660 ISTR=1,NOFK
  660   FSM(ISTR,ILEG) = FSM(ISTR,ILEG)*(ILEG-ZH)*(ISTR-ZH)
C
      IND = 0
      DO 680 JPOT=1,NOFJ
      DO 680 ILEG=1,NLEG
        IND = IND+1
        VPOT = Z0
        DO 670 ISTR=1,NOFK
  670     VPOT = VPOT+XI(ISTR,JPOT)*FSM(ISTR,ILEG)
C
C  SCALING OPTION FOR MATRIX ELEMENTS OFF-DIAGONAL IN V
C  PURPOSE IS TO MAKE RESONANCES WIDER AND EASIER TO FIND
C
        IF (SCAL.LE.1.D0) GOTO 680
        IF (IV(JPOT).EQ.IVD(JPOT)) GOTO 680
        VPOT = VPOT*SCAL
        IF (VPOT.GT.SCALMX) VPOT = SCALMX
        IF (VPOT.LT.-SCALMX) VPOT = -SCALMX
  680   P(IND) = VPOT
      RETURN
      END

C***********************************************************************
      SUBROUTINE VHOM (RR,R0,J,AA,C8,C6,BETA,RM6,RM8,REP,VLAM,IDAMP)
C  SUBROUTINE TO GENERATE THE LEGENDRE RADIAL STRENGTH FUNCTIONS
C  VLAM(LMB) (1/CM) AT DISTANCE RR(ANGSTROMS) FROM BC3(6,8) OR TT3(6,8)
C  FORM FOR EACH LAMBDA.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION AA(2,210),C8(2,210),C6(2,210),BETA(2),VLAM(2),D(4)
      DATA Z1/1.D0/,Z4/4.D0/
      CALL DAMP(IDAMP,RR,R0,BETA(1),D)
      DO 10 ILEG=1,2
   10   VLAM(ILEG)=AA(ILEG,J)*REP-D(4)*C8(ILEG,J)*RM8
     1                           -D(2)*C6(ILEG,J)*RM6
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE FVHOM (RRP,COSTHP,XID,R0,VPRMR,BETA,VPOT,IDAMP)
C  SUBROUTINE TO GENERATE THE FULL HOMONUCLEAR ATOM-DIATOM POTENTIAL
C  VPOT  AT DISTANCE  RRP , ANGLE DEFINED BY  COS(THETA)=COSTHP  AND
C  DIATOM STRETCH OF  XID , FROM THE LE ROY - CARLEY BC3(6,8) OR
C  LE ROY - HUTSON TT3(6,8) POTENTIAL SURFACE.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION VPRMR(4,2,4), BETA(2), PSM(3), PLEG(2), D(4)
      DATA Z0/0.D0/,Z1/1.D0/,Z4/4.D0/
C
      CALL DAMP(IDAMP,RRP,R0,BETA(1),D)
C
      RM2 = Z1/RRP**2
      DO 20 II=1,3
   20   PSM(II) = Z0
      PLEG(1) = Z1
      PLEG(2) = 1.5D0*COSTHP**2-0.5D0
      XIDP = Z1
      DO 40 KK=1,4
      DO 30 LL=1,2
        FCT = PLEG(LL)*XIDP
      DO 30 II=1,3
   30   PSM(II) = PSM(II)+FCT*VPRMR(II,LL,KK)
        XIDP = XIDP*XID
   40 CONTINUE
      VPOT = PSM(1)*EXP(-BETA(1)*RRP)-
     1       (D(4)*PSM(2)*RM2+D(2)*PSM(3))*RM2**3
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE DAMP(IDAMP,R,RM,BETA,D)
C
C  SUBROUTINE TO RETURN THE DAMPING FUNCTIONS AND THEIR
C  DERIVATIVES FOR THE C6 AND C8 TERMS.
C
C  D(1:4) ARE F'(6)/F(6), F(6), F'(8)/F(8) AND F(8) RESPECTIVELY
C
C  IDAMP=1 FOR THE BUCKINGHAM-CORNER DAMPING FUNCTION
C  IDAMP=2 FOR THE TANG-TOENNIES (1984) DAMPING FUNCTION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION D(4)
C
      GOTO(100,200) IDAMP
  100 D(1)=0.D0
      D(2)=1.D0
      D(3)=0.D0
      D(4)=1.D0
      IF (R.GE.RM) RETURN
      X=RM/R
      D(1)=(12.D0/R)*X*(X-1.D0)**2
      D(2)=EXP(-4.D0*(X-1.D0)**3)
      D(3)=D(1)
      D(4)=D(2)
      RETURN
C
  200 KMAX=8
      Y=1.D0
      Z=Y
      BR=BETA*R
      DO 250 K=1,KMAX
        Y=Y*BR/DBLE(K)
        Z=Z+Y
  250   IF (K.GE.5) D(K-4)=Z
C
      Y=EXP(-BR)
      D(1)=BETA*Y*(D(2)-D(1))
      D(2)=1.D0-Y*D(2)
      D(1)=D(1)/D(2)
      D(3)=BETA*Y*(D(4)-D(3))
      D(4)=1.D0-Y*D(4)
      D(3)=D(3)/D(4)
      RETURN
      END
