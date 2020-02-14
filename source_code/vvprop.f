      SUBROUTINE VVPROP(N,NSQ,MXLAM,NHAM,
     1                  RMAT,VECNEW,W,VL,IVL,EINT,CENT,P,
     2                  A1,A1P,B1,B1P,WKS,
     3                  G1,G1P,G2,G2P,COSX,SINX,
     4                  SINE,DIAG,XK,XSQ,TSTORE,W0,
     5                  W1,W2,EYE11,EYE12,EYE22,VECOLD,
     6                  RSTART,RSTOP,NSTEP,DRNOW,DRMAX,TLDIAG,TOFF,
     7                  ERED,EP2RU,CM2RU,RSCALE,IPRINT)
C  This subroutine is part of the MOLSCAT, BOUND and FIELD suite of programs
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C------------------------------------------------------------------
C  MODIFIED FROM NRCC CODE FOR COMPATIBILITY WITH MOLSCAT
C    BY S. GREEN (FEB. 1981) AND J.M. HUTSON (OCT. 1984)
C    APR 87 MODIFY WARNING OUTPUT ASSOC. W/ 1800 FORMAT
C------------------------------------------------------------------
C  ROUTINES USED
C  WAVMAT  -CALCULATES THE POTENTIAL ENERGY INTERACTION MATRIX
C  DERMAT  -CALCULATES THE FIRST AND SECOND DERIVATIVES OF THE POTENTIAL
C  TRNSFM  -TRANSFORMS MATRICES INTO THE NEW BASIS VIA A
C               SIMILARITY TRANSFORMATION
C  PERT1   -CALCULATES THE PERTURBATION CORRECTIONS TO THE
C  PERT2        WAVEFUNCTONS.
C  DGESV   -SOLVES A LINEAR SYSTEMS OF EQUATIONS.
C  DELRD   -PREDICTS THE NEW STEP SIZE.
C  DIAGVC  -DIAGONALIZES A REAL SYMMETRIC MATRIX AND RETURN THE
C              EIGENVALUES AND EIGENVECTORS.
C------------------------------------------------------------------
C  ON ENTERING
C  N      - NUMBER OF CHANNELS
C  NSQ    - N*N
C  DRNOW  - INITIAL STEP SIZE
C  RSTART - MINIMUM RADIAL DISTANCE
C  RSTOP  - MAXIMUM RADIAL DISTANCE
C  DRMAX  - MAXIMUM ALLOWED STEP SIZE
C  TLDIAG - STEP TOLERANCE PARAMETER
C  TOFF   - INTERVAL TOLERANCE PARAMETER
C  ISCRU  - SCRATCH UNIT USED IF IREAD/IWRITE IS TRUE
C------------------------------------------------------------------
C
C  COMMON BLOCK FOR CONTROL OF USE OF PROPAGATION SCRATCH FILE
      LOGICAL IREAD,IWRITE
      COMMON /PRPSCR/ ESHIFT,ISCRU,IREAD,IWRITE

C  CHARACTER VARIABLES
C------------------------------------------------------------------
      CHARACTER(4) LRMAT,LUDP,LUD,LDG2P,LDG2,LDG1P,LDG1,LG2P,LG2,LG1P,
     1             LG1,LW0,LW2,LVECNW,LDIAG,LW1,LEYE11,LEYE12,LEYE22
C-------------------------------------------------------------------
C  LOGICAL VARIABLES
C-------------------------------------------------------------------
      LOGICAL IVD,IVPD,IVPPD,IALFP
      LOGICAL IVECT,IPOTL,IEYE,IGZRO,IGPERT,IWAVE,IRMAT,ITHS
      LOGICAL ITRUE,IFALSE,NEWINT
      LOGICAL IV,IVP,IVPP,ISHIFT,IDIAG,ICRMAT
      LOGICAL IPERT,LAST,ISYM
C-------------------------------------------------------------------
C  LABELLED COMMONS
C      CONTROL VARIABLES PASSED FROM DRIVER
C-------------------------------------------------------------------
      COMMON /LDVVCM/ XSQMAX,ALPHA1,ALPHA2,IALPHA,IALFP,IV,IVP,IVPP,
     1                ISHIFT,IDIAG,IPERT,ISYM
C-------------------------------------------------------------------
C  IF THE LOGICAL VARIABLE IS TRUE THEN
C  IV    - CALCULATES PERTURBATION CORRECTIONS FROM THE CONSTANT
C          TERMS IN THE INTERACTION POTENTIAL.
C  IVP   - CALCULATES PERTURBATION CORRECTIONS FROM THE FIRST
C          DERIVATIVE OF THE INTERACTION POTENTIAL.
C  IVPP  - CALCULATES PERTURBATION CORRECTIONS FROM THE SECOND
C          DERIVATIVE OF THE INTERACTION POTENTIAL.
C  ISHIFT- SHIFTS THE REFERENCE POTENTIAL TO BEST FIT THE TRUE
C          POTENTIAL.
C  NUMDER- CALCULATES POTENTIAL DERIVATIVES NUMERICALLY (PASSED
C          DIRECTLY VIA COMMON BLOCK)
C  IDIAG - INCLUDES ALL OF THE DIAGONAL PERTUBATION CORRECTIONS.
C  ISYM  - SYMMETRIZES THE R-MATRIX AT EACH INTERVAL.
C  IPERT - USES THE PERTURBATIONS CORRECTIONS.
C  IALFP - THE GEOMETRIC PROGRESSION PARAMETER ALPHA IS PREDICTED.
C  ALPHA1- MINIMUM GEOMETRIC PROGRESSION PARAMETER.
C  ALPHA2- MAXIMUM GEOMETRIC PROGRESSION PARAMETER.
C  IALPHA- IF IALPHA.GT.0 THEN THE STEP SIZE IS DETERMINED USING
C          A GEOMETRIC PROGRESSION AND THE INTERVAL IS DIVIDED
C          INTO IALPHA STEPS.
C------------------------------------------------------------------
      COMMON /POPT  / IVECT,IPOTL,IEYE,IGZRO,IGPERT,IWAVE,IRMAT,IOC
      LOGICAL LPOPT(7)
      EQUIVALENCE (LPOPT(1),IVECT)
C  LPOPT CONTAINS PRINTING OPTIONS FROM NRCC VERSION.
C  THESE ARE ALL SET FALSE HERE.  CHANGE TO DEBUG.
C  WHEN THE LOGICAL VARIABLE IS TRUE,
C  IVECT - EIGENVALUES AND EIGENVECTORS.
C  IPOTL - POTENTIAL ENERGY MATRICES AND ITS DERIVATIVES.
C  IEYE  - ACCUMULATED PERTURBATION INTEGRALS.
C  IGZRO - ZERO-TH ORDER WAVEFUNCTIONS.
C  IGPERT- PERTURBED WAVEFUNCTIONS.
C  IWAVE - PERTURBED WAVEFUNCTIONS.
C  IRMAT - R-MATRIX
C  IOC   - INFORMATION PRINTED EVERY IOC-TH STEP
C--------------------------------------------------------------------
C  ARRAYS DIMENSIONED AS VECTORS
C-------------------------------------------------------------------
      DIMENSION G1(N),G1P(N),G2(N),G2P(N)
      DIMENSION A1(N),A1P(N),B1(N),B1P(N)
      DIMENSION XSQ(N),XK(N),COSX(N),SINX(N),SINE(N),DIAG(N)
C-------------------------------------------------------------------
C  ARRAYS DIMENSIONED AS MATRICES
C-------------------------------------------------------------------
      DIMENSION EYE11(NSQ),EYE12(NSQ),EYE22(NSQ)
      DIMENSION W0(NSQ),W1(NSQ),W2(NSQ),W(NSQ)
      DIMENSION RMAT(NSQ)
      DIMENSION TSTORE(NSQ),VECOLD(NSQ),VECNEW(NSQ)
      DIMENSION P(MXLAM),VL(2),IVL(2),EINT(N),CENT(N),WKS(N)
C-------------------------------------------------------------------
C  DATA STATEMENTS FOR PRINTING
C-------------------------------------------------------------------
      DATA LRMAT/'RMAT'/,LUDP/'  UP'/,LUD/'   U'/,LDG2P/'DG2P'/
      DATA LDG2/' DG2'/,LDG1P/'DG1P'/,LDG1/' DG1'/,LG2P/' G2P'/
      DATA LG2/'  G2'/,LG1P/' G1P'/,LG1/'  G1'/,LW0/'  W0'/
      DATA LW2/'  W2'/,LVECNW/'VCNW'/,LDIAG/'DIAG'/
      DATA LW1/'  W1'/,LEYE11/' I11'/,LEYE12/' I12'/,LEYE22/' I22'/
C-------------------------------------------------------------------
C  LOGICAL DATA STATEMENTS
C-------------------------------------------------------------------
      DATA IFALSE/.FALSE./,ITRUE/.TRUE./
C-------------------------------------------------------------------
C
C  SET DEFAULT VALUES FOR PRINTING
      NSGERR=0
      IOC=5
      SOFF = 0.D0
      COFF = 0.D0
      CDIAG = 0.D0
      SDIAG = 0.D0
      DO 100 I=1,7
  100   LPOPT(I)=.FALSE.
      IF (.NOT.(IREAD .AND. IWRITE)) GOTO 101

      WRITE(6,699)
  699 FORMAT(/'  * * * ERROR.  IREAD AND IWRITE CANNOT BOTH BE TRUE.')
      IREAD=.FALSE.
      IWRITE=.FALSE.
  101 CONTINUE
      IVD = IV .OR. IDIAG
      IVPD = IVP .OR. IDIAG
      IVPPD = IVPP .OR. IDIAG
C-------------------------------------------------------------------
C  PRINT CONTROL DATA
      IF (IPRINT.LT.20) GOTO 110

      WRITE(6,1200)
      WRITE(6,1300) IVECT,IPOTL,IEYE,IGZRO,IGPERT,IWAVE,IRMAT,
     X              IWRITE,IREAD,IOC
      WRITE(6,1400)
      WRITE(6,1500) IV,IVP,IVPP,ISHIFT,IDIAG,ISYM,IPERT,IALFP
      WRITE(6,1600) ALPHA1,ALPHA2,IALPHA
      WRITE(6,2500)
      WRITE(6,2600) RSTART,RSTOP,DRNOW,DRMAX,TOFF,TLDIAG
  110 COFFL = 0.D0
      IF (N.EQ.0) RETURN

      NEWINT = .FALSE.
      NP1 = N+1
      ICRMAT = .TRUE.
      TOL = 1.D-11
      LAST = .FALSE.
      ITRANS = 0
      IK = 0
      DO 130 I = 1,N
        G1(I) = 0.D0
        G1P(I) = 1.D0
        G2(I) = 1.D0
        G2P(I) = 0.D0
        DIAG(I)=0.D0
      DO 130 K = 1,N
        IK = IK+1
        VECNEW(IK) = 0.D0
        IF (I.EQ.K) VECNEW(IK) = 1.D0
        EYE11(IK) = 0.D0
        EYE12(IK) = 0.D0
  130   EYE22(IK) = 0.D0
      IF (IPRINT.GE.20) WRITE(6,3100)
      ISTEP = 1
      NTRVL = 0
      DRINT = DRNOW
      DIAGI = RSTART+0.5D0*DRINT
      RMID = RSTART
      RLAST = RSTART
      XBAR = 0.D0
      XSBAR = 0.D0
      EBAR = 0.D0
      EXBAR = 0.D0
      IF (IALPHA.LE.0) GOTO 150

      BALPHA = (ALPHA2-ALPHA1)/(RSTOP-RSTART)
      ALPHA = ALPHA1+BALPHA*(DIAGI-RSTART)
      IF (IALFP) ALPHA = ALPHA1
      IF (ALPHA.NE.1.D0) GOTO 140

      DRNOW = DRINT/IALPHA
      GOTO 150

  140 DRNOW = DRINT*(ALPHA-1.D0)/(ALPHA**IALPHA-1.D0)
  150 RNOW = RSTART+DRNOW
      IF (IWRITE) WRITE(ISCRU) RSTART,RSTOP,DRINT,DIAGI,RMID,IALPHA,
     X                         BALPHA,ALPHA1,ALPHA2,IALPHA
      IF (.NOT. IREAD) GOTO 160
      READ(ISCRU) RSTART,RSTOP,DRINT,DIAGI,RMID,IALPHA,BALPHA,
     X            ALPHA1,ALPHA2,IALPHA
C-------------------------------------------------------------------
C  START OF THE PROPAGATION LOOP
C-------------------------------------------------------------------
  155 READ(ISCRU) ISTEP,RNOW,DRNOW,LAST,N,DIAG,TSTORE,
     X            W0,W1,W2,VECNEW,NTRVL,RMIDI,RLAST,RCENT,ITRANS
      READ(ISCRU) NEWINT
      DO 158 I=1,N
  158   DIAG(I)=DIAG(I)+ESHIFT
  160 RCENT = RNOW-0.5D0*DRNOW
      ITHS = .FALSE.
      IF (((NTRVL+1)/IOC)*IOC.EQ.NTRVL+1) ITHS = .TRUE.
      IF (IREAD) GOTO 300

C-------------------------------------------------------------------
C  EVALUATE THE POTENTIAL AND ITS DERIVATIVES.
C-------------------------------------------------------------------
      CALL WAVMAT(W,N,RCENT,P,VL,IVL,ERED,EINT,CENT,EP2RU,CM2RU,
     1            RSCALE,WKS,MXLAM,NHAM,IPRINT)
      DO 165 I = 1, NSQ
  165   W0(I) = W(I)
      IF (IVPD .AND. IVPPD) GOTO 200

      DO 170 I = 1, NSQ
        W1(I) = 0.D0
  170   W2(I) = 0.D0
  200 IF (IVPPD .OR. ISHIFT) CALL DERMAT(2,W2,N,RCENT,P,VL,IVL,CENT,
     1                                   EP2RU,RSCALE,MXLAM,NHAM,
     2                                   IPRINT)
      IF (IVPD) CALL DERMAT(1,W1,N,RCENT,P,VL,IVL,CENT,
     1                      EP2RU,RSCALE,MXLAM,NHAM,IPRINT)
      FACTOR = DRNOW*DRNOW/24.D0
      IF (.NOT.ISHIFT) FACTOR = 0.D0
      IF (.NOT.ICRMAT) GOTO 270

      RMIDI = DIAGI
C-------------------------------------------------------------------
C  EVALUATE THE POTENTIAL AT THE RMIDI WHERE THE INTERACTION IS TO
C  BE DIAGONALIZED AND SAVE THE OLD EIGENVECTORS.
C-------------------------------------------------------------------
      IF (RMIDI.NE.RCENT)
     1  CALL WAVMAT(W,N,RMIDI,P,VL,IVL,ERED,EINT,CENT,EP2RU,CM2RU,
     2              RSCALE,WKS,MXLAM,NHAM,IPRINT)
      DO 240 I = 1,NSQ
  240   VECOLD(I) = VECNEW(I)
      ITRANS = ITRANS+1
C-------------------------------------------------------------------
C  DIAGONALIZE THE INTERACTION POTENTIAL.
C-------------------------------------------------------------------
      IFAIL=0
      CALL DIAGVC(W,N,N,DIAG,VECNEW)
      IF (.NOT.ITHS) GOTO 270
      IF (.NOT.IVECT) GOTO 270

      WRITE(6,2900) LDIAG
      WRITE(6,2800) (DIAG(I),I = 1,N)
      WRITE(6,2900) LVECNW
      WRITE(6,2800) (VECNEW(I),I = 1,NSQ)
C--------------------------------------------------------------------
C  TRANSFORM THE POTENTIAL AND ITS DERIVATIVES INTO THE LOCAL BASIS.
C--------------------------------------------------------------------
  270 CALL TRNSFM(VECNEW,W0,TSTORE,N,IFALSE,ITRUE)
      IF (IVPD) CALL TRNSFM(VECNEW,W1,TSTORE,N,IFALSE,ITRUE)
      IF (IVPPD .OR. ISHIFT) CALL TRNSFM(VECNEW,W2,TSTORE,N,IFALSE,
     X                                   ITRUE)
C-------------------------------------------------------------------
C  DETERMINE THE NEW TRANSFORMATION MATRIX
C-------------------------------------------------------------------
      IF (ICRMAT) CALL DGEMUL(VECOLD,N,'T',VECNEW,N,'N',TSTORE,N,
     1                        N,N,N)
C-------------------------------------------------------------------
C  TRANSFORM THE R-MATRIX INTO THE NEW BASIS.
C-------------------------------------------------------------------
  300 IF (ICRMAT) CALL TRNSFM(TSTORE,RMAT,W,N,IFALSE,ISYM)
      ICRMAT = .FALSE.
      IF (IREAD) GOTO 350

C-------------------------------------------------------------------
C  SHIFT THE EIGENVALUES AND INITIALIZE FOR CONTRIBUTIONS NOT DESIRED.
C-------------------------------------------------------------------
      INDX = -N
      DO 330 J = 1,N
        INDX = INDX+NP1
        DIAG(J) = -W0(INDX)-FACTOR*W2(INDX)
  330   W0(INDX) = -FACTOR*W2(INDX)
      IF (IVD .AND. IVPD .AND. IVPPD) GOTO 350

      CTERM0 = 1.D0
      CTERM1 = 1.D0
      CTERM2 = 1.D0
      IF (.NOT.IVD) CTERM0 = 0.D0
      IF (.NOT.IVPD) CTERM1 = 0.D0
      IF (.NOT.IVPPD) CTERM2 = 0.D0
      DO 340 I = 1,NSQ
        W0(I) = W0(I)*CTERM0
        W1(I) = W1(I)*CTERM1
  340   W2(I) = W2(I)*CTERM2
C-------------------------------------------------------------------
C  WRITE ON UNIT ISCRU THE INFORMATION NECESSARY FOR SUBSEQUENT ENERGY
C  CALCULATIONS.
C-------------------------------------------------------------------
  350 IF (IWRITE) WRITE(ISCRU) ISTEP,RNOW,DRNOW,LAST,N,DIAG,TSTORE,
     1                         W0,W1,W2,VECNEW,NTRVL,RMIDI,RLAST,RCENT,
     2                         ITRANS
      IF (.NOT.ITHS) GOTO 360
      IF (.NOT.IPOTL) GOTO 360

      WRITE(6,2900) LDIAG
      WRITE(6,2800) (DIAG(I),I = 1,N)
      WRITE(6,2900) LW0
      WRITE(6,2800) (W0(I),I = 1,NSQ)
      WRITE(6,2900) LW1
      WRITE(6,2800) (W1(I),I = 1,NSQ)
      WRITE(6,2900) LW2
      WRITE(6,2800) (W2(I),I = 1,NSQ)
      WRITE(6,2900) IALPHA
C-------------------------------------------------------------------
C  CALCULATE THE ZERO-TH ORDER WAVEFUNCTIONS AND DERIVATIVES.
C-------------------------------------------------------------------
  360 NOPLOC = 0
      DO 390 I = 1,N
        DIF = DIAG(I)
        XSQ(I) = DIF*DRNOW*DRNOW
        XLMBDA = SQRT(ABS(DIF))
        X = XLMBDA*DRNOW
        IF (DIF.LT.0.D0) GOTO 370

        NOPLOC = NOPLOC+1
        SX = SIN(X)/XLMBDA
        CX = COS(X)
        GOTO 380

  370   IF (X.GT.173.D0) WRITE(6,1700) I,DIF,DRNOW,X
        SX = SINH(X)/XLMBDA
        CX = COSH(X)
  380   A = G1P(I)
        SINX(I) = SX
        SINE(I) = SX*XLMBDA
        IF (DIF.LT.0.D0) SINE(I) = -SINE(I)
        COSX(I) = CX
        XK(I) = X
        B = G1(I)
        G1(I) = A*SX+B*CX
        G1P(I) = A*CX-DIF*B*SX
        C = G2P(I)
        D = G2(I)
        A1(I) = B
        A1P(I) = A
        B1(I) = D
        B1P(I) = C
        G2(I) = C*SX+D*CX
  390   G2P(I) = C*CX-DIF*D*SX

C-------------------------------------------------------------------
C  ESTIMATE G2P(N) AT END OF NEXT STEP. IF IT IS TOO LARGE,
C  A NEW INTERVAL WILL BE STARTED.
C-------------------------------------------------------------------
      IF (ABS(G2P(N)).LE.1.D4) GOTO 1801

      IF (IPRINT.GE.10) WRITE(6,1800) RNOW,DRNOW,G2P(N)
      NSGERR=NSGERR+1
 1801 G2PMAX=G2P(N)*CX
      IF (IREAD .AND. .NOT.IPERT) GOTO 410

C-------------------------------------------------------------------
C  CALCULATE THE INTEGRALS NECESSARY FOR THE PERTURBATION CORRECTIONS.
C  THE STEP INTEGRALS ARE STORED IN W0, W1 AND W2, AND THE ACCUMULATED
C  INTEGRALS OVER THE INTERVAL ARE SAVED IN EYE11, EYE12 AND EYE22.
C-------------------------------------------------------------------
      IF (IVPP) CALL PERT2(N,COSX,SINE,XSQ,XK,DRNOW,W0,W1,W2,EYE11,
     X                     EYE12,EYE22,A1,B1,A1P,B1P)
      IF (IVPP) GOTO 400

      IF (IVP) CALL PERT1(N,COSX,SINE,XSQ,XK,DRNOW,W0,W1,W2,EYE11,
     X                    EYE12,EYE22,A1,B1,A1P,B1P)
      IF (IVP) GOTO 400

      IF (IV) CALL PERT2(N,COSX,SINE,XSQ,XK,DRNOW,W0,W1,W2,EYE11,
     X                   EYE12,EYE22,A1,B1,A1P,B1P)
  400 CONTINUE
  410 IF (IREAD .AND. .NOT.NEWINT) GOTO 590

      SOFF = 0.D0
      COFF = 0.D0
      CDIAG = 0.D0
      SDIAG = 0.D0
C-------------------------------------------------------------------
C  THE FOLLOWING IS USED TO DETERMINE THE MAXIMUM PERTURBATION
C  CORRECTIONS TO THE UNPERTURBED WAVEFUNCTIONS. SINCE THE STEP
C  SIZE FOR SUBSEQUENT ENERGIES HAS ALREADY BEEN STORED ON DISK
C  THIS INFORMATION IS NOT NECESSARY FOR SUBSEQUENT ENERGIES.
C-------------------------------------------------------------------
      IF (IREAD .AND. .NOT.IPERT) GOTO 460

      IF (IREAD) GOTO 430

      DO 420 I = 1,N
        A1(I) = 1.D0/SQRT(A1P(I)*A1P(I)/ABS(DIAG(I))+A1(I)*A1(I))
        B1(I) = 1.D0/SQRT(B1P(I)*B1P(I)/ABS(DIAG(I))+B1(I)*B1(I))
        A1P(I) = DRNOW*A1(I)/XK(I)
        B1P(I) = DRNOW*B1(I)/XK(I)
        SINE(I) = 1.D0
        IF (DIAG(I).GT.0.D0) GOTO 420

        EXPX = EXP(-XK(I)*DRINT/DRNOW)
        IF (DIAG(I).LT.-XSQMAX) EXPX=0.D0
        SINE(I) = EXPX
        A1(I) = A1(I)*EXPX
        B1(I) = B1(I)*EXPX
        A1P(I) = A1P(I)*EXPX
        B1P(I) = B1P(I)*EXPX
  420 CONTINUE
  430 IJ = 0

C-------------------------------------------------------------------
C  CALCULATE THE PERTURBATION CORRECTIONS TO THE WAVEFUNCTION AND
C  ITS DERIVATIVE.
C-------------------------------------------------------------------
      DO 450 J = 1,N
        A1J = A1(J)
        A1PJ = A1P(J)
      DO 450 I = 1,N
        JI = J+(I-1)*N
        IJ = IJ+1
        PRT1 = G1(J)*EYE12(IJ)-G2(J)*EYE11(IJ)
        PRT2 = G1(I)*EYE22(IJ)-G2(I)*EYE12(IJ)
        PRT1P = G1P(J)*EYE12(IJ)-G2P(J)*EYE11(IJ)
        PRT2P = G1P(I)*EYE22(IJ)-G2P(I)*EYE12(IJ)
C-------------------------------------------------------------------
C  DON'T DETERMINE THE MAXIMUM PERTURBATION CORRECTION FOR
C  SUBSEQUENT ENERGIES.
C-------------------------------------------------------------------
        IF (IREAD) GOTO 440

        B1I = B1(I)
        B1PI = B1P(I)
        E1 = ABS(PRT1)*A1J
        E2 = ABS(PRT2)*B1I
        E3 = ABS(PRT1P)*A1PJ
        E4 = ABS(PRT2P)*B1PI
        IF (I.NE.J) COFF = MAX(COFF,E1,E2,E3,E4)
        IF (I.EQ.J) CDIAG = MAX(CDIAG,E1,E2,E3,E4)
        IF (J.GT.I) GOTO 440

        CCIJ = W0(IJ)
        CCJI = W0(JI)
        CSIJ = W1(IJ)
        CSJI = W1(JI)
        SSIJ = W2(IJ)
        SSJI = W2(JI)
        E1 = ABS(SINX(J)*CSJI-COSX(J)*SSIJ)*SINE(J)*XK(J)/DRNOW
        E2 = ABS(SINX(I)*CCIJ-COSX(I)*CSJI)*SINE(I)
        E3 = ABS(COSX(J)*CSJI+DIAG(J)*SINX(J)*SSIJ)*SINE(J)
        E4 = ABS(COSX(I)*CCIJ+DIAG(I)*SINX(I)*CSJI)*SINE(I)*DRNOW/XK(I)
        E5 = ABS(SINX(I)*CSIJ-COSX(I)*SSJI)*SINE(I)*XK(I)/DRNOW
        E6 = ABS(SINX(J)*CCJI-COSX(J)*CSIJ)*SINE(J)
        E7 = ABS(COSX(I)*CSIJ+DIAG(I)*SINX(I)*SSJI)*SINE(I)
        E8 = ABS(COSX(J)*CCJI+DIAG(J)*SINX(J)*CSIJ)*SINE(J)*DRNOW/XK(J)
        IF (I.NE.J) SOFF = MAX(SOFF,E1,E2,E3,E4,E5,E6,E7,E8)
        IF (I.EQ.J) SDIAG = MAX(SDIAG,E1,E2,E3,E4,E5,E6,E7,E8)
  440   W2(IJ) = PRT1
        W(JI) = PRT2
        W0(IJ) = PRT1P
  450   W1(JI) = PRT2P

      IF (SOFF.EQ.0.D0) SOFF=1.D-30
      IF (IPERT) GOTO 480

  460 DO 470 I = 1,NSQ
        W2(I) = 0.D0
        W0(I) = 0.D0
        W(I) = 0.D0
  470   W1(I) = 0.D0
  480 IF (LAST) GOTO 500

      IF (IALPHA.LE.0) GOTO 485

      IF ((ISTEP/IALPHA)*IALPHA.EQ.ISTEP) GOTO 500

      GOTO 590

C-------------------------------------------------------------------
C  ARRIVE HERE ONLY FOR IALPHA.EQ.0 OPTION.
C  START NEW INTERVAL IF PREDICTED G2P FOR NEXT STEP IS TOO LARGE
C-------------------------------------------------------------------
  485 IF (.NOT.IREAD .AND. ABS(G2PMAX).GT.1.D4) GOTO 500

      IF (COFFL.EQ.0.D0) GOTO 490

      FACC = COFF/COFFL
      FACS = SOFF/SOFF1
      IF (FACC.GT.2.D0) FACC = 2.D0
      IF (FACS.GT.2.D0) FACS = 2.D0
      IF (FACC*COFF.GT.0.8D0*TOFF) GOTO 500
      IF (FACS*SOFF.GT.0.8D0*TLDIAG) GOTO 500

  490 COFFL = COFF
      SOFF1 = SOFF
      SDIAG1 = SDIAG
      COFFL = COFF
      IF (IREAD .AND. NEWINT) GOTO 500

C-------------------------------------------------------------------
C  CHECK TO SEE IF THE PERTURBATION CORRECTIONS ARE LARGE ENOUGH
C  TO WARRANT A NEW INTERVAL AND BASIS SET TRANSFORMATION.
C-------------------------------------------------------------------
      IF (COFF.LT.0.8D0*TOFF .AND. CDIAG.LT.0.8D0*TOFF) GOTO 590

  500 COFFL = 0.D0
      SOFF1 = SOFF
      SDIAG1 = SDIAG
      ICRMAT = .TRUE.
      IF (.NOT.ITHS) GOTO 510
      IF (.NOT.IEYE) GOTO 510

      WRITE(6,2900) LEYE11
      WRITE(6,2800) (EYE11(I),I = 1,NSQ)
      WRITE(6,2900) LEYE12
      WRITE(6,2800) (EYE12(I),I = 1,NSQ)
      WRITE(6,2900) LEYE22
      WRITE(6,2800) (EYE22(I),I = 1,NSQ)
C-------------------------------------------------------------------
C  MULTIPLY THE OLD R-MATRIX TIMES IRREGULAR WAVEFUNCTION AND ITS
C  PERTURBATION CORRECTION.
C-------------------------------------------------------------------
  510 NP1 = N+1
      II = 1
      DO 520 I = 1,N
        W0(II) = W0(II)+G1P(I)
        W1(II) = W1(II)+G2P(I)
        W2(II) = W2(II)+G1(I)
        W(II) = W(II)+G2(I)
  520   II = II+NP1
      CALL DCOPY(NSQ,W0,1,TSTORE,1)
      CALL DSYMM('L','L',N,N,1.D0,RMAT,N,W1,N,1.D0,TSTORE,N)
      CALL DCOPY(NSQ,W2,1,VECOLD,1)
      CALL DSYMM('L','L',N,N,1.D0,RMAT,N,W ,N,1.D0,VECOLD,N)
      IF (.NOT.ITHS) GOTO 550

      IF (.NOT.IGZRO) GOTO 540

      WRITE(6,2900) LG1
      WRITE(6,2800) (G1(I),I = 1,N)
      WRITE(6,2900) LG1P
      WRITE(6,2800) (G1P(I),I = 1,N)
      WRITE(6,2900) LG2
      WRITE(6,2800) (G2(I),I = 1,N)
      WRITE(6,2900) LG2P
      WRITE(6,2800) (G2P(I),I = 1,N)
  540 IF (.NOT.IGPERT) GOTO 550

      WRITE(6,2900) LDG1
      WRITE(6,2800) (W2(I),I = 1,NSQ)
      WRITE(6,2900) LDG1P
      WRITE(6,2800) (W0(I),I = 1,NSQ)
      WRITE(6,2900) LDG2
      WRITE(6,2800) (W(I),I = 1,NSQ)
      WRITE(6,2900) LDG2P
      WRITE(6,2800) (W1(I),I = 1,NSQ)
  550 IF (.NOT.ITHS) GOTO 560
      IF (.NOT.IWAVE) GOTO 560

      WRITE(6,2900) LUD
      WRITE(6,2800) (EYE12(I),I = 1,NSQ)
      WRITE(6,2900) LUDP
      WRITE(6,2800) (EYE11(I),I = 1,NSQ)
  560 IER = 0
C-------------------------------------------------------------------
C  SOLVE A LINEAR SYSTEM OF EQUATIONS TO DETERMINE THE NEW R-MATRIX
C-------------------------------------------------------------------
      CALL DGESV(N,N,TSTORE,N,WKS,VECOLD,N,IER)
C-------------------------------------------------------------------
C  REINITIALIZE FOR THE NEXT INTERVAL. STORE THE NEW R-MATRIX IN RMAT.
C-------------------------------------------------------------------
      DO 570 I = 1,N
        G1(I) = 0.D0
        G2(I) = 1.D0
        G1P(I) = 1.D0
        G2P(I) = 0.D0
      DO 570 J = 1,N
        IJ = I+(J-1)*N
        JI = J+(I-1)*N
        RMAT(JI) = VECOLD(IJ)
        IF (ISYM) RMAT(JI) = 0.5D0*(VECOLD(IJ)+VECOLD(JI))
        EYE11(IJ) = 0.D0
        EYE12(IJ) = 0.D0
  570   EYE22(IJ) = 0.D0

      NTRVL = NTRVL+1
      IF (.NOT.ITHS) GOTO 580
      IF (.NOT.IRMAT) GOTO 580

      WRITE(6,2900) LRMAT
      WRITE(6,2800) (RMAT(I),I = 1,NSQ)
      GOTO 590

  580 CONTINUE
  590 CONTINUE
      IF (.NOT.IWRITE) GOTO 600

      NEWINT=ICRMAT
      WRITE(ISCRU) NEWINT
C-------------------------------------------------------------------
C  WRITE THE MINIMAL STEP INFORMATION AND DETERMINE THE NEW STEP SIZE
C  AND PREDICT THE NEW INTERVAL SIZE.
C-------------------------------------------------------------------
  600 IF (IPRINT.GE.20) WRITE(6,2700) ITRANS,RMIDI,DRNOW,RLAST,RNOW,
     X                                DIAG(1),DIAG(N),CDIAG,COFF,SDIAG,
     &                                SOFF,ALPHA,NOPLOC,ISTEP
      IF (.NOT.LAST) ISTEP = ISTEP+1
      IF (LAST) GOTO 670

      IF (IREAD) GOTO 155

      XBAR = XBAR+RNOW
      XSBAR = XSBAR+RNOW*RNOW
      EBAR = EBAR+SDIAG
      EXBAR = EXBAR+RNOW*SDIAG
      TMXX = 0.5D0*TOFF
      IF (TLDIAG.GT.TMXX) TMXX = TLDIAG
      IF (DRINT.NE.DRNOW) SOFF = 0.8D0*TLDIAG*(SOFF/(0.8D0*TMXX))**1.5D0
      IF (IALPHA.GT.0) DRNOW = DRNOW*ALPHA
      IF (IALPHA.LE.0) CALL DELRD(DRNOW,SDIAG,SOFF,TLDIAG,DRMAX,
     1                            DIAG(1),DIAG(N),RNOW,RSTOP)
      IF (.NOT.ICRMAT) GOTO 650

      DRINT = RNOW-RMID
      DRINT1 = DRINT
      IF (IALPHA.LE.0) GOTO 630

      XBAR = XBAR/IALPHA
      XSBAR = XSBAR/IALPHA
      EBAR = EBAR/IALPHA
      EXBAR = EXBAR/IALPHA
      IF (IALPHA.EQ.1) SLOPE=0.D0
      IF (IALPHA.NE.1) SLOPE = (EXBAR-XBAR*EBAR)/(XSBAR-XBAR*XBAR)
      BINT = EBAR-XBAR*SLOPE
      EMAX = BINT+SLOPE*RNOW
      EMIN = BINT+SLOPE*(RNOW-DRINT)
      ALFNEW = ALPHA
      IF (IALPHA.LE.1) GOTO 630

      IF (EMAX.EQ.0.D0) EMAX=1.D-30
      FAC = EMIN/EMAX
      IF (FAC.LE.0.D0) GOTO 620

      FAC = (FAC)**(1.D0/DBLE(3*IALPHA-3))
      IF (FAC.GT.1.1D0) FAC = 1.1D0
      IF (FAC.LT.0.9D0) FAC = 0.9D0
      ALFNEW = ALPHA*FAC
      GOTO 630

  620 FAC = 1.1D0
      IF (EMIN.LE.0.D0) FAC = 0.9D0
      ALFNEW = ALPHA*FAC
  630 XBAR = 0.D0
      XSBAR = 0.D0
      EBAR = 0.D0
      EXBAR = 0.D0
      TMXX = TOFF
      IF (TLDIAG.GT.TOFF) TMXX = TLDIAG
      IF (DRINT.NE.DRNOW) COFF = 0.8D0*TMXX*(COFF/(0.8D0*TMXX))**1.5D0
      CALL DELRD(DRINT,CDIAG,COFF,TMXX,DRMAX,
     1           DIAG(1),DIAG(N),RNOW,RSTOP)

      IF (DRINT1.NE.DRNOW) SOFF1 = TLDIAG*(2.D0*SOFF1/TLDIAG)**1.5D0
      CALL DELRD(DRINT1,SDIAG1,SOFF1,TLDIAG,DRMAX,
     1           DIAG(1),DIAG(N),RNOW,RSTOP)

      IF (ABS(DRINT1).LT.ABS(DRINT)) DRINT = DRINT1
      IF (DRINT.LT.DRNOW) DRINT = DRNOW
      IF (ABS(RSTOP-RNOW-DRINT).LT.ABS(0.01D0*DRINT)) DRINT = RSTOP-RNOW
      IF ((RSTOP-RNOW-DRINT)*DRINT.LT.0.D0) DRINT = RSTOP-RNOW
      RMID = RNOW
      DIAGI = RNOW+0.5D0*DRINT
      IF (IALPHA.LE.0) GOTO 650

      ALPHA = ALPHA1+BALPHA*(DIAGI-RSTART)
      IF (IALFP) ALPHA = ALFNEW
      IF (ALPHA.NE.1.D0) GOTO 640

      DRNOW = DRINT/IALPHA
      GOTO 650

  640 DRNOW = DRINT*(ALPHA-1.D0)/(ALPHA**IALPHA-1.D0)
  650 IF (ABS(RSTOP-RNOW-DRNOW).LT.ABS(0.01D0*DRNOW)) DRNOW = RSTOP-RNOW
      IF ((RSTOP-RNOW-DRNOW)*DRNOW.LT.0.D0) DRNOW = RSTOP-RNOW
      RLAST = RNOW
      RNOW = RNOW+DRNOW
      DEL = (RNOW-RSTOP)/DRNOW
      IF (ABS(DEL).LT.0.005D0) LAST = .TRUE.
      GOTO 160

C-------------------------------------------------------------------
C  THE PROPAGATION IS NOW COMPLETE. TRANSFORM THE R-MATRIX INTO THE
C  ORIGINAL BASIS.
C-------------------------------------------------------------------
  670 NCOL = 1
      NLAST = N
      DO 690 IR = 1,N
        NORIG = IR
        DO 680 NTRANS = NCOL,NLAST
          VECOLD(NTRANS) = VECNEW(NORIG)
  680     NORIG = NORIG+N
        NLAST = NLAST+N
  690   NCOL = NCOL+N

      CALL TRNSFM(VECOLD,RMAT,TSTORE,N,IFALSE,ISYM)
      IF (IPRINT.LT.10 .AND. NSGERR.GT.0) WRITE(6,1802) NSGERR
      NSTEP=ISTEP
      RETURN

C-------------------------------------------------------------------
C  FORMAT STATEMENTS
C-------------------------------------------------------------------
 1200 FORMAT(/'  IVECT IPOTL  IEYE IGZRO IGPERT IWAVE IRMAT IWRITE',
     1       ' IREAD IOC')
 1300 FORMAT(1X, 9L6,I4)
 1400 FORMAT(/2X,'IV IVP IVPP ISHIFT IDIAG ISYM IPERT IALFP')
 1500 FORMAT(1X,L3,L4,L5,L7,L6,L5,4L6)
 1600 FORMAT(/1X,' ALPHA1 ALPHA2 IALPHA'/1X,2F7.2,I7)
 1700 FORMAT(/'  *** ERROR IN VVPROP.  FOR CHANNEL',I3,
     1       ', REDUCED V-E =',E13.5/6X,'FOR STEP SIZE',E13.5,
     2       ', COSH ARGUMENT OF',E13.5,' WILL CAUSE OVERFLOW.'/
     3       6X,'USE A SMALLER STEP SIZE TO AVOID THIS ERROR.')
 1800 FORMAT(/'  *** WARNING IN VVPROP. INTERVAL SIZE TOO LARGE, SO',
     1       ' CLOSED CHANNEL GROWTH MAY CAUSE NUMERICAL INSTABILITY.'/
     2       24X,'RNOW =',F8.3,',  DRNOW =',F8.3,', G2P(N) =',E13.5)
 1802 FORMAT(/'  *** WARNING IN VVPROP.  INTERVAL SIZE POSSIBLY TOO',
     1       ' LARGE FOR',I5,' STEPS.  INCREASE IPRINT FOR DETAILS')
 2500 FORMAT(/' RSTART    RSTOP     DRNOW     DRMAX     TOFF     ',
     1       'TLDIAG')
 2600 FORMAT(F9.5,10F10.5,I10)
 2700 FORMAT(1X,I5,10E11.4,F6.3,2I5)
 2800 FORMAT(1X,9E14.7)
 2900 FORMAT(/2X,A10)
 3000 FORMAT(/'  VVPROP.   R-MATRIX PROPAGATED FROM',F12.4,'  TO',
     &       F12.4,'  IN',I6,'  STEPS.')
 3100 FORMAT(/' NTRVL    RCENT      DRNOW      RLAST      RNOW      ',
     X       'DIAG(1)    DIAG(N)     CDIAG      COFF       SDIAG      ',
     X       'SOFF     ALFP NOPN ISTP')
C----------------***END-VVPROP***-------------------------------------
      END
