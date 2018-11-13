      SUBROUTINE SCPROP(NORIG,MXLAM,NV,NPOTL,
     1                  JSINDX,SR,SI,U,VL,
     2                  IV,EINT,CENT,WVEC,
     3                  L,NB,P,DEGTOL,NOPMAX,DEEP,
     4                  IK,ICODE,IPRINT,IBOUND,ICHAN,WAVE,WMAXIN)
C  Copyright (C) 2018 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE potential
C
C  THIS SUBROUTINE SETS UP THE STORAGE REQUIREMENTS FOR ALL THE
C  DIFFERENT PROPAGATORS IMPLEMENTED
C
C  03-12-15 CR Le Sueur:
C  LONG RANGE AND SHORT RANGE PROPAGATORS HAVE BEEN DECOUPLED.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IPRINT
      LOGICAL IREAD,IWRITE
      DIMENSION JSINDX(*),SR(*),SI(*),U(*),VL(*),IV(*),EINT(*),CENT(*),
     1          WVEC(*),L(*),NB(*),P(*)
C
C  DYNAMIC STORAGE COMMON BLOCK ...
      COMMON /MEMORY/ MX,IXNEXT,NIPR,IDUMMY,X(1)
C
      COMMON /DRIVE / STEST,STEPS,STABIL,CONV,RMIN,RMAX,XEPS,DR,
     1                DRMAX,RMID,RMATCH,TOLHI,RTURN,VTOL,ESHIFT,
     2                ERED,RMLMDA,NOPEN,JKEEP,ISCRU,MAXSTP,ILDSVU
      COMMON /HIBRIN/ POWRX,DRAIRY
      LOGICAL WAVE

C  COMMON BLOCK FOR PROPAGATOR CONTROL
      INTEGER MXSEG,IDIR,IFLG,NSEGS
      DOUBLE PRECISION RST,REND,DRPR
      PARAMETER (MXSEG=3)
      COMMON /PRPDTA/ RST(MXSEG),REND(MXSEG),DRPR(MXSEG),
     1                IDIR(MXSEG),IFLG(MXSEG),NSEGS,NSTEPS(MXSEG)

C  COMMON BLOCK FOR DERIVATIVES
      LOGICAL NUMDER
      COMMON /DERIVS/ NUMDER

C  COMMON BLOCK FOR INPUT/OUTPUT CHANNEL NUMBERS
      LOGICAL PSIFMT
      INTEGER IPSISC,IWAVSC,IPSI
      COMMON /IOCHAN/ IPSISC,IWAVSC,IPSI,NWVCOL,PSIFMT
C
      N=NORIG ! BECAUSE RMTPRP CAN CHANGE THE VALUE OF N
      NSQ=N*N
      IREC=1
      IF (WAVE) CALL WVPROP(N)
C
C  IC2 IS NEXT AVAILABLE LOCATION ...
      IC2=IXNEXT
      NUSED=0
C
C  COUNT THE NUMBER OF OPEN CHANNELS AND SET UP WVEC ARRAY
      CALL WVCALC(WVEC,WMAX,ERED,EINT,NOPEN,N)
      IF (NOPEN.EQ.0) RETURN

C  IF NO PROPAGATION TO BE DONE, SET S MATRIX TO IDENTITY
      IF (NSEGS.EQ.0) THEN
        IJ=0
        DO I=1,NOPEN
        DO J=1,NOPEN
          IJ=IJ+1
          SR(IJ)=0.D0
          SI(IJ)=0.D0
          IF (I.EQ.J) SR(IJ)=1.D0
        ENDDO
        ENDDO
        RETURN
      ENDIF

C  SET UP TO USE UNIT  (ISCRU)
      IREAD  = ICODE.EQ.2 .AND. ISCRU.GT.0
      IWRITE = ICODE.EQ.1 .AND. ISCRU.GT.0
C     ---------------------------------------------------------------

      ISTART=0
      DRTEMP=DR
      RMNTMP=RMIN
      RMXTMP=RMAX

      PRPLEN=REND(1)-RST(1)
      IF (STEPS.GT.0.D0) THEN
        NSTEP=WMAX*STEPS*PRPLEN/ACOS(-1.D0)
        IF (IWRITE) NSTEP=WMAXIN*STEPS*PRPLEN/ACOS(-1.D0)
        NSTEP=MAX(1,NSTEP) ! THOUGH IF IT IS ZERO, THERE'S PROBABLY SOMETHING WRONG
      ELSE
        NSTEP=NINT(PRPLEN/DRPR(1))
      ENDIF

      DO ISEG=1,NSEGS
        IF (ISEG.GT.1) RST(ISEG)=RMAX
        IPROP=IFLG(ISEG)
        RMIN=RST(ISEG)
        RMAX=REND(ISEG)

        IF (ISEG.EQ.1) THEN
          IF (STEPS.GT.0.D0 .AND. IPROP.NE.1 .AND. IPROP.NE.14) THEN
            DR=PRPLEN/DBLE(NSTEP)
          ELSEIF (DRPR(ISEG).NE.-1.D0) THEN
            DR=DRPR(ISEG)
C         ELSE ! THIS BRANCH PREVENTED IN PRCHCK
          ENDIF
        ELSEIF (ISEG.GT.1 .AND. DRPR(ISEG).NE.-1.D0) THEN
          DR=DRPR(ISEG)
        ENDIF
        IF (IPROP.GT.-1) THEN
          IF (.NOT.IREAD) THEN
            IF (IWRITE) THEN
              WRITE(ISCRU) RMIN,RMAX,DR,ERED,NSTEP
            ENDIF
          ENDIF

          IF (IREAD) THEN
            READ(ISCRU) RMIN,RMAX,DR,EFIRST,NSTEP
          ENDIF
        ENDIF
C
C  INITIALISE Y MATRIX (HERE CALLED SR)
C
        IF (IPROP.EQ.5 .OR. IPROP.EQ.6 .OR. IPROP.EQ.7 .OR.
     1      IPROP.EQ.9 .OR. IPROP.EQ.1 .OR. IPROP.EQ.14) THEN
          IY2=IC2     ! EVALS
          IXNEXT=IY2+N
          CALL CHKSTR(NUSED)

          CALL YINIT(SR,U,VL,IV,P,CENT,EINT,X(IY2),SI,
     &               IPROP,N,MXLAM,NPOTL,ISCRU,ISTART,
     &               ERED,RMIN,RMLMDA,.TRUE.,IREAD,IWRITE,
     &               IPRINT)
          IXNEXT=IC2
        ENDIF
C
        IF (IPROP.EQ.1) THEN
C  AIRY ALGORITHM OF ALEXANDER AND MANOLOPOULOS
C
          IF (RMIN.LT.RMAX) THEN
            XF=RMIN
            ITWO=ICODE-1
            IF (ISCRU.EQ.0) ITWO=-1
            IT7=IC2     ! EIGOLD
            IT1=IT7+N   ! Y1
            IT2=IT1+N   ! Y2
            IT3=IT2+N   ! Y3
            IT4=IT3+N   ! Y4
            IT5=IT4+N   ! VECNOW
            IT6=IT5+NSQ ! VECNEW
            IT8=IT6+NSQ ! EIGNOW
            IT9=IT8+N   ! HP
            IXNEXT=IT9+N
            CALL CHKSTR(NUSED)
            CALL AIRPRP(N,MXLAM,NPOTL,
     1                  SR,SI,U,VL,IV,EINT,CENT,P,
     3                  X(IT1),X(IT2),X(IT3),X(IT4),X(IT5),X(IT6),
     4                  X(IT7),X(IT8),X(IT9),TOLHI,POWRX,
     5                  RMIN,XF,RMAX,DR,NSTEP,ERED,RMLMDA,
     6                  ITWO,ISCRU,IPRINT,ISTART,NODES)
            IF (XF.LT.RMAX) RMAX=XF
            IF (IPRINT.GE.3) WRITE(6,2000) 'AIRPRP',RMIN,RMAX,
     1                                     NSTEP
          ENDIF
C
        ELSEIF (IPROP.EQ.2) THEN
C  SOLVE COUPLED EQUATIONS BY METHOD OF DEVOGELAERE
C
          IT1=IC2        ! Y
          IT2=IT1+4*NSQ  ! YP
          IT3=IT2+2*NSQ  ! F
          IT4=IT3+4*NSQ  ! XM
          IT5=IT4+NSQ    ! YM
          IT6=IT5+NSQ    ! DIAG
          IXNEXT=IT6+N
          CALL CHKSTR(NUSED)
C
          CALL DVPROP(N,NSQ,MXLAM,NPOTL,
     1                SR,SI,U,VL,IV,EINT,CENT,WVEC,L,NB,P,
     2                X(IT1),X(IT2),X(IT3),X(IT4),X(IT5),X(IT6),
     3                RMIN,RMAX,DR,ICODE,IPRINT)
C
        ELSEIF (IPROP.EQ.3) THEN
C  SOLVE COUPLED EQUATIONS BY WALKER-LIGHT R-MATRIX PROPAGATOR METHOD
C
          IT1=IC2      ! W
          IT2=IT1+NSQ  ! R
          IT3=IT2+NSQ  ! EIGOLD
          IT4=IT3+N    ! EIGNOW
          IT5=IT4+N    ! DIAG
          IT6=IT5+N    ! R1/N1
          IT7=IT6+N    ! R2/N2
          IT8=IT7+N    ! R3
          IT9=IT8+N    ! R4
          IT10=IT9+N   ! CLOSED
          IXNEXT=IT10+N
          CALL CHKSTR(NUSED)
          RSTOP=RMAX
          CALL RMTPRP(N,NSQ,MXLAM,NPOTL,
     1                SR,SI,U,VL,IV,EINT,CENT,WVEC,JSINDX,L,NB,P,
     2                X(IT1),X(IT2),X(IT3),X(IT4),X(IT5),X(IT6),X(IT7),
     3                X(IT8),X(IT9),X(IT6),X(IT7),X(IT10),RSTOP,
     4                NOPMAX,DEEP,IK,ICODE,IPRINT,NV)
          RMAX=RSTOP
C
        ELSEIF (IPROP.EQ.5) THEN
C  SOLVE COUPLED EQUATIONS BY JOHNSONS LOG DERIVATIVE PROPAGATION METHOD
C
          IT1=IC2              ! DIAG
          IXNEXT=IT1+N
          CALL CHKSTR(NUSED)
          CALL LDPROP(N,MXLAM,NPOTL,
     1                SR,U,VL,IV,EINT,CENT,P,X(IT1),
     2                RMIN,RMAX,NSTEP,ERED,ESHIFT,RMLMDA,
     3                IREAD,IWRITE,ISCRU,NODES,IPRINT)
          IF (IPRINT.GE.3) WRITE(6,2000) 'LDPROP',RMIN,RMAX,NSTEP
C
        ELSEIF (IPROP.EQ.6) THEN
C  DIABATIC MODIFIED LOG DERIVATIVE ALGORITHM.
C
          ESHIFT=ERED-EFIRST
          IF (N.GT.1 .AND. NSTEP.GT.0 .AND. RMIN.LT.RMAX) THEN
            IT1=IC2       ! Y14
            IT2=IT1+N     ! Y23
            IT3=IT2+N     ! DIAG
            IXNEXT=IT3+N
            IT4=IT3+N     ! W
            IT5=IT4+NSQ   ! W2
            IT6=IT5+NSQ   ! W3
            IF (WAVE) IXNEXT=IT6+NSQ
            CALL CHKSTR(NUSED)
            CALL DWPROP(N,MXLAM,NPOTL,
     1                  SR,U,VL,IV,EINT,CENT,P,
     2                  X(IT1),X(IT2),X(IT3),X(IT4),X(IT5),X(IT6),
     3                  RMIN,RMAX,NSTEP,ERED,ESHIFT,RMLMDA,
     4                  IREAD,IWRITE,ISCRU,NODES,IREC,WAVE,IPRINT)
            IF (IPRINT.GE.3) WRITE(6,2000) 'DWPROP',RMIN,RMAX,
     1                                     NSTEP
          ELSEIF (N.EQ.1) THEN
            NPT=NSTEP+1
            IT1=IC2           ! U
            IT2=IT1+NPT       ! W
            IT3=IT2+NPT       ! Q
            IT4=IT3+NPT       ! Y1
            IT5=IT4+NPT       ! Y2
            IF1=IT5+NPT
            ITP=IT3           ! P
            IF2=ITP+NPT*MXLAM
            IXNEXT=MAX(IF1,IF2)
            NUSED=0
            CALL CHKSTR(NUSED)
            CALL ODPROP(MXLAM,NPOTL,
     1                  SR(1),VL,IV,EINT,CENT,X(ITP),
     3                  X(IT1),X(IT2),X(IT3),X(IT4),X(IT5),
     4                  RMIN,RMAX,NPT,ERED,RMLMDA,
     5                  IREAD,IWRITE,ISCRU,NODES,IPRINT)
            IF (IPRINT.GE.3) WRITE(6,2000) 'ODPROP',RMIN,RMAX,NPT-1
          ENDIF
C
        ELSEIF (IPROP.EQ.7) THEN
C  QUASIADIABATIC MODIFIED LOG DERIVATIVE ALGORITHM.
C
          IT4=IY2      ! EIVAL
          IT9=IT4+N    ! DIAG
          IT2=IT9+N    ! Q
          IT3=IT2+NSQ  ! W
          IT5=IT3+NSQ  ! Y1
          IT6=IT5+N    ! Y2
          IT7=IT6+N    ! Y3
          IT8=IT7+N    ! Y4
          IXNEXT=IT8+N
          CALL CHKSTR(NUSED)
C
          CALL QAPROP(N,NSQ,MXLAM,NPOTL,
     1                SR,SI,U,VL,IV,EINT,CENT,P,
     2                X(IT2),X(IT3),X(IT4),X(IT5),X(IT6),X(IT7),X(IT8),
     3                X(IT9),RMIN,RMAX,NSTEP,ERED,RMLMDA,
     4                IREAD,IWRITE,ISCRU,IPRINT,NODES)
          IF (IPRINT.GE.3) WRITE(6,2000) 'QAPROP',RMIN,RMAX,NSTEP
 2000     FORMAT(/2X,A,'. LOG DERIVATIVE MATRIX PROPAGATED FROM ',
     &           F12.4,'  TO',1PG12.5,'  IN ',I6,'  STEPS.')

C
        ELSEIF (IPROP.EQ.9) THEN
C  SYMPLECTIC LOG-DERIVATIVE PROPAGATOR
C
          ESHIFT=ERED-EFIRST
          IT1=IC2      ! DIAG
          IXNEXT=IT1+N
          CALL CHKSTR(NUSED)
          CALL MGPROP(N,MXLAM,NPOTL,
     1                SR,U,VL,IV,EINT,CENT,P,X(IT1),
     3                RMIN,RMAX,NSTEP,ERED,ESHIFT,RMLMDA,
     4                IREAD,IWRITE,ISCRU,NODES,IPRINT)
          IF (IPRINT.GE.3) WRITE(6,2000) 'MGPROP',RMIN,RMAX,NSTEP
C
        ELSEIF (IPROP.EQ.14) THEN
C  SOLVE COUPLED EQUATIONS BY VIVAS METHOD
C
          IF (ISTART.EQ.1) THEN
C  GET R-MATRIX BY INVERTING LOG-DERIVATIVE MATRIX OR BY DIRECT INITIALISATION
            CALL SYMINV(SR,N,N,IFAIL)
            CALL DSYFIL('U',N,SR,N)
          ELSE
            DO I=1,NSQ
              SR(I)=0.0D0
            ENDDO
            DO I=1,NSQ,N+1
              SR(I)=1D30
            ENDDO
          ENDIF
          TLDIAG=0.064D0*SQRT(TOLHI/0.001D0)
          TOFF=TLDIAG
          IT1=IC2              ! A1
          IT2=IT1+N            ! A1P
          IT3=IT2+N            ! B1
          IT4=IT3+N            ! B1P
          IT5=IT4+N            ! WKS
          IT6=IT5+N            ! G1
          IT7=IT6+N            ! G1P
          IT8=IT7+N            ! G2
          IT9=IT8+N            ! G2P
          IT10=IT9+N           ! COSX
          IT11=IT10+N          ! SINX
          IT12=IT11+N          ! SINE
          IT13=IT12+N          ! DIAG
          IT14=IT13+N          ! XK
          IT15=IT14+N          ! XSQ
          IT16=IT15+N          ! TSTORE
          IT17=IT16+NSQ        ! W0
          IT18=IT17+NSQ        ! W1
          IT19=IT18+NSQ        ! W2
          IT20=IT19+NSQ        ! EYE11
          IT21=IT20+NSQ        ! EYE12
          IT22=IT21+NSQ        ! EYE22
          IT23=IT22+NSQ        ! VECOLD
          IXNEXT=IT23+NSQ
          CALL CHKSTR(NUSED)
C
          DRNOW=DR
          CALL VIVAS(N,NSQ,MXLAM,NPOTL,
     1               SR,SI,U,VL,IV,EINT,CENT,P,
     2               X(IT1),X(IT2),X(IT3),X(IT4),X(IT5),
     3               X(IT6),X(IT7),X(IT8),X(IT9),X(IT10),X(IT11),
     4               X(IT12),X(IT13),X(IT14),X(IT15),X(IT16),X(IT17),
     5               X(IT18),X(IT19),X(IT20),X(IT21),X(IT22),X(IT23),
     6               RMIN,RMAX,DRNOW,DRMAX,TLDIAG,TOFF,ERED,ESHIFT,
     7               RMLMDA,IPRINT,IREAD,IWRITE,ISCRU)
          CALL SYMINV(SR,N,N,IFAIL)
C
        ELSEIF (IPROP.EQ.-1) THEN
C
C  SOLVE EQUATIONS BY WKB USING GAUSS-MEHLER INTEGRATION.
C  ONLY GOOD FOR ONE-CHANNEL CASES
C
          IF (N.EQ.1) GOTO 810
          WRITE(6,601) N
  601     FORMAT(/' ***** ERROR.  WKB ONLY IMPLEMENTED FOR',
     1           ' ONE-CHANNEL CASE.  TERMINATED WITH N =',I4)
          STOP
  810     IT1=IC2   ! W
          IT2=IT1+1 ! DIAG
          IXNEXT=IT2+1
          IF (NUMDER) IXNEXT=IXNEXT+2*IPDIM
          CALL CHKSTR(NUSED)
          CALL WKB(N,MXLAM,NPOTL,SR,SI,VL,IV,EINT,CENT,
     1             WVEC,L,P,X(IT1),X(IT2),IPRINT)
C
        ELSE
          WRITE(6,699) IPROP
  699     FORMAT(/' SCPROP CALLED WITH AN ILLEGAL IPROP=',I4)
          STOP
        ENDIF
        RSTOP=RMAX
        ISTART=1
        NSTEPS(ISEG)=NSTEP
      ENDDO
      DR=DRTEMP
      RMIN=RMNTMP
      RMAX=RMXTMP
C
C  -----------------------------------------------------------------------
C  END OF PROPAGATION
C
      IF (IFLG(NSEGS).GE.4 .OR. IFLG(NSEGS).EQ.1) THEN
        CALL YTRANS(SR,SI,EINT,WVEC,
     1              JSINDX,L,N,P,VL,IV,
     2              MXLAM,NPOTL,ERED,RMLMDA,DEGTOL,NOPEN,
     3              IBOUND,CENT,IPRINT,.TRUE.)
        IF (NOPEN.EQ.0) RETURN
C  IF LOG-DERIVATIVE WRITE OUT REQUESTED, HERE IS WHERE THE
C  PROPAGATION VECTORS (LDRWPV) AND THE MATRIX DATA (LDRWMD)
C  ARE WRITTEN
        IF (ILDSVU.GT.0) THEN
          IDUM=LDRWPV(ILDSVU,.TRUE.,N,JSINDX,L,EINT)
          IDUM=LDRWMD(ILDSVU,.TRUE.,N,RMAX,SR)
        ENDIF
      ENDIF
      IF (IFLG(NSEGS).GE.3 .OR. IFLG(NSEGS).EQ.1) THEN
        IT1=IC2      ! SJ
        IT2=IT1+N    ! SJP
        IT3=IT2+N    ! SN
        IT4=IT3+N    ! SNP
        IT5=IT4+N    ! WORKSPACE TO SAVE SI IN FOR WAVEFUNCTIONS
        IXNEXT=IT5
        IF (WAVE) IXNEXT=IT5+NSQ
        CALL CHKSTR(NUSED)

C  LOG-DERIVATIVE MATRIX STORED IN SR ON EXIT FROM PROPAGATOR ROUTINES
        IF (WAVE) X(IT5:IT5+NSQ-1)=SI(1:NSQ)
C
C  CONVERT LOG-DERIVATIVE MATRIX TO K MATRIX AND THEN TO S MATRIX
        CALL YTOK(NB,WVEC,L,N,NOPEN,X(IT1),X(IT2),X(IT3),X(IT4),
     2            SR,SI,U,RSTOP,CENT)
        CALL KTOS(U,SR,SI,NOPEN)
      ENDIF

      IF (WAVE) THEN
C  CALCULATE SCATTERING WAVEFUNCTION
        IT6=IT5+NSQ
        IXNEXT=IT6+NSQ
        CALL CHKSTR(NUSED)
        NTSTPS=SUM(NSTEPS(1:NSEGS))
        CALL WVSTPS(NTSTPS+1)
        CALL SCWAVE(RMIN,RMAX,WVEC,X(IT1),X(IT2),X(IT3),X(IT4),
     1              X(IT6),SR,SI,X(IT5),U,L,N,NSQ,NOPEN,NB,NTSTPS,
     2              ICHAN,IREC,IPRINT)
      ENDIF
C
C  WE ARE FINISHED WITH THIS TEMPORARY STORAGE; RESTORE IXNEXT.
C  THIS IS CONSISTENT W/ V11 WHICH DID NOT MODIFY SCPROP IC2 ARGUMENT
C  HOWEVER, THIS MEANS THAT ONE CANNOT EXPECT ALLOCATED STORAGE
C  TO BE RETAINED BEYOND A SCATTERING CALL
      IXNEXT=IC2
      IF (WAVE) THEN
        CLOSE(IPSISC)
        CLOSE(IWAVSC)
      ENDIF
      RETURN
      END
