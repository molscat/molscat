      SUBROUTINE COUPLE(N,ITYPE,MXLAM,NHAM,LAM,NSTATE,JSTATE,JSINDX,L,
     1                  JTOT,VL,IV,IEX,IPRINT,ATAU)
C  Copyright (C) 2025 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS SUBROUTINE CALCULATES COUPLING COEFFICIENTS FOR DIFFERENT ITYPES
C
C  MODIFIED FOR ITYPE=4 BY S Green 29 JUN 94
C
C  CRLS 10-05-22: SAVES COUPLING COEFFICIENTS USING ALLOCATABLE ARRAYS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE LFIRST
      INTEGER IPRINT
      DIMENSION LAM(*),JSTATE(NSTATE,3),JSINDX(N),L(N),ATAU(*)
      DIMENSION VL(*),IV(*)
      LOGICAL LFIRST,LODD
      ALLOCATABLE CPLJ(:),CPLL(:),CPL6J(:)
C
      COMMON /VLFLAG/ IVLFL
C
      PARAMETER (EPS=1.D-10,PIFCT=(4.D0*ACOS(-1.D0))**(-1.5D0),
     1           SQRTHF=SQRT(0.5D0))
C  PIFCT IS FACTOR (4.*PI)**(-3/2)
C
C  STATEMENT FUNCTIONS
      Z(I)=DBLE(I+I+1)
      LODD(I)=I-2*(I/2).NE.0
C
      IF (ITYPE.EQ.1 .OR. ITYPE.EQ.31) GOTO 8001
      IF (ITYPE.EQ.2 .OR. ITYPE.EQ.32) GOTO 8002
      IF (ITYPE.EQ.7 .OR. ITYPE.EQ.37) GOTO 8007
      IF (ITYPE.EQ.3)  GOTO 8003
      IF (ITYPE.EQ.4)  GOTO 8004
      IF (ITYPE.EQ.5)  GOTO 8005
      IF (ITYPE.EQ.6)  GOTO 8006
      IF (ITYPE.EQ.11) GOTO 6001
      IF (ITYPE.EQ.12) GOTO 6002
      IF (ITYPE.EQ.13) GOTO 6003
      IF (ITYPE.EQ.15) GOTO 6005
      IF (ITYPE.EQ.16) GOTO 6006
      IF (ITYPE.EQ.17) GOTO 6007
      WRITE(6,698) ITYPE
  698 FORMAT(/' * * * ERROR.  COUPLING MATRIX ELEMENTS NOT IMPLEMENTED',
     1       ' FOR ITYPE =',I12)
      STOP
C
C  COUPLING FOR ATOM + RIGID LINEAR ROTOR
C  THIS VERSION BY JM Hutson, JUNE 93, TO REDUCE NON-VECTORIZABLE CODE
C
 8001 IF (IVLFL.NE.0) GOTO 9999

      NZERO=0
!  Start of long DO loop #1
      DO LL=1,MXLAM
        NNZ=0
        I=LL
        JSAV=-1
        LSAV=-1
        LLL=LAM(LL)
        ALLOCATE (CPLJ(2*LLL+1),CPLL(2*LLL+1),CPL6J(2*LLL+1))
        J6JMAX=2*LLL+1
!  Start of long DO loop #2
        DO ICOL=1,N
          JCOL=JSTATE(JSINDX(ICOL),1)
C
C  GET ALL ZERO-PROJECTION 3J SYMBOLS FOR THIS JCOL, LAMBDA
C
          IF (JCOL.NE.JSAV) THEN
            CALL J3J000(DBLE(JCOL),DBLE(LLL),IVALJ,CPLJ,XJMIN)
            JMIN=ABS(JCOL-LLL)
            JMAX=JCOL+LLL
            JSAV=JCOL
          ENDIF
C
C  GET ALL ZERO-PROJECTION 3J SYMBOLS FOR THIS LCOL, LAMBDA
C
          IF (L(ICOL).NE.LSAV) THEN
            CALL J3J000(DBLE(L(ICOL)),DBLE(LLL),IVALL,CPLL,XLMIN)
            LMIN=ABS(L(ICOL)-LLL)
            LMAX=L(ICOL)+LLL
            LSAV=L(ICOL)
          ENDIF
          LSAV6=-1
C
!  Start of long DO loop #3
        DO IROW=1,ICOL
          JROW=JSTATE(JSINDX(IROW),1)
!  Start of long IF block #1
          IF (JROW.LT.JMIN .OR. JROW.GT.JMAX .OR.
     1        L(IROW).LT.LMIN .OR. L(IROW).GT.LMAX .OR.
     2        LODD(JROW+JMAX) .OR. LODD(L(IROW)+LMAX)) THEN
            VL(I)=0.D0
          ELSE
C
C  GET ALL 6J SYMBOLS FOR THIS JCOL, LCOL, JROW, LAMBDA, JTOT,
C  CHECKING WHETHER LROW HAS CHANGED SINCE THE LAST CALL TO J6J
C
            IF (L(IROW).NE.LSAV6) THEN
              IVAL6=J6JMAX
              CALL J6J(DBLE(L(IROW)),DBLE(JTOT),DBLE(L(ICOL)),
     1                 DBLE(JCOL),DBLE(LLL),IVAL6,XJMIN6,CPL6J)
              JMIN6=INT(XJMIN6)
              LSAV6=L(IROW)
            ENDIF
            IF (JROW.LT.JMIN6 .OR. JROW.GE.JMIN6+IVAL6) THEN
              VL(I)=0.D0
            ELSE
C
C  CALCULATE THE PERCIVAL-SEATON COEFFICIENT USING THE STORED
C  3-J AND 6-J SYMBOLS.
C
C  ARRIVE HERE ONLY IF THE TRIANGLE RELATIONSHIPS ARE SATISFIED,
C  AND IF JCOL+LAMBDA+JROW AND LCOL+LAMBDA+LROW ARE EVEN.
C  NOTE THAT ONLY 3-J SYMBOLS FOR WHICH THIS IS TRUE ARE STORED.
C
              INDJ=1+(JROW-JMIN)/2
              INDL=1+(L(IROW)-LMIN)/2
              IND6=1+JROW-JMIN6
              VL(I)=SQRT(Z(JCOL)*Z(JROW)*Z(L(ICOL))*Z(L(IROW)))
     2              *CPLJ(INDJ)*CPLL(INDL)*CPL6J(IND6)
              IF (LODD(JCOL+JROW+JTOT)) VL(I)=-VL(I)
              IF (VL(I).NE.0.D0) NNZ=NNZ+1
            ENDIF
          ENDIF
!  End of long IF block #1
          I=I+NHAM
        ENDDO
!  End of long DO loop #3
        ENDDO
!  End of long DO loop #2
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) JTOT,LL
        ENDIF
  612   FORMAT('  * * * NOTE.  FOR JTOT =',I4,',  ALL COUPLING ',
     1         'COEFFICIENTS ARE 0.0 FOR SYMMETRY',I4)
        DEALLOCATE (CPLJ,CPLL,CPL6J)
      ENDDO
!  End of long DO loop #1
      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) JTOT,NZERO
      RETURN
C
C  COUPLING FOR VIBROTOR - ATOM MAKING USE OF IV() (S Green, JAN 94)
C
 8002 IF (IVLFL.LE.0) GOTO 9999

      II=NHAM*N*(N+1)/2
      DO I=1,II
        VL(I)=0.D0
        IV(I)=0
      ENDDO
C
      NZERO=0
!  Start of long DO loop #4
      DO LL=1,MXLAM
        LLL=LAM(3*LL-2)
        NV=LAM(3*LL-1)
        NV1=LAM(3*LL)
        NNZ=0
        JSAV=-1
        LSAV=-1
        J6JMAX=2*LLL+1
        ALLOCATE (CPLJ(2*LLL+1),CPLL(2*LLL+1),CPL6J(2*LLL+1))
C
        II=0
!  Start of long DO loop #5
        DO ICOL=1,N
          JVCOL=JSTATE(JSINDX(ICOL),2)
          IF (JVCOL.NE.NV .AND. JVCOL.NE.NV1) THEN
            II=II+ICOL
            CYCLE
          ENDIF
          JCOL=JSTATE(JSINDX(ICOL),1)
C
C  GET ALL ZERO-PROJECTION 3J SYMBOLS FOR THIS JCOL, LAMBDA
C
          IF (JCOL.NE.JSAV) THEN
            CALL J3J000(DBLE(JCOL),DBLE(LLL),IVALJ,CPLJ,XJMIN)
            JMIN=ABS(JCOL-LLL)
            JMAX=JCOL+LLL
            JSAV=JCOL
          ENDIF
C
C  GET ALL ZERO-PROJECTION 3J SYMBOLS FOR THIS LCOL, LAMBDA
C
          IF (L(ICOL).NE.LSAV) THEN
            CALL J3J000(DBLE(L(ICOL)),DBLE(LLL),IVALL,CPLL,XLMIN)
            LMIN=ABS(L(ICOL)-LLL)
            LMAX=L(ICOL)+LLL
            LSAV=L(ICOL)
          ENDIF
          LSAV6=-1
C
!  Start of long DO loop #6
          DO IROW=1,ICOL
            JVROW=JSTATE(JSINDX(IROW),2)
            JROW =JSTATE(JSINDX(IROW),1)
            II=II+1
            I=(II-1)*NHAM+LLL+1
!  Start of long IF block #2
            IF (.NOT. (JROW.LT.JMIN .OR. JROW.GT.JMAX .OR.
     1                 L(IROW).LT.LMIN .OR. L(IROW).GT.LMAX .OR.
     2                 LODD(JROW+JMAX) .OR. LODD(L(IROW)+LMAX))
     3          .AND. ((JVCOL.EQ.NV  .AND. JVROW.EQ.NV1) .OR.
     5                 (JVCOL.EQ.NV1 .AND. JVROW.EQ.NV))) THEN
C
C  GET ALL 6J SYMBOLS FOR THIS JCOL, LCOL, JROW, LAMBDA, JTOT,
C  CHECKING WHETHER LROW HAS CHANGED SINCE THE LAST CALL TO J6J
C
              IF (L(IROW).NE.LSAV6) THEN
                IVAL6=J6JMAX
                CALL J6J(DBLE(L(IROW)),DBLE(JTOT),DBLE(L(ICOL)),
     1                   DBLE(JCOL),DBLE(LLL),IVAL6,XJMIN6,CPL6J)
                JMIN6=INT(XJMIN6)
                LSAV6=L(IROW)
              ENDIF
C
C  CALCULATE THE PERCIVAL-SEATON COEFFICIENT USING THE STORED
C  3-J AND 6-J SYMBOLS.
C
C  ARRIVE HERE ONLY IF THE TRIANGLE RELATIONSHIPS ARE SATISFIED,
C  AND IF JCOL+LAMBDA+JROW AND LCOL+LAMBDA+LROW ARE EVEN.
C  NOTE THAT ONLY 3-J SYMBOLS FOR WHICH THIS IS TRUE ARE STORED.
C
              IF (JROW.GE.JMIN6 .AND. JROW.LT.JMIN6+IVAL6) THEN
                INDJ=1+(JROW-JMIN)/2
                INDL=1+(L(IROW)-LMIN)/2
                IND6=1+JROW-JMIN6
                VL(I)=SQRT(Z(JCOL)*Z(JROW)*Z(L(ICOL))*Z(L(IROW)))
     2               *CPLJ(INDJ)*CPLL(INDL)*CPL6J(IND6)
                IF (LODD(JCOL+JROW+JTOT)) VL(I)=-VL(I)
                IF (VL(I).NE.0.D0) THEN
                  IV(I)=LL
                  NNZ=NNZ+1
                ELSE
                  IV(I)=0
                ENDIF
              ENDIF
            ENDIF
!  End of long IF block #2
          ENDDO
!  End of long DO loop #6
        ENDDO
!  End of long DO loop #5
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) JTOT,LL
        ENDIF
        DEALLOCATE (CPLJ,CPLL,CPL6J)
      ENDDO
!  End of long DO loop #4
      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) JTOT,NZERO
      RETURN
C
C  COUPLING MATRIX ELEMENTS FOR LINEAR ROTOR - LINEAR ROTOR.
C  THESE ARE EVALUATED BY CPL3 USING STORED JTOT-INDEPENDENT
C  PARTS.  LFIRST INDICATES WHETHER THESE ARE ALREADY STORED.
C  TO ALLOW STACKING &INPUT DECKS W/LASTIN=0, LFIRST MUST BE
C  RESET BY CALL TO ENTRY COUPLX FOR EACH SET OF INPUT.
C
 8003 IF (IVLFL.NE.0) GOTO 9999

      CALL CPL3(N,MXLAM,LAM,NSTATE,JSTATE,JSINDX,L,JTOT,IEX,
     1          VL,IPRINT,LFIRST)
      RETURN
C
C *** ITYPE4 OBTAINED VIA CALL TO CPL4
 8004 CALL CPL4(N,MXLAM,LAM,NSTATE,JSTATE,JSINDX,L,JTOT,ATAU,
     1          VL,IPRINT,LFIRST)
      RETURN
C
C *** ITYPE = 5  - NEAR SYMMETRIC TOP CODE
C    N.B. JSTATE(I,) HAS J, ABS(K), PARITY.
C *** MODIFIED SEPT. 75 FOR ODD MU VALUES . . .
C
 8005 IF (IVLFL.NE.0) GOTO 9999

      NZERO=0
!  Start of long DO loop #7
      DO LL=1,MXLAM
        NNZ=0
        I=LL
        LM=LAM(2*LL-1)
        MU=LAM(2*LL)
!  Start of long DO loop #8
        DO ICOL=1,N
          J1 =JSTATE(JSINDX(ICOL),1)
          K1 =JSTATE(JSINDX(ICOL),2)
          IS1=JSTATE(JSINDX(ICOL),3)
!  Start of long DO loop #9
        DO IROW=1,ICOL
          J2 =JSTATE(JSINDX(IROW),1)
          K2 =JSTATE(JSINDX(IROW),2)
          IS2=JSTATE(JSINDX(IROW),3)
          VL(I)=0.D0
          PARFCT=(1.D0+PARSGN(J1+J2+IS1+IS2+LM+MU))*0.5D0
          IF (PARFCT.LT.EPS) THEN
            I=I+NHAM
            CYCLE
          ENDIF
C  SPECIAL NORMALIZATION FOR K1 AND/OR K2 =0.
          IF (K1.EQ.0) PARFCT=PARFCT*SQRTHF
          IF (K2.EQ.0) PARFCT=PARFCT*SQRTHF
          KDIF=K2-K1
          IF (ABS(KDIF).EQ.MU) THEN

            WPAR=1.D0
            IF (KDIF.LT.0) WPAR=PARSGN(MU)
C  CONTRIBUTION FROM (J1,K1,L1/Y(LM,MU)/J2,K2,L2).
            VL(I)=VL(I)+WPAR*PARFCT*
     &                  FSYMTP(J1,K1,L(ICOL),J2,K2,L(IROW),JTOT,LM,KDIF)
          ENDIF
          KSUM=K2+K1
          IF (ABS(KSUM).EQ.MU) 

C  CONTRIBUTION FROM (J1,-K1,L1/ Y(LM,MU) / J2,K2,L2)
C  N.B. FOR K1=0 AND/OR K2=0, WE RECOMPUTE SAME FSYMTP.
     1      VL(I)=VL(I)+PARFCT*PARSGN(IS1)*
     2                 FSYMTP(J1,-K1,L(ICOL),J2,K2,L(IROW),JTOT,LM,KSUM)
          IF (VL(I).NE.0.D0) NNZ=NNZ+1
          I=I+NHAM
        ENDDO
!  End of long DO loop #9
        ENDDO
!  End of long DO loop #8
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) JTOT,LL
        ENDIF
      ENDDO
!  End of long DO loop #7
      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) JTOT,NZERO
      RETURN
C
C *** ITYPE=6 OBTAINED VIA CALL TO SET6/CPL6
C
 8006 CALL CPL6(N,MXLAM,LAM,NSTATE,JSINDX,L,JTOT,
     1          VL,IPRINT,LFIRST)
      RETURN
C
C *** ITYPE=7 MAKES NON-TRIVIAL USE OF THE IV ARRAY
C
 8007 IF (IVLFL.LE.0) GOTO 9999

      II=NHAM*N*(N+1)/2
      DO I=1,II
        VL(I)=0.D0
        IV(I)=0
      ENDDO
C
      NZERO=0
!  Start of long DO loop #10
      DO LL=1,MXLAM
        LLL=LAM(5*LL-4)
        NV=LAM(5*LL-3)
        NJ=LAM(5*LL-2)
        NV1=LAM(5*LL-1)
        NJ1=LAM(5*LL)
        NNZ=0
        JSAV=-1
        LSAV=-1
        ALLOCATE (CPLJ(2*LLL+1),CPLL(2*LLL+1),CPL6J(2*LLL+1))
        J6JMAX=2*LLL+1
C
        II=0
!  Start of long DO loop #11
        DO ICOL=1,N
          JVCOL=JSTATE(JSINDX(ICOL),2)
          IF (JVCOL.NE.NV .AND. JVCOL.NE.NV1) THEN
            II=II+ICOL
            CYCLE
          ENDIF
          JCOL=JSTATE(JSINDX(ICOL),1)
C
C  GET ALL ZERO-PROJECTION 3J SYMBOLS FOR THIS JCOL, LAMBDA
C
          IF (JCOL.NE.JSAV) THEN
            CALL J3J000(DBLE(JCOL),DBLE(LLL),IVALJ,CPLJ,XJMIN)
            JMIN=ABS(JCOL-LLL)
            JMAX=JCOL+LLL
            JSAV=JCOL
          ENDIF
C
C  GET ALL ZERO-PROJECTION 3J SYMBOLS FOR THIS LCOL, LAMBDA
C
          IF (L(ICOL).NE.LSAV) THEN
            CALL J3J000(DBLE(L(ICOL)),DBLE(LLL),IVALL,CPLL,XLMIN)
            LMIN=ABS(L(ICOL)-LLL)
            LMAX=L(ICOL)+LLL
            LSAV=L(ICOL)
          ENDIF
          LSAV6=-1
C
!  Start of long DO loop #12
          DO IROW=1,ICOL
            JVROW=JSTATE(JSINDX(IROW),2)
            JROW =JSTATE(JSINDX(IROW),1)
            II=II+1
            I=(II-1)*NHAM+LLL+1
!  Start of long IF block #3
            IF (.NOT. (JROW.LT.JMIN .OR. JROW.GT.JMAX .OR.
     1                 L(IROW).LT.LMIN .OR. L(IROW).GT.LMAX .OR.
     2                 LODD(JROW+JMAX) .OR. LODD(L(IROW)+LMAX))
     3          .AND. ((JVCOL.EQ.NV  .AND. JCOL.EQ.NJ .AND.
     5                  JVROW.EQ.NV1 .AND. JROW.EQ.NJ1) .OR.
     6                 (JVCOL.EQ.NV1 .AND. JCOL.EQ.NJ1 .AND.
     7                  JVROW.EQ.NV  .AND. JROW.EQ.NJ))) THEN
C
C  GET ALL 6J SYMBOLS FOR THIS JCOL, LCOL, JROW, LAMBDA, JTOT,
C  CHECKING WHETHER LROW HAS CHANGED SINCE THE LAST CALL TO J6J
C
              IF (L(IROW).NE.LSAV6) THEN
                IVAL6=J6JMAX
                CALL J6J(DBLE(L(IROW)),DBLE(JTOT),DBLE(L(ICOL)),
     1                   DBLE(JCOL),DBLE(LLL),IVAL6,XJMIN6,CPL6J)
                JMIN6=INT(XJMIN6)
                LSAV6=L(IROW)
              ENDIF
C
C  CALCULATE THE PERCIVAL-SEATON COEFFICIENT USING THE STORED
C  3-J AND 6-J SYMBOLS.
C
C  ARRIVE HERE ONLY IF THE TRIANGLE RELATIONSHIPS ARE SATISFIED,
C  AND IF JCOL+LAMBDA+JROW AND LCOL+LAMBDA+LROW ARE EVEN.
C  NOTE THAT ONLY 3-J SYMBOLS FOR WHICH THIS IS TRUE ARE STORED.
C
              IF (JROW.GE.JMIN6 .AND. JROW.LT.JMIN6+IVAL6) THEN
                INDJ=1+(JROW-JMIN)/2
                INDL=1+(L(IROW)-LMIN)/2
                IND6=1+JROW-JMIN6
                VL(I)=SQRT(Z(JCOL)*Z(JROW)*Z(L(ICOL))*Z(L(IROW)))
     2                *CPLJ(INDJ)*CPLL(INDL)*CPL6J(IND6)
                IF (LODD(JCOL+JROW+JTOT)) VL(I)=-VL(I)
                IF (VL(I).NE.0.D0) THEN
                  IV(I)=LL
                  NNZ=NNZ+1
                ELSE
                  IV(I)=0
                ENDIF
              ENDIF
            ENDIF
!  End of long IF block #3
          ENDDO
!  End of long DO loop #12
        ENDDO
!  End of long DO loop #11
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) JTOT,LL
        ENDIF
        DEALLOCATE (CPLJ,CPLL,CPL6J)
      ENDDO
!  End of long DO loop #10
      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) JTOT,NZERO
  620 FORMAT('  * * * NOTE.  FOR JTOT =',I4,',  ALL COUPLING ',
     1       'COEFFICIENTS ARE 0.0 FOR',I5,' POTENTIAL EXPANSION TERMS')
      RETURN
C
C  CODING BELOW IS FOR EFFECTIVE POTENTIAL METHOD OF H. RABITZ.
C  N.B. MATRIX ELEMENTS ARE INDEPENDENT OF JTOT (PARTIAL WAVE)
C  AND COULD BE COMPUTED ONCE AND SAVED.
C
 6001 IF (IVLFL.NE.0) GOTO 9999

      DO LL=1,MXLAM
        NNZ=0
        I=LL
        DO ICOL=1,N
          J1P=JSTATE(JSINDX(ICOL),1)
        DO IROW=1,ICOL
          J1 =JSTATE(JSINDX(IROW),1)
          VL(I)=PARSGN((ABS(J1P-J1)+J1P+J1)/2) *
     1          SQRT(SQRT(Z(J1P)*Z(J1))/Z(LAM(LL))) *
     2          THREEJ(J1P,LAM(LL),J1)
          IF (VL(I).NE.0.D0) NNZ=NNZ+1
          I=I+NHAM
        ENDDO
        ENDDO
        IF (NNZ.EQ.0) WRITE(6,612) JTOT,LL
      ENDDO
      RETURN

 6002 IF (IVLFL.LE.0) GOTO 9999

      NZERO=NHAM*N*(N+1)/2
      DO I=1,NZERO
        VL(I)=0.D0
        IV(I)=0
      ENDDO

      NZERO=0
      DO LL=1,MXLAM
        LLL=LAM(3*LL-2)
        NV =LAM(3*LL-1)
        NV1=LAM(3*LL)
        NNZ=0
        II=0
        DO ICOL=1,N
          NVC=JSTATE(JSINDX(ICOL),2)
          NJC=JSTATE(JSINDX(ICOL),1)
        DO IROW=1,ICOL
          NVR=JSTATE(JSINDX(IROW),2)
          NJR=JSTATE(JSINDX(IROW),1)
          II=II+1
          IF ((NV.EQ.NVC .AND. NV1.EQ.NVR) .OR.
     1        (NV.EQ.NVR .AND. NV1.EQ.NVC)) THEN
            I=(II-1)*NHAM+LLL+1
            VL(I)=PARSGN((ABS(NJC-NJR)+NJC+NJR)/2) *
     1            SQRT(SQRT(Z(NJC)*Z(NJR))/Z(LLL))*THREEJ(NJC,LLL,NJR)
            IV(I)=LL
            IF (VL(I).NE.0.D0) NNZ=NNZ+1
          ENDIF
        ENDDO
        ENDDO
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) JTOT,LL
        ENDIF
      ENDDO
      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) JTOT,NZERO
      RETURN

 6003 IF (IVLFL.NE.0) GOTO 9999

      NZERO=0
!  Start of long DO loop #13
      DO LL=1,MXLAM
        NNZ=0
        I=LL
        LM1=LAM(3*LL-2)
        LM2=LAM(3*LL-1)
        LM=LAM(3*LL)
        DO ICOL=1,N
          J1P=JSTATE(JSINDX(ICOL),1)
          J2P=JSTATE(JSINDX(ICOL),2)
        DO IROW=1,ICOL
          J1=JSTATE(JSINDX(IROW),1)
          J2=JSTATE(JSINDX(IROW),2)
          PARFCT=PARSGN((ABS(J1+J2-J1P-J2P)+J1+J2+J1P+J2P)/2)
     1           *PIFCT*SQRT(Z(LM)*SQRT(Z(J1)*Z(J2)*Z(J1P)*Z(J2P)))
          VL(I) = PARFCT*THREEJ(J1,LM1,J1P)*THREEJ(J2,LM2,J2P)
          IF (IEX.NE.0) THEN
C
C *** N.B. THE FORMULATION BELOW ASSUMES THAT POTENTIAL IS SYMMETRIC TO
C ***   INTERCHANGE OF L1, L2.  I.E. A(L1,L2,L) = A(L2,L1,L) MUST BOTH
C ***   BE PRESENT IN INTERACTION POTENTIAL.
            VL(I)=VL(I)+PARSGN(IEX+JTOT)*PARFCT
     1                 *THREEJ(J1,LM1,J2P)*THREEJ(J2,LM2,J1P)
            IF (J1 .EQ.J2)  VL(I)=VL(I)*SQRTHF
            IF (J1P.EQ.J2P) VL(I)=VL(I)*SQRTHF
          ENDIF
 6093     IF (VL(I).NE.0.D0) NNZ=NNZ+1
          I=I+NHAM
        ENDDO
        ENDDO
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) JTOT,LL
        ENDIF
      ENDDO
!  End of long DO loop #13
      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) JTOT,NZERO
      RETURN

 6005 IF (IVLFL.NE.0) GOTO 9999

      NZERO=0
!  Start of long DO loop #14
      DO LL=1,MXLAM
        NNZ=0
        I=LL
        LM=LAM(2*LL-1)
        MU=LAM(2*LL)
        DO ICOL=1,N
          J1 =JSTATE(JSINDX(ICOL),1)
          K1 =JSTATE(JSINDX(ICOL),2)
          IS1=JSTATE(JSINDX(ICOL),3)
        DO IROW=1,ICOL
          J2 =JSTATE(JSINDX(IROW),1)
          K2 =JSTATE(JSINDX(IROW),2)
          IS2=JSTATE(JSINDX(IROW),3)
          VL(I)=0.D0
          PARFCT=(1.D0+PARSGN(J1+J2+IS1+IS2+LM+MU))*0.5D0
          IF (PARFCT.LT.EPS) THEN
            I=I+NHAM
            CYCLE
          ENDIF
          IF (K1.EQ.0) PARFCT=PARFCT*SQRTHF
          IF (K2.EQ.0) PARFCT=PARFCT*SQRTHF
          KDIF=K2-K1
          IF (ABS(KDIF).EQ.MU) THEN

            WPAR=1.D0
            IF (KDIF.LT.0) WPAR=PARSGN(MU)
            VL(I)=VL(I)+WPAR*PARFCT*ESYMTP(J1,K1,J2,K2,LM,KDIF)
          ENDIF
          KSUM=K2+K1
          IF (ABS(KSUM).EQ.MU)
C (J1, -K1 / Y(LM,MU) / J2, K2) - - - - -
     1      VL(I)=VL(I)+PARFCT*PARSGN(IS1)*ESYMTP(J1,-K1,J2,K2,LM,KSUM)
          IF (VL(I).NE.0.D0) NNZ=NNZ+1
          I=I+NHAM
        ENDDO
        ENDDO
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) JTOT,LL
        ENDIF
      ENDDO
!  End of long DO loop #14
      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) JTOT,NZERO
      RETURN
C
 6006 CALL CPL16(N,MXLAM,LAM,NSTATE,JSINDX,L,VL,IV)
      RETURN
C
 6007 IF (IVLFL.LE.0) GOTO 9999

      NZERO=NHAM*N*(N+1)/2
      DO I=1,NZERO
        VL(I)=0.D0
        IV(I)=0
      ENDDO

      NZERO=0
!  Start of long DO loop #15
      DO LL=1,MXLAM
        LLL=LAM(5*LL-4)
        NV =LAM(5*LL-3)
        NJ =LAM(5*LL-2)
        NV1=LAM(5*LL-1)
        NJ1=LAM(5*LL)
        NNZ=0
        II=0
        DO ICOL=1,N
          NVC=JSTATE(JSINDX(ICOL),2)
          NJC=JSTATE(JSINDX(ICOL),1)
        DO IROW=1,ICOL
          NVR=JSTATE(JSINDX(IROW),2)
          NJR=JSTATE(JSINDX(IROW),1)
          II=II+1
          IF (.NOT.((NV.EQ.NVC .AND. NJ.EQ.NJC .AND.
     1               NV1.EQ.NVR .AND. NJ1.EQ.NJR) .OR.
     2              (NV.EQ.NVR .AND. NJ.EQ.NJR .AND.
     3               NV1.EQ.NVC .AND. NJ1.EQ.NJC))) CYCLE
          I=(II-1)*NHAM+LLL+1
          VL(I)=PARSGN((ABS(NJC-NJR)+NJC+NJR)/2) *
     1          SQRT(SQRT(Z(NJC)*Z(NJR))/Z(LLL))*THREEJ(NJC,LLL,NJR)
          IV(I)=LL
          NNZ=NNZ+1
        ENDDO
        ENDDO
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) JTOT,LL
        ENDIF
        NZERO=NZERO+1
      ENDDO
!  End of long DO loop #15
      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) JTOT,NZERO
      RETURN
C
 9999 WRITE(6,699) IVLFL,ITYPE
  699 FORMAT(/'  COUPLE (JAN 93).  IVLFL =',I6,
     1        '  INCONSISTENT WITH ITYPE =',I6)
      STOP
C
      ENTRY COUPLX
      LFIRST=.TRUE.
      RETURN
      END
