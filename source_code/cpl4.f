      SUBROUTINE CPL4(N,MXLAM,LAM,NSTATE,JSTATE,JSINDX,L,JTOT,ATAU,
     1                VL,IPRINT,LFIRST)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS SUBROUTINE CALCULATES COUPLING MATRIX ELEMENTS FOR ITYPE=4 (CPL4)
C  & ITYPE=24 (CPL24)
C
C  CRLS 10-05-22: SAVES COUPLING COEFFICIENTS USING ALLOCATABLE ARRAYS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE IFIRST,NL12,IXMX,ISTART,CPL
C  SPECIFICATIONS FOR PARAMETER LIST
      INTEGER JSINDX(N),L(N),LAM(*),JSTATE(*)
      INTEGER IPRINT
      DIMENSION ATAU(*),VL(*)
      LOGICAL LFIRST
      ALLOCATABLE CPL(:),CPLTMP(:)
C
      INTEGER P1,Q1,P2,P
      LOGICAL LODD,L20
C
      COMMON /VLSAVE/ IVLU
C
      PARAMETER (PIFCT=4.D0*ACOS(-1.D0), EPS=1.D-8, Z0=0.D0)
C
C  STATEMENT FUNCTIONS ...
      F(NN) = DBLE(NN+NN+1)
      LODD(I) = I-2*(I/2).NE.0
C
      L20=.FALSE.
      IADD=0
      IF (LFIRST) THEN
        IFIRST=-1
        LFIRST=.FALSE.
        IF (ALLOCATED(CPL)) DEALLOCATE (CPL)
      ENDIF
      IF (IFIRST.GT.-1) GOTO 1000

      NL12=NSTATE*(NSTATE+1)/2
      IXMX=NL12*MXLAM
      ALLOCATE (CPL(IXMX))

      IX=0
!  Start of long DO loop #1
      DO LL=1,MXLAM
        P1 = LAM(4*LL-3)
        Q1 = LAM(4*LL-2)
        P2 = LAM(4*LL-1)
        P  = LAM(4*LL)
        XP1 = P1
        XQ1 = Q1
        FACL=SQRT(F(P)*F(P2))/PIFCT
!  Start of long DO loop #2
        DO IC=1,NSTATE
          JC   = JSTATE(IC)
          J1C  = JSTATE(IC + 2*NSTATE)
          J2C  = JSTATE(IC +   NSTATE)
          XJC  = JC
          XJ1C = J1C
          XJ2C = J2C
          ISTC = JSTATE(IC + 5*NSTATE)
          NKC  = JSTATE(IC + 6*NSTATE)
          FACLC=FACL*SQRT(F(J2C)*F(JC)*F(J1C))
!  Start of long DO loop #3
        DO IR=1,IC
          IX=IX+1
          JR   = JSTATE(IR)
          J1R  = JSTATE(IR + 2*NSTATE)
          J2R  = JSTATE(IR +   NSTATE)
          XJR  = JR
          XJ1R = J1R
          XJ2R = J2R
          ISTR = JSTATE(IR + 5*NSTATE)
          NKR  = JSTATE(IR + 6*NSTATE)
          XCPL=Z0
          FACRLC=FACLC*SQRT(F(J2R)*F(JR)*F(J1R))
          DO KC=1,NKC
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
            IF (ABS(ATAU(ISTC+KC)).LE.EPS) CYCLE

            KKC=KC-JC-1
            XKC=KKC
            DO KR=1,NKR
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
              IF (ABS(ATAU(ISTR+KR)).LE.EPS) CYCLE

              KKR=KR-JR-1
              XKR=KKR
              AF=ATAU(ISTR+KR)*ATAU(ISTC+KC)
              IF (LODD(KKR)) AF=-AF
              IF (KKR-KKC.EQ.Q1) THEN
                XCPL=XCPL+AF*THRJ(XJ1R,XP1,XJ1C,-XKR,XQ1,XKC)
                IF (Q1.EQ.0) CYCLE
              ENDIF
              IF (KKC-KKR.NE.Q1) CYCLE

C  ADJUST FOR (-1)**MU IN POTENTIAL. . .
              AF=AF*PARSGN(P1+Q1+P2+P)
              XCPL=XCPL+AF*THRJ(XJ1R,XP1,XJ1C,-XKR,-XQ1,XKC)
            ENDDO
          ENDDO
C  NOW GET 'CONSTANT FACTORS'
          XFCT=FACRLC*PARSGN(JR-J1C+J2C)
     1         *THREEJ(J2R,P2,J2C)*XNINEJ(JC,JR,P,J1C,J1R,P1,J2C,J2R,P2)
          CPL(IX)=XCPL*XFCT
        ENDDO
!  End of long DO loop #3
        ENDDO
!  End of long DO loop #2
      ENDDO
!  End of long DO loop #1

      IFIRST=0
C
C  NOW GET COUPLING MATRIX ELEMENTS FROM STORED PARTS
 1000 PJT=PARSGN(JTOT)
      IF (IVLU.GT.0) REWIND IVLU
      NZERO=0
!  Start of long DO loop #4
      DO LL=1,MXLAM
        P1 = LAM(4*LL-3)
        Q1 = LAM(4*LL-2)
        P2 = LAM(4*LL-1)
        P  = LAM(4*LL)
C
        PPP = PARSGN(P)
        IX1=(LL-1)*NL12
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
C
!  Start of long DO loop #5
        DO IC=1,N
          INJ12P = JSINDX(IC)
          JC=JSTATE(INJ12P)
          LC=L(IC)
C
        DO IR=1,IC
          INJ12=JSINDX(IR)
          JR=JSTATE(INJ12)
          LR=L(IR)
C
          XFACT = PJT*PPP*THREEJ(LR,P,LC)*SIXJ(LR,JR,LC,JC,JTOT,P)
     1               *SQRT(F(LR)*F(LC))
          IF (INJ12.GE.INJ12P) THEN
            IX2=INJ12*(INJ12-1)/2+INJ12P
          ELSE
            IX2=INJ12P*(INJ12P-1)/2+INJ12
          ENDIF
          INDX=IX1+IX2
C
          IF (CPL(INDX).EQ.0.D0) THEN
            VL(IX) = 0.D0
          ELSE
            VL(IX)=XFACT*CPL(INDX)
          ENDIF
          IF (VL(IX).NE.0.D0) NNZ=NNZ+1
          IF (IVLU.EQ.0) THEN
            IX=IX+MXLAM
          ELSE
            IX=IX+1
          ENDIF
        ENDDO
        ENDDO
!  End of long DO loop #5
        IF (NNZ.EQ.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,697) P1,Q1,P2,P
        ENDIF
        IF (IVLU.GT.0) WRITE(IVLU) (VL(I),I=1,N*(N+1)/2)
      ENDDO
!  End of long DO loop #4

      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) 'JTOT',JTOT,NZERO
  620 FORMAT('  * * * NOTE.  FOR ',A,' =',I4,',  ALL COUPLING ',
     1       'COEFFICIENTS ARE 0.0 FOR',I5,' POTENTIAL EXPANSION TERMS')

      RETURN
C----------------------------------------------------------------------------
      ENTRY CPL24(N,MXLAM,LAM,NSTATE,JSTATE,JSINDX,MVAL,ATAU,
     1            VL,IPRINT,LFIRST)
      L20=.TRUE.
      IADD=20
C
C  IF LFIRST IS TRUE (FIRST CALL), DO SOME INITIALIZATION
      IF (LFIRST) THEN
        IFIRST=-1
        LFIRST=.FALSE.
        IF (ALLOCATED(CPL)) DEALLOCATE (CPL)
      ENDIF
C
      IF (IFIRST.GT.-1) GOTO 1200

C  FIRST TIME THROUGH SET UP SOME STORAGE POINTERS
      NL12=NSTATE*(NSTATE+1)/2
      IXMX=NL12*MXLAM
      ALLOCATE (CPL(IXMX))
C
 1200 MVABS=ABS(MVAL)
C  SEE IF VALUES ARE STORED FOR THIS HIGH AN MVALUE
C  IF NOT, TRY TO STORE THEM IN XCPL().
      IF (MVABS.LE.IFIRST) GOTO 2200

      MV=IFIRST+1

C  TEST FOR AVAILABLE STORAGE; NEED IXMX FOR THIS MVAL
      NCPL=SIZE(CPL)
      IF (NCPL/IXMX.LE.MVABS) THEN
        ALLOCATE (CPLTMP(NCPL))
        CPLTMP=CPL
        DEALLOCATE (CPL)
        ALLOCATE (CPL(IXMX*(MVABS+1)))
        CPL(1:NCPL)=CPLTMP
        DEALLOCATE (CPLTMP)
      ENDIF

!  Start of long DO loop #6
      DO MV=IFIRST+1,MVABS
C
        IX=MV*IXMX
!  Start of long DO loop #7
        DO LL=1,MXLAM
          P1 = LAM(4*LL-3)
          Q1 = LAM(4*LL-2)
          P2 = LAM(4*LL-1)
          P =  LAM(4*LL)
!  Start of long DO loop #8
          DO IC=1,NSTATE
            JC  = JSTATE(IC)
            J1C = JSTATE(IC+2*NSTATE)
            J2C = JSTATE(IC+  NSTATE)
            ISTC= JSTATE(IC+5*NSTATE)
            NKC = JSTATE(IC+6*NSTATE)
!  Start of long DO loop #9
          DO IR=1,IC
            JR  = JSTATE(IR)
            J1R = JSTATE(IR+2*NSTATE)
            J2R = JSTATE(IR+  NSTATE)
            ISTR= JSTATE(IR+5*NSTATE)
            NKR = JSTATE(IR+6*NSTATE)
            IX=IX+1
            XCPL=Z0
            DO KC=1,NKC
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
              IF (ABS(ATAU(ISTC+KC)).LE.EPS) CYCLE

              KKC=KC-JC-1
              DO KR=1,NKR
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
                IF (ABS(ATAU(ISTR+KR)).LE.EPS) CYCLE

                KKR=KR-JR-1
                AF=ATAU(ISTR+KR)*ATAU(ISTC+KC)
                IF (KKR-KKC.EQ.Q1) THEN

                  XCPL=XCPL+AF*RSYMTP(J1R,KKR,J2R,J1C,KKC,J2C,
     1                                JR,JC,MVAL,P1,Q1,P2,P)
                  IF (Q1.EQ.0) CYCLE
                ENDIF
                IF (KKC-KKR.NE.Q1) CYCLE

C  ADJUST FOR (-1)**MU IN POTENTIAL. . .
                IF (LODD(P1+Q1+P2+P)) AF = -AF
                XCPL=XCPL+AF*RSYMTP(J1R,KKR,J2R,J1C,KKC,J2C,
     1                              JR,JC,MVAL,P1,-Q1,P2,P)
              ENDDO
            ENDDO
            CPL(IX)=XCPL
          ENDDO
!  End of long DO loop #9
          ENDDO
!  End of long DO loop #8
        ENDDO
!  End of long DO loop #7
      ENDDO
!  End of long DO loop #6

      IFIRST = MVABS
C
C  MV.GT.IFIRST ==> VALUES NOT STORED.  CALCULATE THEM VIA OLD CODE
C
 2200 IF (MVABS.GT.IFIRST) GOTO 9000

C  MVABS.LE.IFIRST.  COEFFS STORED.  FILL VL() FROM XCPL
      IXM=MVABS*IXMX
      IF (IVLU.GT.0) REWIND IVLU
      NZERO=0
!  Start of long DO loop #10
      DO LL=1,MXLAM
        P1 = LAM(4*LL-3)
        Q1 = LAM(4*LL-2)
        P2 = LAM(4*LL-1)
        P  = LAM(4*LL)
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
        DO ICOL=1,N
          I1=JSINDX(ICOL)
          JC=JSTATE(I1)
        DO IROW=1,ICOL
          I2=JSINDX(IROW)
          JR=JSTATE(I2)
          IF (I1.GT.I2) THEN
            IX12=I1*(I1-1)/2+I2
          ELSE
            IX12=I2*(I2-1)/2+I1
          ENDIF
          IXX=IXM+(LL-1)*NL12+IX12
          VL(IX)=CPL(IXX)
C  WE HAVE STORED COUPLING FOR POSITIVE MVALUES; CORRECT IF NECESSARY
C  FOR PARITY OF THRJ(JR, P ,JC, MVAL, 0, -MVAL)
          IF (MVAL.LT.0 .AND. LODD(JC+JR+P)) VL(IX)=-VL(IX)
          IF (VL(IX).NE.Z0) NNZ=NNZ+1
          IF (IVLU.EQ.0) THEN
            IX=IX+MXLAM
          ELSE
            IX=IX+1
          ENDIF
        ENDDO
        ENDDO
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) 'MVAL',MVAL,LL
        ENDIF
  612   FORMAT('  * * * NOTE.  FOR ',A,' =',I4,',  ALL COUPLING '
     1         'COEFFICIENTS ARE 0.0 FOR EXPANSION TERM',I4)
        IF (IVLU.GT.0) WRITE(IVLU) (VL(I),I=1,N*(N+1)/2)
      ENDDO
!  End of long DO loop #10

      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) 'MVAL',MVAL,NZERO

      RETURN
C
C  -------------------- OLD CODE REJOINS HERE ---------------------
C
 9000 IF (IVLU.GT.0) REWIND IVLU
C
C  ----- LOOP OVER RADIAL SURFACES -----
C
      NZERO=0
!  Start of long DO loop #11
      DO LL=1,MXLAM
        P1 = LAM(4*LL-3)
        Q1 = LAM(4*LL-2)
        P2 = LAM(4*LL-1)
        P = LAM(4*LL)
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
C
!  Start of long DO loop #12
        DO IC=1,N
          JC   = JSTATE(JSINDX(IC))
          J1C  = JSTATE(JSINDX(IC) + 2*NSTATE)
          J2C  = JSTATE(JSINDX(IC) +   NSTATE)
          ISTC = JSTATE(JSINDX(IC) + 5*NSTATE)
          NKC  = JSTATE(JSINDX(IC) + 6*NSTATE)
C
!  Start of long DO loop #13
        DO IR=1,IC
          JR   = JSTATE(JSINDX(IR))
          J1R  = JSTATE(JSINDX(IR) + 2*NSTATE)
          J2R  = JSTATE(JSINDX(IR) +   NSTATE)
          ISTR = JSTATE(JSINDX(IR) + 5*NSTATE)
          NKR  = JSTATE(JSINDX(IR) + 6*NSTATE)
C
          VL(IX)=0.D0
C
C ----- LOOP OVER EXPANSION COEFFICIENTS. -----
C ----- SKIP IMMEDIATELY IF COEFFICIENT IS SMALL. -----
C
!  Start of long DO loop #14
          DO KC=1,NKC
            IF (ABS(ATAU(ISTC+KC)).LE.EPS) CYCLE
            KKC=KC-JC-1
C
            DO KR=1,NKR
              IF (ABS(ATAU(ISTR+KR)).LE.EPS) CYCLE
              KKR=KR-JR-1
              AF=ATAU(ISTR+KR)*ATAU(ISTC+KC)
              IF (KKR-KKC.EQ.Q1) THEN
                IF (IADD.EQ.0) THEN
                  VL(IX)=VL(IX)
     1                   +AF*QSYMTP(J1R,KKR,J1C,KKC,J2R,J2C,L(IR),L(IC),
     2                              JR,JC,JTOT,P1,Q1,P2,P)
                ELSEIF (IADD.EQ.20) THEN
                  VL(IX)=VL(IX)+AF*RSYMTP(J1R,KKR,J2R,J1C,KKC,J2C,
     1                                    JR,JC,MVAL,P1,Q1,P2,P)
                ENDIF
                IF (Q1.EQ.0) CYCLE
              ENDIF
              IF (KKC-KKR.NE.Q1) CYCLE
              AF = AF*PARSGN(P1+P2+P+Q1)
              IF (IADD.EQ.0) THEN
                VL(IX)=VL(IX)
     1                 +AF*QSYMTP(J1R,KKR,J1C,KKC,J2R,J2C,L(IR),L(IC),
     2                            JR,JC,JTOT,P1,-Q1,P2,P)
              ELSEIF (IADD.EQ.20) THEN
                VL(IX)=VL(IX)+AF*RSYMTP(J1R,KKR,J2R,J1C,KKC,J2C,
     1                                  JR,JC,MVAL,P1,-Q1,P2,P)
              ENDIF
            ENDDO
C
          ENDDO
!  End of long DO loop #14
C
          IF (VL(IX).NE.0.D0) NNZ=NNZ+1
          IF (IVLU.EQ.0) THEN
            IX = IX + MXLAM
          ELSE
            IX = IX + 1
          ENDIF
        ENDDO
!  End of long DO loop #13
        ENDDO
!  End of long DO loop #12
        IF (NNZ.EQ.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,697) P1,Q1,P2,P
        ENDIF
  697   FORMAT('  * * *  NOTE.  ALL COUPLING COEFFICIENTS ARE ZERO ',
     1         'FOR P1, Q1, P2, P = ', 4I4)
      ENDDO
!  End of long DO loop #11

      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14) THEN
        IF (L20) WRITE(6,620) 'MVAL',MVAL,NZERO
        IF (.NOT.L20) WRITE(6,620) 'JTOT',JTOT,NZERO
      ENDIF

      RETURN
      END
