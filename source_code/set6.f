      SUBROUTINE SET6(LEVIN,EIN,NSTATE,EFACT,EUNAME,IASYMU,IPRINT)
C  Copyright (C) 2025 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE basis_data, ONLY: ELEVEL, EMAX, JLEVEL, JMAX, JMIN,
     b                      JSTEP, NLEVEL, MXELVL, ROTI
      USE pair_state, ONLY: ATAU, JSTATE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  CRLS 06-2022: CONVERTED TO USE ALLOCATABLE ARRAYS, INCLUDING THE
C                JSTATE ARRAY WHICH IS NOW IN pair_state MODULE
C  CRLS 08-2022: REARRANGED SO THAT CODE COMMON TO SET4 AND SET6 IS IN
C                SUBROUTINES
C
C  REVISED FOR VERSION 14:
C  POSSIBLE METHODS OF SPECIFYING ASYMMETRIC TOP LEVELS
C    1. A,B,C .GT.0 IMPLIES GENERATE VIA SET6C
C    2. A=B=C=0 IMPLIES READ FROM IASYMU
C      a) NLEVEL.GE.0 IMPLIES READ FROM IASYMU (FILTER ON JMIN,JMAX,
C         JSTEP,EMAX); IF (NLEVEL.EQ.0) READ TO END-OF-FILE
C      b) NLEVEL.LT.0 IMPLIES READ FROM IASYMU BUT ACCEPT ONLY THOSE
C         J,ITAU CORRESPONDING TO JLEVEL(2*I-1),JLEVEL(2*I),
C         I=1,ABS(NLEVEL)
C
C  BELOW REPLACES GENERIC SAVE IN V11, WHICH APPEARED UNNECESSARY.
      SAVE IFIRST,NOMEM,NL12,IXMX,ISTART,CPL,LAMMAX
C
C  THIS ROUTINE HANDLES INPUT, ALSO MATRIX ELEMENTS FOR ITYPE=6.
C  LATTER ARE OBTAINED VIA ENTRIES CPL16, CPL6, CPL26.
C  FIRST VERSION WRITTEN AT MPI, MUNCHEN, JULY 1976.
C    CPL16 (EFF. POTL) COULD BE CHANGED, BUT PROBABLY NO LONGER USED
C
C  N.B. NKVAL HERE COULD BE OBTAINED AS NEEDED - NK=2*J+1.
C       THIS CODE IS MORE FLEXIBLE AS NOT ALL NEED BE STORED BUT
C       K-VALUE COULD BE OBTAINED VIA ADDITIONAL VECTOR KVAL(IST+1).
C
      CHARACTER(8) EUNAME
      LOGICAL LEVIN,EIN
      LOGICAL LIN
      LOGICAL NOMEM,LODD
      INTEGER IV(*)
      ALLOCATABLE JTEMP(:),ATEMP(:),CPL(:),CPL3J(:),CPL6J(:),CPLTMP(:)
C
C  SPECIFICATIONS FOR CPL16, CPL6, CPL26 ENTRIES.
      INTEGER JSINDX(N),L(N),LAM(*)
      DIMENSION VL(*)
      INTEGER IPRINT
      LOGICAL LFIRST,L20
      CHARACTER(1) PLUR(2)
C
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
C  NOTE: EIN=.TRUE. IMPLIES AT LEAST ONE VALUE IN ELEVEL IS NON-ZERO
      IF (A.GT.0.D0 .AND. B.GT.0.D0 .AND. C.GT.0.D0) THEN
        CALL SET6C(NSTATE,EIN,EFACT,EUNAME,IPRINT,-IASYMU)
        RETURN

      ENDIF
C
C  OTHERWISE, INPUT FROM UNIT IASYMU
      IF (IASYMU.LE.0 .OR. IASYMU.GE.100) THEN

        WRITE(6,601) IASYMU,IDU
  601   FORMAT(/'  ILLEGAL UNIT =',I12,'  SPECIFIED FOR IASYMU, ',
     1         'DEFAULTED TO ',I4)
        IASYMU=IDU
      ENDIF
C
      ALLOCATE (ATEMP(1))
      CALL RASYMU(NSTATE,IASYMU,ATEMP,EFACT,EUNAME,IPRINT,.TRUE.)
      DEALLOCATE (ATEMP)

C  CHECK THAT ENERGIES ARE NOT ALL IDENTICALLY ZERO.
      DO I=1,NSTATE
        IF (ELEVEL(I).NE.0.D0) RETURN
      ENDDO

      IF (NLEVEL.GT.1) THEN
        IF (IPRINT.GE.1) WRITE(6,609)
  609   FORMAT('  *** WARNING.  SET6. ENERGIES ARE ALL ZERO')
      ENDIF

      RETURN
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * END OF SET6 * *
C
C  THESE ENTRY POINTS COMPUTE COUPLING MATRIX ELEMENTS . . .
C
      ENTRY CPL16(N,MXLAM,LAM,NSTATE,JSINDX,L,VL,IV)
      IADD=10
      L20=.FALSE.
      GOTO 3000
C
C ------------------------------------------------------ END OF CPL16
      ENTRY CPL6(N,MXLAM,LAM,NSTATE,JSINDX,L,JTOT,
     1           VL,IPRINT,LFIRST)
      IF (LFIRST) THEN
        IFIRST=-1
        LFIRST=.FALSE.
        NOMEM=.FALSE.
        IF (ALLOCATED(CPL)) DEALLOCATE (CPL)
        LAMMAX=0
        DO LL=1,MXLAM
          LAMMAX=MAX(LAMMAX,LAM(2*LL-1))
        ENDDO
      ENDIF
      L20=.FALSE.
      IF (IFIRST.GT.-1) GOTO 5500

      IF (NOMEM) GOTO 5900

      NL12=NSTATE*(NSTATE+1)/2
      IXMX=NL12*MXLAM

      ALLOCATE (CPL(IXMX))

 5100 IX=0
!  Start of long DO loop #1
      DO LL=1,MXLAM
        LM=LAM(2*LL-1)
        XLM=LM
        MU=LAM(2*LL)
        XMU=MU
!  Start of long DO loop #2
        DO IC=1,NSTATE
          JC  =JSTATE(IC)
          XJC=JC
          ISTC=JSTATE(IC+3*NSTATE)
          NKC =JSTATE(IC+4*NSTATE)
!  Start of long DO loop #3
        DO IR=1,IC
          IX=IX+1
          JR  =JSTATE(IR)
          XJR=JR
          ISTR=JSTATE(IR+3*NSTATE)
          NKR =JSTATE(IR+4*NSTATE)
          XCPL=Z0
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
              IF (KKR-KKC.EQ.MU) THEN
                XCPL=XCPL+AF*THRJ(XJC,XJR,XLM,XKC,-XKR,XMU)
                IF (MU.EQ.0) CYCLE
              ENDIF
              IF (KKC-KKR.NE.MU) CYCLE

C  ADJUST FOR (-1)**MU IN POTENTIAL. . .
              IF (LODD(MU)) AF=-AF
              XCPL=XCPL+AF*THRJ(XJC,XJR,XLM,XKC,-XKR,-XMU)
            ENDDO
          ENDDO
C  NOW GET 'CONSTANT FACTORS'
          XFCT=PARSGN(JC+JR)*SQRT((F(JC)*F(JR)*F(LM))/(4.D0*PI))
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
 5500 PJT=PARSGN(JTOT)
      IF (IVLU.GT.0) REWIND IVLU
C
C  STORAGE FOR 3J AND 6J SYMBOLS
C
      ALLOCATE (CPL3J(2*LAMMAX+1),CPL6J(2*LAMMAX+1))
!  Start of long DO loop #4
      DO LL=1,MXLAM
        LM=LAM(2*LL-1)
        MU=LAM(2*LL)
        IX1=(LL-1)*NL12
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
C
        LSAV=-1
!  Start of long DO loop #5
        DO IC=1,N
          LEVC=JSINDX(IC)
          JC=JSTATE(LEVC)
          LC=L(IC)
          IF (LC.NE.LSAV) THEN
            CALL J3J000(DBLE(LC),DBLE(LM),IVALL,CPL3J,XLMIN)
            LMIN=ABS(LC-LM)
            LMAX=LC+LM
            LSAV=LC
          ENDIF
C
          LSAV6=-1
!  Start of long DO loop #6
        DO IR=1,IC
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
          IF (CPL(INDX).EQ.0.D0  .OR. LR.LT.LMIN .OR. LR.GT.LMAX
     2                               .OR. LODD(LR+LMAX)) THEN
            VL(IX)=0.D0

          ELSE
            IF (LR.NE.LSAV6) THEN
              IVAL6=2*LAMMAX+1
              CALL J6J(DBLE(LR),DBLE(JTOT),DBLE(LC),DBLE(JC),DBLE(LM),
     1                 IVAL6,XJMIN6,CPL6J)
              JMIN6=NINT(XJMIN6)
              LSAV6=LR
            ENDIF
            IF (JR.LT.JMIN6 .OR. JR.GE.JMIN6+IVAL6) THEN
              VL(IX)=0.D0
            ELSE
              INDL=(LR-LMIN)/2+1
              IND6=JR-JMIN6+1
              VL(IX)=PJT*SQRT(F(LC)*F(LR))*CPL(INDX)*
     &                   CPL3J(INDL)*CPL6J(IND6)
            ENDIF
          ENDIF
          IF (VL(IX).NE.0.D0) NNZ=NNZ+1
          IF (IVLU.EQ.0) THEN
            IX=IX+MXLAM
          ELSE
            IX=IX+1
          ENDIF
        ENDDO
!  End of long DO loop #6
        ENDDO
!  End of long DO loop #5
        IF (NNZ.EQ.0) WRITE(6,697) LM,MU
        IF (IVLU.GT.0) WRITE(IVLU) (VL(I),I=1,N*(N+1)/2)
      ENDDO
!  End of long DO loop #4
      DEALLOCATE (CPL3J,CPL6J)
      RETURN
C
C  IF WE CANNOT STORE PARTIAL COUPLING MATRIX, RECALCULATE.
 5900 IADD=0
      GOTO 3000
C--------------------------------------------------------- END OF CPL6
C
      ENTRY CPL26(N,MXLAM,LAM,NSTATE,JSINDX,MVAL,
     1            VL,IPRINT,LFIRST)
C
C  IF LFIRST IS TRUE (FIRST CALL), DO SOME INITIALIZATION
      IF (LFIRST) THEN
        IFIRST=-1
        LFIRST=.FALSE.
        NOMEM=.FALSE.
        IF (ALLOCATED(CPL)) DEALLOCATE (CPL)
      ENDIF
      L20=.TRUE.
      IF (NOMEM) GOTO 4900
C
      IF (IFIRST.GT.-1) GOTO 4500

C  FIRST TIME THROUGH SET UP SOME STORAGE POINTERS
      NL12=NSTATE*(NSTATE+1)/2
      IXMX=NL12*MXLAM
      ALLOCATE (CPL(IXMX))
C
 4500 MVABS=ABS(MVAL)
C  SEE IF VALUES ARE STORED FOR THIS HIGH AN MVALUE
C  IF NOT, TRY TO STORE THEM IN CPL().
      IF (MVABS.LE.IFIRST) GOTO 4900

C  CRLS 06-2022: EVERY TIME WE ARE ABOUT TO OVERFLOW SPACE ALLOWED FOR CPL	
C                1) ALLOCATE A TEMPORARY ARRAY THE SAME SIZE AS CPL
C                2) COPY CPL TO THE TEMPORARY ARRAY,
C                3) DEALLOCATE CPL
C                4) REALLOCATE A LARGER-SIZED CPL ARRAY
C                5) COPY THE CONTENTS OF THE TEMPORARY ARRAY BACK TO CPL
C                6) DEALLOCATE THE TEMPORARY ARRAY
      NCPL=SIZE(CPL)
      IF (NCPL/IXMX.LE.MVABS) THEN
        ALLOCATE (CPLTMP(NCPL))
        CPLTMP=CPL
        DEALLOCATE (CPL)
        ALLOCATE (CPL(IXMX*(MVABS+1)))
        CPL(1:NCPL)=CPLTMP
        DEALLOCATE (CPLTMP)
      ENDIF

!  Start of long DO loop #7
      DO MV=IFIRST+1,MVABS

        IX=MV*IXMX
        NZERO=0
!  Start of long DO loop #8
        DO LL=1,MXLAM
          LM=LAM(2*LL-1)
          MU=LAM(2*LL)
!  Start of long DO loop #9
          DO IC=1,NSTATE
            JC  =JSTATE(IC)
            ISTC=JSTATE(IC+3*NSTATE)
            NKC =JSTATE(IC+4*NSTATE)
          DO IR=1,IC
            JR  =JSTATE(IR)
            ISTR=JSTATE(IR+3*NSTATE)
            NKR =JSTATE(IR+4*NSTATE)
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
                IF (KKR-KKC.EQ.MU) THEN
                  XCPL=XCPL+AF*GSYMTP(JC,KKC,JR,KKR,MV  ,LM,MU)
                  IF (MU.EQ.0) CYCLE
                ENDIF
                IF (KKC-KKR.NE.MU) CYCLE

C  ADJUST FOR (-1)**MU IN POTENTIAL. . .
                IF (LODD(MU)) AF=-AF
                XCPL=XCPL+AF*GSYMTP(JC,KKC,JR,KKR,MV,LM,-MU)
              ENDDO
            ENDDO
            CPL(IX)=XCPL
          ENDDO
          ENDDO
!  End of long DO loop #9
        ENDDO
!  End of long DO loop #8

      ENDDO
!  End of long DO loop #7
C  RESET IFIRST TO REFLECT HIGHEST M-VALUE STORED.
      IFIRST=MVABS
C
 4900 IF (MVABS.GT.IFIRST) GOTO 4800

C  MVABS.LE.IFIRST.  COEFFS STORED.  FILL VL() FROM CPL
      IXM=MVABS*IXMX
      IF (IVLU.GT.0) REWIND IVLU
!  Start of long DO loop #10
      DO LL=1,MXLAM
        LM=LAM(2*LL-1)
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
        DO ICOL=1,N
          I1=JSINDX(ICOL)
          J1=JSTATE(I1)
        DO IROW=1,ICOL
          I2=JSINDX(IROW)
          J2=JSTATE(I2)
          IF (I1.GT.I2) THEN
            IX12=I1*(I1-1)/2+I2
          ELSE
            IX12=I2*(I2-1)/2+I1
          ENDIF
          IXX=IXM+(LL-1)*NL12+IX12
          VL(IX)=CPL(IXX)
C  WE HAVE STORED COUPLING FOR POSITIVE MVALUES; CORRECT IF NECESSARY
C  FOR PARITY OF THRJ(J1,LM,J2,-MVAL,0,MVAL)
          IF (MVAL.LT.0 .AND. LODD(J1+J2+LM)) VL(IX)=-VL(IX)
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
          IF (IPRINT.GE.14) WRITE(6,612) MVAL,LL
        ENDIF
  612   FORMAT('  * * * NOTE.  FOR MVALUE =',I4,',  ALL COUPLING '
     1         'COEFFICIENTS ARE 0.0 FOR EXPANSION TERM',I4)
        IF (IVLU.GT.0) WRITE(IVLU) (VL(I),I=1,N*(N+1)/2)
      ENDDO
!  End of long DO loop #10

      RETURN
C
C  MV.GT.IFIRST ==> VALUES NOT STORED.  CALCULATE THEM VIA OLD CODE
 4800 IADD=20
      GOTO 3000
C
C  -------------------- OLD CODE REJOINS HERE ---------------------
C
 3000 IF (IVLU.GT.0) REWIND IVLU
      NZERO=0
!  Start of long DO loop #11
      DO LL=1,MXLAM
        LM=LAM(2*LL-1)
        MU=LAM(2*LL)
        NNZ=0
        IF (IVLU.EQ.0) THEN
          IX=LL
        ELSE
          IX=1
        ENDIF
C
!  Start of long DO loop #12
        DO IC=1,N
          JC  =JSTATE(JSINDX(IC))
          ISTC=JSTATE(JSINDX(IC)+3*NSTATE)
          NKC =JSTATE(JSINDX(IC)+4*NSTATE)
!  Start of long DO loop #13
        DO IR=1,IC
          JR  =JSTATE(JSINDX(IR))
          ISTR=JSTATE(JSINDX(IR)+3*NSTATE)
          NKR =JSTATE(JSINDX(IR)+4*NSTATE)
C
          VL(IX)=0.D0
!  Start of long DO loop #14
          DO KC=1,NKC
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
            IF (ABS(ATAU(ISTC+KC)).LE.EPS) CYCLE
            KKC=KC-JC-1

!  Start of long DO loop #15
            DO KR=1,NKR
C  SKIP IMMEDIATELY IF COEFFICIENT IS SMALL.
              IF (ABS(ATAU(ISTR+KR)).LE.EPS) CYCLE
              KKR=KR-JR-1

              AF=ATAU(ISTR+KR)*ATAU(ISTC+KC)
              IF (KKR-KKC.EQ.MU) THEN
                IF (IADD.EQ.0) THEN
                  VL(IX)=VL(IX)+
     1                   AF*FSYMTP(JC,KKC,L(IC),JR,KKR,L(IR),JTOT,LM,MU)
                ELSEIF (IADD.EQ.20) THEN
                  VL(IX)=VL(IX)+AF*GSYMTP(JC,KKC,JR,KKR,MVAL,LM,MU)
                ELSEIF (IADD.EQ.10) THEN
                  VL(IX)=VL(IX)+AF*ESYMTP(JC,KKC,JR,KKR,LM,MU)
                ENDIF
C
                IF (MU.EQ.0) CYCLE
              ENDIF
              IF (KKC-KKR.NE.MU) CYCLE

C  ADJUST FOR (-1)**MU IN POTENTIAL. . .
              AF=AF*PARSGN(MU)
              IF (IADD.EQ.0) THEN
                VL(IX)=VL(IX)+
     &                  AF*FSYMTP(JC,KKC,L(IC),JR,KKR,L(IR),JTOT,LM,-MU)

              ELSEIF (IADD.EQ.20) THEN
                VL(IX)=VL(IX)+AF*GSYMTP(JC,KKC,JR,KKR,MVAL,LM,-MU)

              ELSEIF (IADD.EQ.10) THEN
                VL(IX)=VL(IX)+AF*ESYMTP(JC,KKC,JR,KKR,LM,-MU)
C
              ENDIF
            ENDDO
!  End of long DO loop #15
          ENDDO
!  End of long DO loop #14
          IF (VL(IX).NE.0.D0) NNZ=NNZ+1
          IF (IVLU.EQ.0) THEN
            IX=IX+MXLAM
          ELSE
            IX=IX+1
          ENDIF
        ENDDO
!  End of long DO loop #13
        ENDDO
!  End of long DO loop #12
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,697) LM,MU
        ENDIF
  697   FORMAT('  * * * NOTE.  ALL COUPLING COEFFICIENTS ARE ZERO '
     1         ' FOR LM, MU = ',2I4)
        IF (IVLU.GT.0) WRITE(IVLU) (VL(I),I=1,N*(N+1)/2)
      ENDDO
!  End of long DO loop #11

      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14) THEN
        IF (L20) WRITE(6,620) 'MVAL',MVAL,NZERO
        IF (.NOT.L20) WRITE(6,620) 'JTOT',JTOT,NZERO
      ENDIF
  620 FORMAT('  * * * NOTE.  FOR ',A,' =',I4,',  ALL COUPLING ',
     1       'COEFFICIENTS ARE 0.0 FOR',I5,' POTENTIAL EXPANSION TERMS')

      RETURN
      END
C -------------------------------------------------------- END OF CPL26
