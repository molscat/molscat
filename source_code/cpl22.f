      SUBROUTINE CPL22(N,MXLAM,NHAM,LAM,NSTATE,JSTATE,JSINDX,MVALUE,IV,
     1                 VL,IPRINT,LFIRST)
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  CS COUPLING MATRIX FOR VIBRATING ROTOR-ATOM (ITYPE=22)
C  SEE (FOR EXAMPLE) EQN 30 OF MCGUIRE AND KOURI JCP (1974) 60 2488
C  S Green (MAR 94) USES IV(), I.E., IVLFL=1
C                   SAVES COUPLING MATRIX FOR MV=0,MX  IN ALLOCATABLE
C                   ARRAY
C                   USES J3J000 ROUTINE AS PER JMH CPL21 CODE
C                   STORES ON J OR NSTATE, DEPENDING ON WHICH IS SMALLER
C
C  CRLS 10-05-22: SAVES COUPLING COEFFICIENTS USING ALLOCATABLE ARRAYS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE NL12,IXMX,ISTART,IFIRST,LOGIX,JTOP,CPL
C  SPECIFICATIONS FOR ARGUMENTS
      DIMENSION LAM(3,MXLAM),JSTATE(NSTATE),JSINDX(N),VL(*),IV(*)
      INTEGER IPRINT
      LOGICAL LFIRST
      ALLOCATABLE CPL(:),CPLTMP(:)
C
      LOGICAL LODD,LOGIX
      PARAMETER (Z0=0.D0)
C
C
C  STATEMENT FUNCTION DEFINITIONS
      Z(I)=DBLE(I+I+1)
      LODD(I)=I-2*(I/2).NE.0
C
C  IF LFIRST IS TRUE (FIRST CALL), DO SOME INITIALIZATION
      IF (LFIRST) THEN
        IFIRST=-1
        LFIRST=.FALSE.
        IF (ALLOCATED(CPL)) DEALLOCATE (CPL)
      ENDIF
C
      XM=MVALUE
      PM=1.D0
      IF (LODD(MVALUE)) PM=-1.D0
C
      IF (IFIRST.GT.-1) GOTO 1000

C  FIRST TIME THROUGH SET UP SOME STORAGE POINTERS.
C  LOGIX=.TRUE. IF JTOP IS SMALLER THAN NSTATE (SO STORE ON J)
      JTOP=0
      DO I=1,NSTATE
        JTOP=MAX(JTOP,JSTATE(I))
      ENDDO
      LOGIX=JTOP.LT.NSTATE
      IF (LOGIX) THEN
        NL12=(JTOP+1)*(JTOP+2)/2
      ELSE
        NL12=NSTATE*(NSTATE+1)/2
      ENDIF
      IXMX=NL12*NHAM
      ALLOCATE (CPL(IXMX))
C
 1000 MVABS=ABS(MVALUE)
C  SEE IF VALUES ARE STORED FOR THIS HIGH AN MVALUE
C  IF NOT, TRY TO STORE THEM IN CPL().
      IF (MVABS.LE.IFIRST) GOTO 2000

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

!  Start of long DO loop #1
      DO MV=IFIRST+1,MVABS

        IX=MV*IXMX
C
C  THIS SECTION OF CODE CALCULATES THE COUPLING MATRIX ELEMENTS USING
C  THE SUBROUTINE J3J000 WHICH CALCULATES A WHOLE SET OF 3J
C  COEFFICIENTS (J1 J2 J3).  THE MATRIX ELEMENTS ARE STORED IN CPLTMP
C               ( 0  0  0)
        PMV=1.D0
        ALLOCATE(CPLTMP(2*NHAM-1))
        IF (LODD(MV)) PMV=-1.D0
C  CODE BELOW FROM V12 (DEC 94) CPL21 CODE
C  EXCEPT LIMIT ON IL LOOP AND VALUE OF LM
        IF (LOGIX) THEN
          ITOP=JTOP+1
        ELSE
          ITOP=NSTATE
        ENDIF
!  Start of long DO loop #2
        DO IL=1,NHAM
          LM=IL-1
          JSAV=-1
!  Start of long DO loop #3
          DO I1=1,ITOP
            IF (LOGIX) THEN
              J1=I1-1
            ELSE
              J1=JSTATE(I1)
            ENDIF
            IF (J1.NE.JSAV) THEN
              CALL J3J000(DBLE(J1),DBLE(LM),IVALJ,CPLTMP,XJMIN)
              JMIN=ABS(J1-LM)
              JMAX=J1+LM
              JSAV=J1
            ENDIF
          DO I2=1,I1
            IF (LOGIX) THEN
              J2=I2-1
            ELSE
              J2=JSTATE(I2)
            ENDIF
            IX=IX+1
            IF (J2.LT.JMIN .OR. J2.GT.JMAX .OR. J1.LT.MV .OR.
     1          J2.LT.MV .OR. LODD(J2+JMAX)) THEN
              CPL(IX)=0.D0
            ELSE
              INDJ=(J2-JMIN)/2+1
              IF (MV.EQ.0) THEN
                CPL(IX)=PMV*SQRT(Z(J1)*Z(J2))*CPLTMP(INDJ)**2
              ELSE
                CPL(IX)=PMV*SQRT(Z(J1)*Z(J2))*CPLTMP(INDJ)*
     1                  THRJ(DBLE(J1),DBLE(LM),DBLE(J2),
     2                       -DBLE(MV),0.D0,DBLE(MV))
              ENDIF
            ENDIF
          ENDDO
          ENDDO
!  End of long DO loop #3
        ENDDO
!  End of long DO loop #2
        DEALLOCATE (CPLTMP)
C
C  RESET IFIRST TO REFLECT HIGHEST M-VALUE STORED.
      ENDDO
!  End of long DO loop #1
      IFIRST=MVABS
C
C  THIS SECTION OF CODE TRANSFERS THE COUPLING MATRIX ELEMENTS FROM THE
C  ARRAY CPL WHERE THEY WERE PREVIOUSLY STORED INTO THE VL ARRAY, CHANGING
C  SIGNS AS NECESSARY.
C
C  START BY ZEROING VL, IV ARRAYS
 2000 NTOP=NHAM*N*(N+1)/2
      DO I=1,NTOP
        VL(I)=0.D0
        IV(I)=0
      ENDDO

      IXM=MVABS*IXMX
      NZERO=0
!  Start of long DO loop #4
      DO LL=1,MXLAM
        NNZ=0
        LM=LAM(1,LL)
        XLM=LM
        NV=LAM(2,LL)
        NV1=LAM(3,LL)
C  ICR COUNTS ICOL,IROW LOOP; NEEDED FOR IXVL (VL INDEX)
        ICR=0
!  Start of long DO loop #5
        DO ICOL=1,N
          I1=JSINDX(ICOL)
          JCOL =JSTATE(I1)
          XJCOL=JCOL
          NVC=JSTATE(NSTATE+I1)
!  Start of long DO loop #6
        DO IROW=1,ICOL
          I2=JSINDX(IROW)
          JROW =JSTATE(I2)
          XJROW=JROW
          NVR=JSTATE(NSTATE+I2)
          ICR=ICR+1
!  Start of long IF block #1
          IF ((NV.EQ.NVC .AND. NV1.EQ.NVR) .OR.
     1        (NV.EQ.NVR .AND. NV1.EQ.NVC)) THEN
C  FIRST GET INDEX IN VL, IV
            IXVL=(ICR-1)*NHAM+LM+1
            IV(IXVL)=LL
C           IF (NOMEM) THEN
C
C  THIS SECTION OF CODE CALCULATES COUPLING MATRIX ELEMENTS USING
C  FUNCTION THREEJ WHICH PRODUCES A SINGLE 3J COEFFICIENT (J1 J2 J3)
C                                                         ( 0  0  0).
C             VL(IXVL)=PM*SQRT(Z(JROW)*Z(JCOL))*THREEJ(JROW,LM,JCOL)*
C    &                                   THRJ(XJROW,XLM,XJCOL,-XM,Z0,XM)
C           ELSE
C  THEN GET INDEX OF STORED COUPLING COEFFICIENT, DEPENDING ON LOGIX
              IF (LOGIX) THEN
                IF (J1.GT.J2) THEN
                  IX12=(J1+1)*J1/2+J2+1
                ELSE
                  IX12=(J2+1)*J2/2+J1+1
                ENDIF
              ELSE
                IF (I1.GT.I2) THEN
                  IX12=I1*(I1-1)/2+I2
                ELSE
                  IX12=I2*(I2-1)/2+I1
                ENDIF
              ENDIF
              IX=IXM+LM*NL12+IX12
              VL(IXVL)=CPL(IX)
              IF (MVALUE.LT.0 .AND. LODD(JCOL+JROW+LM))
     1          VL(IXVL)=-VL(IXVL)
C           ENDIF
            IF (VL(IXVL).NE.0.D0) NNZ=NNZ+1
C  WE HAVE STORED COUPLING FOR POSITIVE MVALUES; CORRECT IF NEC
          ENDIF
!  End of long IF block #1
        ENDDO
!  End of long DO loop #6
        ENDDO
!  End of long DO loop #5
        IF (NNZ.LE.0) THEN
          NZERO=NZERO+1
          IF (IPRINT.GE.14) WRITE(6,612) MVALUE,LL
        ENDIF
  612   FORMAT('  * * * NOTE.  FOR MVALUE =',I4,',  ALL COUPLING '
     1         'COEFFICIENTS ARE 0.0 FOR EXPANSION TERM',I4)
      ENDDO
!  End of long DO loop #4

      IF (NZERO.GT.0 .AND. IPRINT.GE.10 .AND. IPRINT.LT.14)
     1  WRITE(6,620) MVALUE,NZERO
  620 FORMAT('  * * * NOTE.  FOR MVALUE =',I4,',  ALL COUPLING ',
     1       'COEFFICIENTS ARE 0.0 FOR',I5,' POTENTIAL EXPANSION TERMS')

      RETURN
      END
