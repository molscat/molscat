      SUBROUTINE SIG6(NSTATE,JSTATE,LI,LF,SIG,S,IMSG,QL,IXQL,
     2                NIXQL,NQL,LM,LMAX)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  ROUTINE TO EVALUATE SIG(J,TAU->J',TAU') FROM IOS Q(L,M1,M2)
C  VALUE FOR LEVEL LI TO LF RETURNED IN SIG
C
      USE pair_state, ONLY: ATAU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JSTATE(4,NSTATE),IXQL(NIXQL,NQL),LM(3,LMAX)
      DIMENSION QL(*)
      CHARACTER(1) S,STAR
      ALLOCATABLE COEFFR(:),COEFFI(:)
C
      DATA STAR/'*'/
      DATA EPS/1.D-8/
C
C  STATEMENT FUNCTION FOR INDEX M1.GE.M2, M STARTING AT ZERO.
      IX(M1,M2)=M1*(M1+1)/2+M2+1
C
      SIG=0.D0
      JI=JSTATE(1,LI)
      XJI=JI
      NKI=2*JI+1
      ISTAI=JSTATE(4,LI)
      JF=JSTATE(1,LF)
      XJF=JF
      NKF=2*JF+1
      ISTAF=JSTATE(4,LF)
      LMN=ABS(JI-JF)
      LMX=JI+JF

!  Start of long DO loop #1
      DO L=LMN,LMX
        XL=L
        MMAX=L
        IXMX=IX(MMAX,MMAX)

C  SET STORAGE POINTERS AND ZERO TEMP STORAGE.
        ALLOCATE (COEFFR(IXMX),COEFFI(IXMX))
        DO II=1,IXMX
          COEFFR(II)=0.D0
          COEFFI(II)=0.D0
        ENDDO
C  -------------LOOP OVER IPI,IPF  IQI,IQF -----------
        IPI=-JI-1
!  Start of long DO loop #2
        DO IIPI=1,NKI
          IPI=IPI+1
          API=ATAU(ISTAI+IIPI)
          IF (ABS(API).LE.EPS) CYCLE

          PI=IPI
          IPF=-JF-1
!  Start of long DO loop #3
          DO IIPF=1,NKF
            IPF=IPF+1
            APF=ATAU(ISTAF+IIPF)
            IF (ABS(APF).LE.EPS) CYCLE

            PF=IPF
            IF (ABS(IPI-IPF).GT.MMAX) CYCLE

            IQI=-JI-1
!  Start of long DO loop #4
            DO IIQI=1,NKI
              IQI=IQI+1
              AQI=ATAU(ISTAI+IIQI)
              IF (ABS(AQI).LE.EPS) CYCLE

              QI=IQI
              IQF=-JF-1
!  Start of long DO loop #5
              DO IIQF=1,NKF
                IQF=IQF+1
                AQF=ATAU(ISTAF+IIQF)
                IF (ABS(AQF).LE.EPS) CYCLE

                QF=IQF
                IF (ABS(IQI-IQF).GT.MMAX) CYCLE

C  CALCULATE FACTOR
                TJ1 = THRJ(XJI,XL,XJF,-PI,PI-PF,PF)
                IF (ABS(TJ1).LE.EPS) CYCLE

                TJ2 = THRJ(XJI,XL,XJF,-QI,QI-QF,QF)
                IF (ABS(TJ2).LE.EPS) CYCLE

                FACT=API*AQI*APF*AQF  *TJ1*TJ2
C  RECALCULATE MP,MQ AS THEY MIGHT HAVE BEEN SWAPPED IN LAST LOOP.
                MP=IPI-IPF
                MQ=IQI-IQF
                SIGNR=1.D0
                SIGNI=1.D0
                IF (MP.GE.0) GOTO 1401

                P=PARSGN(MP)
                SIGNR=P*SIGNR
                SIGNI=P*SIGNI
                MP=ABS(MP)
 1401           IF (MQ.GE.0) GOTO 1402

                P=PARSGN(MQ)
                SIGNR=P*SIGNR
                SIGNI=P*SIGNI
                MQ=ABS(MQ)
 1402           IF (MP.GE.MQ) GOTO 1403

                SIGNI=-SIGNI
                MT=MP
                MP=MQ
                MQ=MT
 1403           INDX=IX(MP,MQ)
                IF (MP.EQ.MQ) SIGNI=0.D0
                COEFFR(INDX)=COEFFR(INDX)+SIGNR*FACT
                COEFFI(INDX)=COEFFI(INDX)+SIGNI*FACT
              ENDDO
!  End of long DO loop #5
            ENDDO
!  End of long DO loop #4
C  ----------  THIS ENDS LOOP OVER IQI,IQF
          ENDDO
!  End of long DO loop #3
        ENDDO
!  End of long DO loop #2
C  ----------  THIS ENDS LOOP OVER IPI,IPF

C  MATCH CONTRIBUTING (I.E., NON-ZERO) CR WITH QL VALUES
        IZERO=0
        INDX=0
        DO MP=IZERO,MMAX
        DO MQ=IZERO,MP
          INDX=INDX+1
C  N.B. IMAGINARY PART SHOULD VANISH; ERROR MESSAGE IF ANY SURVIVE.
          IF (ABS(COEFFI(INDX)).LE.EPS) GOTO 1501

          WRITE(6,694) L,MP,MQ,COEFFI(INDX),LI,LF
  694     FORMAT(/' *** ERROR.  NON-ZERO IMAGINARY COEFF QL(',
     &            3I4,' ) =',F12.6,'  FOR LI,LF =',2I4)
 1501     IF (ABS(COEFFR(INDX)).LE.EPS) CYCLE

C  CALL IXQLF TO GET INDEX OF L,MP,MQ IN QL
C  AND ACCUMULATE IN CROSS SECTION
          CALL IXQLF(LM,LMAX,L,MP,MQ,1,INDEX,IXQL,NIXQL,NQL)
C  N.B. 6TH ARG (1) ASKS FOR REAL PART; SHOULD WORK OK FOR MP.EQ.MQ
          IF (INDEX.GT.0) GOTO 1502

          IF (INDEX.EQ.-1) CYCLE

          S=STAR
          IMSG=1
          CYCLE

 1502     SIG=SIG + COEFFR(INDX)*QL(INDEX)
  602     FORMAT(2X,'I/F=',2I3,'    QL(',3I3,' )  COEFF/QL =',2F10.5)
        ENDDO
        ENDDO

C  RECOVER TEMPORARY STORAGE ...
        DEALLOCATE (COEFFR,COEFFI)
      ENDDO
!  End of long DO loop #1
C  ----------  THIS ENDS LOOP OVER L - VALUES

C  MULTIPLY FINALLY BY 2*JF+1
      SIG = SIG * (2*JF+1)
      RETURN
      END
