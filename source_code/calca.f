      SUBROUTINE CALCA(NOPEN,NBASIS,L,WVEC,SREAL,SIMAG,AWVMAX,
     1                 SCLEN,ICHAN,RUNAME,LRPOW,IPRINT)
C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT NONE
C  ROUTINE TO CALCULATE SCATTERING LENGTHS/VOLUMES ADDED SEPT 06
C  THIS USES EXACT EXPRESSION A=(1/IK)*(1-S)/(1+S)
C  FROM J. M. HUTSON, NEW J PHYS 9, 152 (2007),
C  GENERALISED TO A**M=(1/I*K**M)*(1-S)/(1+S)
C  WHERE M = MIN(LRPOW,2*LL+1)
C  AND THE POTENTIAL IS 1/R**LRPOW AT LONG RANGE
C  THE POWERS USED ARE DISCUSSED BY SADEGHPOUR ET AL.,
C  J PHYS B 33, R93 (2000).
C  ROUTINE SPINIT SETS LRPOW=6 BUT ANOTHER VALUE CAN BE SET IF NEEDED

      INTEGER,          INTENT(IN)  :: NOPEN,NBASIS(1),L(1),ICHAN,
     1                                 LRPOW,IPRINT
      DOUBLE PRECISION, INTENT(IN)  :: WVEC(1),SREAL(1),SIMAG(1),AWVMAX
      CHARACTER(10),    INTENT(IN)  :: RUNAME
      DOUBLE COMPLEX,   INTENT(OUT) :: SCLEN

      INTEGER IROW,LL,IJ,NLOW,I
      DOUBLE PRECISION DD,PHASE,AREALD,AREALI,AIMAGR
      DOUBLE COMPLEX SCATLN

      NLOW=0
      DO IROW=1,NOPEN
        DD=WVEC(NBASIS(IROW))
        IF (DD.LT.AWVMAX) NLOW=NLOW+1
      ENDDO

      IF (NLOW.GT.0 .AND. IPRINT.GE.6)
     1  WRITE(6,605) (ADJUSTL(RUNAME),I=1,3)
  605 FORMAT(/'  K-DEPENDENT SCATTERING LENGTHS/VOLUMES/HYPERVOLUMES',
     1       ' FOR CHANNELS WITH LOW KINETIC ENERGY'/
     2       '  CHAN   L POW',5X,'WVEC*',A10,7X,'RE(A)/',A10,6X,
     3       'IM(A)/',A)
  604 FORMAT(2I5,I4,4(ES21.13E3,2X))

      DO 2000 IROW=1,NOPEN
        DD=WVEC(NBASIS(IROW))
        LL=L(NBASIS(IROW))
        IJ=IROW+NOPEN*(IROW-1)
        IF (SREAL(IJ).LE.-1.D0) THEN
C  pgf90 COMPILER DOES NOT RECOGNISE CMPLX AS GENERIC
          SCATLN=DCMPLX(1.D30,0.D0)
        ELSE
          SCATLN=DCMPLX(1.D0-SREAL(IJ),-SIMAG(IJ))
     1         /(DCMPLX(1.D0+SREAL(IJ),+SIMAG(IJ))
     2           *DCMPLX(0.D0,DD**(MIN(LRPOW-2,LL+LL+1))))
        ENDIF
        IF (ICHAN.EQ.IROW) SCLEN=SCATLN
        IF (DD.GE.AWVMAX) GOTO 2000
        IF (IPRINT.GE.6) WRITE(6,604) IROW,LL,MIN(LRPOW-2,LL+LL+1),DD,
     1                                DBLE(SCATLN),IMAG(SCATLN)
        IF (IPRINT.LT.26) GOTO 2000

C  THE LIMITING FORMS THAT APPEAR IN MANY PAPERS WORK ONLY IF
C  THE PHASE SHIFT DELTA=-K*A IS SMALL, WHICH IS NOT TRUE NEAR
C  A RESONANCE EVEN AT VERY LOW ENERGY
        PHASE=ATAN2(SIMAG(IJ),SREAL(IJ))
        AREALD=-0.5D0*PHASE/DD**MIN(LRPOW-2,LL+LL+1)
        AREALI=-0.5D0*SIMAG(IJ)/DD**MIN(LRPOW-2,LL+LL+1)
        AIMAGR= 0.5D0*(1.D0-SREAL(IJ))/DD**MIN(LRPOW-2,LL+LL+1)
        WRITE(6,604) IROW,LL,MIN(LRPOW-2,LL+LL+1),DD,AREALI,AIMAGR,
     1               AREALD
 2000 CONTINUE

      RETURN
      END
