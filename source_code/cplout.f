      SUBROUTINE CPLOUT(IV,VL,N,NHAM)
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS ROUTINE PRINTS OUT THE COUPLING MATRIX ELEMENTS.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IV(1), VL(1)
      PARAMETER (NCOLS=10)
      CHARACTER(20) F11,F20
      CHARACTER(2) CNCOLS
C
      COMMON /VLFLAG/ IVLFL

      WRITE(CNCOLS,'(I2)') NCOLS
      F11='(/'//CNCOLS//'(9X,I3))'
      F20='(I4,'//CNCOLS//'(1X,F11.5))'
      WRITE(6,602) NHAM
  602 FORMAT(/'  COUPLING MATRIX ELEMENTS BETWEEN CHANNELS FOR',I4,
     1       ' EXPANSION TERMS.')

      IF (IVLFL.GT.0) THEN
        IMAX=0
        DO 1000 I=1,N
        DO 1000 J=1,I
          IMIN=IMAX+1
          IMAX=IMAX+NHAM
          WRITE(6,600) I,J
  600     FORMAT(/'  FOR CHANNEL ',I3,'  TO  CHANNEL',I4)
          WRITE(6,601) (IV(IJ),VL(IJ),IJ=IMIN,IMAX)
  601     FORMAT(' ',7(I3,1X,F12.5))
 1000   CONTINUE
      ELSE
        DO 900 LL=1,NHAM
          WRITE(6,*) ' POTENTIAL TERM',LL
          ITIME=(N+NCOLS-1)/NCOLS
          KF=0
          KF1=KF
          WRITE(6,16)
   16     FORMAT(/)
          J1=1
          DO 800 KK=1,ITIME
            KS=(KK-1)*NCOLS + 1
            KF=KS + NCOLS-1
            KF1=KF
            IF (KK.EQ.ITIME) KF1=MIN(N,KF)
            WRITE(6,FMT=F11) (I,I=KS,KF1)
   11       FORMAT(/10(9X,I3))
            J1=KS
            DO 700 JJ=J1,N
              KF1=MIN(JJ,KF)
              II1=JJ*(JJ-1)/2
              WRITE(6,FMT=F20) JJ,(VL(NHAM*(II1+I-1)+LL),I=KS,KF1)
   20         FORMAT(I4,10(1X,F11.5))
  700       CONTINUE
            WRITE(6,16)
  800     CONTINUE
  900   CONTINUE
      ENDIF
      RETURN
      END
