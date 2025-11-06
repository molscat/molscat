      SUBROUTINE CPLOUT(IV,VL,N,NHAM)
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS ROUTINE PRINTS OUT THE COUPLING MATRIX ELEMENTS.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IV(*), VL(*)
      PARAMETER (NCOLS=10)
      CHARACTER(30) F11,F20
      CHARACTER(2) CNCOLS
C
      COMMON /VLFLAG/ IVLFL

      WRITE(CNCOLS,'(I2)') NCOLS
      F11='(/'//CNCOLS//'(9X,I3))'
      WRITE(6,602) NHAM
  602 FORMAT(/'  COUPLING MATRIX ELEMENTS BETWEEN CHANNELS FOR',I4,
     1       ' EXPANSION TERMS.')

!  Start of long IF block #1
      IF (IVLFL.GT.0) THEN
        IMAX=0
        DO I=1,N
        DO J=1,I
          IMIN=IMAX+1
          IMAX=IMAX+NHAM
          WRITE(6,600) I,J
  600     FORMAT(/'  FOR CHANNEL ',I3,'  TO  CHANNEL',I4)
          WRITE(6,601) (IV(IJ),VL(IJ),IJ=IMIN,IMAX)
  601     FORMAT(' ',7(I3,1X,F12.5))
        ENDDO
        ENDDO
      ELSE
!  Start of long DO loop #1
        DO LL=1,NHAM
          F20='(I4,'//CNCOLS//'(1X,F11.5))'
          VLMAX=0.D0
          DO K=1,N*(N+1)/2
            VLMAX=MAX(VLMAX,ABS(VL(NHAM*(K-1)+LL)))
          ENDDO
          IF (VLMAX.EQ.0D0) THEN
            WRITE(6,*) ' ALL ELEMENTS FOR POTENTIAL TERM',LL,' ARE 0'
            WRITE(6,16)
            CYCLE
          ENDIF
          WRITE(6,*) ' POTENTIAL TERM',LL
          ITIME=(N+NCOLS-1)/NCOLS
          IF (VLMAX.LT.1D-2) F20='(0P,I4,1P,'//CNCOLS//'(1X,G11.4))'
          KF=0
          KF1=KF
          WRITE(6,16)
   16     FORMAT(/)
          J1=1
          DO KK=1,ITIME
            KS=(KK-1)*NCOLS + 1
            KF=KS + NCOLS-1
            KF1=KF
            IF (KK.EQ.ITIME) KF1=MIN(N,KF)
            WRITE(6,FMT=F11) (I,I=KS,KF1)
            J1=KS
            DO JJ=J1,N
              KF1=MIN(JJ,KF)
              II1=JJ*(JJ-1)/2
              WRITE(6,FMT=F20) JJ,(VL(NHAM*(II1+I-1)+LL),I=KS,KF1)
            ENDDO
            WRITE(6,16)
          ENDDO
        ENDDO
!  End of long DO loop #1
      ENDIF
!  End of long IF block #1
      RETURN
      END
