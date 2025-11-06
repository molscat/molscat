      SUBROUTINE WVBKSB(N, RBEGIN, REND, NSTEP, DR, 
     1                  PSIA, IWREC, SUMPSI, IPRINT, IPREC,
     2                  LFIRST, LREGSP)
C  Copyright (C) 2025 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ROUTINE TO PROPAGATE THE EIGENFUNCTION (WAVEFUNCTION) AT THE MATCHING
C  POINT GENERATING A VECTOR FOR IT AT EACH STEP IN THE PROPAGATION.
C  THE COUPLING MATRIX EVALUATED AT THE MIDPOINT OF EACH SECTOR
C  IS USED AS A REFERENCE POTENTIAL FOR THE SECTOR. THE MODIFIED
C  LOG DERIVATIVE PROPAGATOR OF MANOLOPOULOS IS USED. THIS VERSION
C  ADAPTED FROM DAPROP BY AE Thornley NOV 92
C
C  INTEGRATION BY SIMPSON'S RULE CORRECTED CR Le Sueur MAY 12
C  ROUTINE SIMPLIFIED TOO
C
C  FURTHER ADAPTIONS TO DEAL WITH THE POSSIBILITY OF MULTIPLE PROPAGATION
C  SEGMENTS BY CR Le Sueur NOV 2016
C
C  NOV 2016 CR Le Sueur
C  CHANGES PREPARATORY TO USING QUADRATURE WITH UNEQUAL SPACES:
C  R IS WRITTEN TO PROPAGATION FILE BY MDPROP, AND READ IN HERE, AND IS
C  ALSO WRITTEN OUT TO WAVEFUNCTION FILE.
C
C  JAN 2020 CR Le Sueur
C  NOW USES EXTENDED ALTERNATIVE SIMPSON'S RULE FOR CALCULATING
C  NORMALISATION FACTOR
C
C  SEPT 2022 CR Le Sueur
C  ROUTINE ADDED TO DO COMPOSITE SIMPSON'S RULE FOR UNEQUAL INTERVAL LENGTHS
C  BUT COMMENTED OUT AS SEEMS TO BE LESS ACCURATE THAN TRAPEZIUM RULE FOR EQUAL
C  INTERVAL LENGTHS
C
C  ON ENTRY: RBEGIN}
C            REND  } ARE LIMITS OF RADIAL PROPAGATION VARIABLE,
C            N       IS THE SIZE OF THE BASIS,
C            NSTEP   IS THE NUMBER OF STEPS USED IN THE RADIAL PROPAGATION
C            DR      IS THE STEP LENGTH
C            PSIA    IS THE WAVEFUNCTION TO START PROPAGATING FROM
C            IWREC   IS THE ADDRESS OF THE LAST W (RAB) MATRIX WRITTEN TO
C                    CHANNEL IWAVSC
C            IPREC   IS THE ADDRESS OF WHERE THE FIRST PSI IS TO BE WRITTEN
C                    TO ON CHANNEL IPSISC
C            LFIRST  IS A LOGICAL INDICATING THAT THIS IS THE FIRST
C                    SEGMENT AND SO TO WRITE OUT PSI AT THE FIRST STEP
C            LREGSP  IS A LOGICAL INDICATING THAT THE SPACING IS REGULAR
C                    SO THAT SIMPSON'S RULE CAN BE USED.  OTHERWISE THE
C                    TRAPEZIUM RULE IS USED
C  DURING:   PSIA}
C            PSIB}   ARE WORKSPACE ARRAYS USED FOR THE WAVEFUNCTION
C            RAB     IS THE PROPAGATION MATRIX
C
C  ON EXIT:  SUMPSI  CONTAINS THE ACCUMULATED NORMALISATION INTEGRALS
C
      DIMENSION PSIA(N),SUMPSI(N)
      COMMON /IOCHAN/ IPSISC,IWAVSC,IWAVE,NWVCOL,IWVSTP,IWAVEF
      LOGICAL IWAVEF
      DOUBLE PRECISION, PARAMETER :: C1=17.D0/48.D0,C2=59.D0/48.D0,
     1                               C3=43.D0/48.D0,C4=49.D0/48.D0
      ALLOCATABLE PSIB(:),RAB(:),PSIINT(:,:),HINT(:),SUMTMP(:)

      LOGICAL LFIRST,LREGSP
C
C  THIS VERSION USES A CONSTANT STEP SIZE THROUGHOUT THE
C  INTEGRATION RANGE, WITH NSTEP STEPS BETWEEN RBEGIN AND REND.
C
      H=DR
      IF (REND.GT.RBEGIN) THEN
        IDIR=1
      ELSE
        IDIR=-1
      ENDIF
      R=REND
      ROLD=R
C
      ALLOCATE (PSIB(N),RAB(N*N),PSIINT(NSTEP+1,N),HINT(NSTEP),
     1          SUMTMP(N))
      PSIINT=0.D0
      SUMTMP=0.D0
      IF (LREGSP) THEN
        E1=C1*ABS(H)
        E2=C2*ABS(H)
        E3=C3*ABS(H)
        E4=C4*ABS(H)
        E5=ABS(H)
      ENDIF

!  Start of long DO loop #1
      DO ISTEP=1,NSTEP+1
C
C  WRITE OUT PSI TO WAVEFN SCRATCH FILE
C
        IF (ISTEP.NE.1 .OR. LFIRST) THEN
          ROLD=R
          IF (.NOT.LREGSP .AND. ISTEP.NE.1) THEN
            HINT(ISTEP-1)=ABS(R-RP)
            R=RP
          ENDIF
          WRITE(IPSISC,REC=IPREC) R,PSIA
        ENDIF


        IF (ISTEP.EQ.NSTEP+1) GOTO 105

        IF (LREGSP) R=R-H
C
C  READ IN R(A,B) FROM WAVESCRATCH
C
        IWREC=IWREC-1
        READ(IWAVSC,REC=IWREC,ERR=900) RP,RAB
        IPREC=IPREC-IDIR
C
C  GENERATE PSI(B) FROM PSI(A) AND R(A,B)
C
        CALL DGEMUL(RAB,N,'N',PSIA,N,'N',PSIB,N,N,N,1)
 105    CONTINUE

        DO I=1,N
C  COLLECT TERMS FOR THE NORMALISATION CONSTANT
          PSIINT(ISTEP,I)=PSIINT(ISTEP,I)+PSIA(I)*PSIA(I)
C  COPY PSI(B) TO PSI(A) FOR NEXT SECTOR
          PSIA(I)=PSIB(I)
        ENDDO
        IF (IPRINT.GE.30) THEN
          WRITE(6,110) RP,(PSIA(K),K=1,N)
 110      FORMAT(E14.7,7X,10(E14.7,5X))
        ENDIF
      ENDDO
!  End of long DO loop #1

      DEALLOCATE (PSIB,RAB)
      IF (.NOT.LREGSP) THEN
        FAC=1D0
        IF (MOD(NSTEP,2).EQ.1) FAC=2D0
      ENDIF

!  Start of long DO loop #2
      DO ISTEP=1,(NSTEP+2)/2
        IRIGHT=ISTEP
        ILEFT=NSTEP+2-ISTEP
        IF (LREGSP) THEN
C  SET UP PARAMETERS FOR INTEGRATION BY ALTERNATIVE EXTENDED SIMPSON'S RULE
C  (SEE NUMERICAL RECIPES EQN 4.1.14)
          IF (ISTEP.EQ.1) THEN
            ALEFT=E1
            ARIGHT=E1
          ELSEIF (ISTEP.EQ.2) THEN
            ALEFT=E2
            ARIGHT=E2
          ELSEIF (ISTEP.EQ.3) THEN
            ALEFT=E3
            ARIGHT=E3
          ELSEIF (ISTEP.EQ.4) THEN
            ALEFT=E4
            ARIGHT=E4
          ELSEIF (ISTEP.EQ.5) THEN
            ALEFT=E5
            ARIGHT=E5
          ENDIF
        ELSEIF (.NOT.LREGSP) THEN
C  COMPOSITE SIMPSON'S RULE FOR IRREGULARLY SPACED DATA
C         CALL CSIMP(HINT,ARIGHT,ALEFT,IRIGHT,ILEFT,NSTEP,FAC)
C  TRAPEZIUM RULE EMPLOYED WITH UNEQUAL INTERVALS
          IF (ISTEP.EQ.1) THEN
            ARIGHT=HINT(IRIGHT)/2D0
            ALEFT=HINT(ILEFT-1)/2D0
          ELSE
            ARIGHT=(HINT(IRIGHT)+HINT(IRIGHT-1))/2D0
            ALEFT=(HINT(ILEFT)+HINT(ILEFT-1))/2D0
          ENDIF
        ENDIF
C  TO AVOID DOUBLE COUNTING OF CENTRAL POINT
        IF (ILEFT.EQ.IRIGHT) ALEFT=0.D0
        DO I=1,N
          SUMTMP(I)=SUMTMP(I)+PSIINT(IRIGHT,I)*ARIGHT
     1                       +PSIINT(ILEFT,I)*ALEFT
        ENDDO
      ENDDO
!  End of long DO loop #2
      SUMPSI=SUMPSI+SUMTMP
      DEALLOCATE (HINT,SUMTMP,PSIINT)
C
      RETURN
C
 900  WRITE(6,910) IWREC
 910  FORMAT('  *** ERROR - CRASHED ON READ OF RAB FILE - REC=',I6)
      END
C  =======================================================================
      SUBROUTINE CSIMP(HINT,ARIGHT,ALEFT,IRIGHT,ILEFT,NSTEP,FAC)
      IMPLICIT NONE
c  Written by CR Le Sueur September 2022

c  This calculates the quadrature weights for numerically integrating the
c  overlap using composite Simpson's rule for irregularly spaced data
c  (see Wiki article on Simpson's rule), but takes the mean of using that
c  rule in both directions.

c  The weights are calculated from the interval sizes around the current point:
c         -h1r=h_(c-2)-       -h2r=h_(c-1)-     -h3r=h_(c)-       -h4r=h_(c+1)
c   p_(c-2)            p_(c-1)             p_(c)           p_(c+1)
c         -h4l=h_(c-2)-       -h3l=h_(c-1)-     -h2l=h_(c)-       -h1l=h_(c+1)
c
c  Suffix right refers to summing rightwards (i.e., in the normal sense)
c  Suffix left refers to summing leftwards (i.e. backwards)

c  The number of intervals is nstep and the number of points is nstep+1.

c  The original rule has different formulae for the weights of even-numbered
c  points (starting at i=0) and odd-numbered points and is based on Simpson's
c  composite rule (with weights that go 1,4,2,4,2,..,4,1).  If there are an
c  odd number of intervals then correction terms are applied to the weights for
c  the last 3 points.

c  If nstep is even then the weights are symmetric about the middle point, so
c  it is not necessary to take the mean.

c  If nstep is odd then all weights have contributions from the formula for
c  even-numbered points and the formula for odd-numbered points, and the
c  corrections to the weights for the last 3 points mentioned above are
c  applied to both ends.
      DOUBLE PRECISION, INTENT(IN) :: HINT(NSTEP),FAC
      INTEGER, INTENT(IN) :: IRIGHT,ILEFT,NSTEP

      DOUBLE PRECISION, INTENT(OUT) :: ARIGHT,ALEFT

      DOUBLE PRECISION BRIGHT,BLEFT,CRIGHT,CLEFT,DRIGHT,DLEFT,
     1                 H1R,H2R,H3R,H4R,H1L,H2L,H3L,H4L

      H3R=HINT(IRIGHT)
      H4R=HINT(IRIGHT+1)
      H3L=HINT(ILEFT-1)
      H4L=HINT(ILEFT-2)
      DLEFT=0D0
      DRIGHT=0D0
      IF (IRIGHT.GT.1) THEN
        H2R=HINT(IRIGHT-1)
        H2L=HINT(ILEFT)
      ENDIF
      IF (IRIGHT.GT.2) THEN
        H1R=HINT(IRIGHT-2)
        H1L=HINT(ILEFT+1)
      ENDIF
      IF (MOD(IRIGHT,2).EQ.1 .OR. MOD(ILEFT,2).EQ.1) THEN
        ARIGHT=(H3R+H4R)*(2D0-H4R/H3R)/6D0
        ALEFT=(H3L+H4L)*(2D0-H4L/H3L)/6D0
      ELSE
        ARIGHT=0D0
        ALEFT=0D0
      ENDIF
      IF (IRIGHT.GT.1 .AND.
     1    (MOD(IRIGHT,2).EQ.0 .OR. MOD(ILEFT,2).EQ.0)) THEN
        BRIGHT=(H2R+H3R)**3/(6D0*H2R*H3R)
        BLEFT=(H2L+H3L)**3/(6D0*H2L*H3L)
      ELSE
        BRIGHT=0D0
        BLEFT=0D0
      ENDIF
      IF (IRIGHT.GT.2 .AND.
     1    (MOD(IRIGHT,2).EQ.1 .OR. MOD(ILEFT,2).EQ.1)) THEN
        CRIGHT=(H1R+H2R)*(2D0-H1R/H2R)/6D0
        CLEFT=(H1L+H2L)*(2D0-H1L/H2L)/6D0
      ELSE
        CRIGHT=0D0
        CLEFT=0D0
      ENDIF
      IF (MOD(NSTEP,2).EQ.1) THEN
        IF (IRIGHT.EQ.1) THEN
          DLEFT=(2D0*H3L**2+3D0*H3L*H4L)/(6D0*(H3L+H4L))
          DRIGHT=(2D0*H3R**2+3D0*H3R*H4R)/(6D0*(H3R+H4R))
        ELSEIF (IRIGHT.EQ.2) THEN
          DLEFT=(H3L**2+3D0*H3L*H4L)/(6D0*H4L)
          DRIGHT=(H3R**2+3D0*H3R*H4R)/(6D0*H4R)
        ELSEIF (IRIGHT.EQ.3) THEN
          DLEFT=-H3L**3/(6D0*H4L*(H3L+H4L))
          DRIGHT=-H3R**3/(6D0*H4R*(H3R+H4R))
        ENDIF
      ENDIF


      ARIGHT=(ARIGHT+BRIGHT+CRIGHT+DRIGHT)/FAC
      ALEFT=(ALEFT+BLEFT+CLEFT+DLEFT)/FAC
        
      RETURN
      END
