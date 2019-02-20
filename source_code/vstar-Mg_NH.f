      SUBROUTINE  VINIT(II,RM,RESULT)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C -------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DOUBLE PRECISION, ALLOCATABLE:: ANGLS(:),ALPHAS(:,:),BETA(:)
      DOUBLE PRECISION, ALLOCATABLE:: POL(:,:), POLN(:,:),WT(:)
      DOUBLE PRECISION, ALLOCATABLE:: PTST(:)
      DOUBLE PRECISION, ALLOCATABLE:: RPT(:,:), VPT(:,:), RR(:),P(:)
      INTEGER, ALLOCATABLE:: NRA(:)
      LOGICAL LFIRST
      CHARACTER(80) FILNAM

      DATA LFIRST /.TRUE./
      DATA FILNAM /'data/pot-Mg_NH.data'/

C
C  MGNH - POTENTIAL OBTAINED IN AVQZ (AVTZ ON CR) BASIS
C  SET WITH  MIDBOND

      IF (II.NE.1) GOTO 15

C  READ THE FILE data/pot-Mg_NH.data
C  THE STRUCTURE OF FILE IS FOLLOWING:
C  NANGLE  ICASEM: NUMBER OF ANGULAR CUTS,
C          ICASE = IF THE MORPHING IS SWITCHED ON
C  ANGLE1
C      XX1  EINT1 EINT2  EEX  EIND  EDISP
C      ...
C      XXNR1   EINT1 EINT2  EEX  EIND  EDISP
C  ANGLE2
C  ETC...

      IF (.NOT.LFIRST) DEALLOCATE(ANGLS,NRA,VPT,ALPHAS,RPT,P,RR,POL,
     1                            BETA,WT)
      LFIRST=.FALSE.
      OPEN(9,FILE=FILNAM,FORM='FORMATTED',STATUS='OLD',ERR=2010)
      READ(9,*) NANGS
      ALLOCATE(ANGLS(NANGS))
      ALLOCATE(NRA(NANGS))
C  DETERMINE MAXIMUM NRA(I) TO ALLOCATE PROPERLY ALPHA
C  AND VPT ARRAYS
      NRAMAX = 0
      DO I=1,NANGS
        READ(9,*) ANG_I,NRA_I
        ANGLS(I) = ANG_I
        NRA(I) = NRA_I
        DO J=1,NRA(I)
          READ(9,*)
          NRAMAX = MAX(NRAMAX,J)
        ENDDO
      ENDDO

      REWIND(9)
      READ(9,*)

      ALLOCATE(VPT(NRAMAX,NANGS))
      ALLOCATE(ALPHAS(NRAMAX,NANGS))
      ALLOCATE(RPT(NRAMAX,NANGS))
      ALLOCATE(P(NANGS))
      ALLOCATE(RR(NANGS))

      DO I=1,NANGS
        READ(9,*) ANGLS(I), NRA(I)
        DO J=1,NRA(I)
           READ(9,*)  RPT0, EINT
           RPT(J,I) = RPT0
           VPT(J,I) = EINT
        ENDDO
      ENDDO
      CLOSE(9)

      WRITE(6,*) ' NUMBER OF ANGLES: ',NANGS
C     WRITE(6,'(20F8.2)') (ANGLS(I),I=1,NANGS)
C     WRITE(6,*) ' NUMBER OF CALCULATED POINTS FOR EACH ANGLE: '
C     WRITE(6,'(20I8)') (NRA(I),I=1,NANGS)
      LMAX = NANGS-1

C  INITIALIZE ARRAY WITH PN FOR EACH ANGLE
C  CALCULATE ASSOCIATED LEGENDRE POLYNOMIALS
C  IT WILL BE NEEDED LATER
C
      ALLOCATE(POL(NANGS,NANGS))
      PI = ACOS(-1.0D0)
      LAMMAX = LMAX
      DO IA=1,NANGS
        COSTH = COS(ANGLS(IA)*PI/180.0D0)
        CALL CALCPL(NANGS, COSTH, POL(1,IA))
      ENDDO

      NP = 3
      MP = 5
      LP = 1
      NFIX = 0
      ASYM = 0.0D0
      BETA1= 3.D0/56.D0
      BETA2=-4.D0*BETA1/3.D0
      BETA3= 7.D0*BETA1/15.D0
      ALLOCATE(BETA(NP))
      BETA(1) =  BETA1
      BETA(2) =  BETA2
      BETA(3) =  BETA3

C  CALCULATE RKHS EXPANSION COEFFICIENTS 'ALPHA'
C  NOTE THAT FOR THIS SURFACE WE USE N=3, M=5
C  THE ASYMPTOTICS IS NOT FORCED

      DO I=1,NANGS
C       WRITE(*,*) (RPT(JJ,I),JJ=1,NRA(I))
        CALL RK_INIT(NRA(I),RPT(1,I),VPT(1,I),ALPHAS(1,I),BETA,
     O               NP, MP, LP, NFIX, ASYM)
      ENDDO

      IP=0
      MXLAM=LMAX
      ALLOCATE(WT(NANGS))
      DO 70 IA=1,NANGS
        WT(IA)=2.D0/(DBLE(NANGS*(NANGS-1))
     1         *POL(NANGS,IA)*POL(NANGS,IA))
   70 CONTINUE

  15  RETURN

C -------------------------------------------------------
      ENTRY VSTAR(II,X,TOTAL)
C -------------------------------------------------------
      IF (II.GT.LMAX) GOTO 17

      TOTAL = 0.0D0

      CALL GETRAD(X, RR, VPT, NANGS, NRA, RPT, WT,
     O            ALPHAS, BETA, NP, MP, LP, NRAMAX )

C  RR - OUTPUT - RR FOR EACH ANGLE AT GIVEN X
C  INTEGRATION - LOOP OVER QUADRATURES WT TIMES VALUE

      DO IA=1, NANGS
        TOTAL = TOTAL + POL(II,IA)*RR(IA)*WT(IA)
      ENDDO
      TOTAL = TOTAL * (2.0D0*(II-1) + 1.0D0 )/2.0D0

   16 RETURN

   17 STOP 'WRONG LAMBDA IN VSTAR ROUTINE '
      RETURN
C --------------------------------------------------------
C  DUMMY ROUTINES
C --------------------------------------------------------
      ENTRY VSTAR1(II,X,TOTAL)
C --------------------------------------------------------
C --------------------------------------------------------
      ENTRY VSTAR2(II,X,TOTAL)
C --------------------------------------------------------
      WRITE(6,*) ' VSTAR: DERIVATIVES NOT IMPLEMENTED'
      RETURN
 2010 WRITE(6,*) ' *** ERROR: '//TRIM(FILNAM)//' NOT FOUND'
      STOP
      END
C ----------------------------------------------------
      SUBROUTINE RK_INIT(NR,RPT,VPT,ALPHA,BETA,N,M,L,NFIX,ASYM)
C ----------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IPIV(:),WORK(:),Q(:,:),RPT(NR),VPT(NR)
      DIMENSION BETA(1),ALPHA(NR)
      DIMENSION WORK2(:)
      ALLOCATABLE IPIV,WORK,Q,WORK2

      ALLOCATE(IPIV(NR))
      ALLOCATE(WORK(NR))
      ALLOCATE(WORK2(NR))
      ALLOCATE(Q(NR,NR))

      DO I=1,NR
        DO J=1,NR
          RS = MIN( RPT(I),RPT(J))
          RL = MAX( RPT(I),RPT(J))

          Q(I,J)=0.0D0
          DO K=0,N-1
            Q(I,J)=Q(I,J)+BETA(K+1)*(RS/RL)**K
          ENDDO
          Q(I,J)=Q(I,J)*(RL)**(-(M+1))
        ENDDO
      ENDDO

      CALL DCOPY (NR,VPT,1,WORK2,1)
      CALL DGESV (NR, 1, Q, NR, IPIV, VPT, NR, INFO)
      CALL DCOPY (NR,VPT,1,ALPHA,1)
      CALL DCOPY (NR,WORK2,1,VPT,1)
      DEALLOCATE(IPIV,WORK,WORK2,Q)

      END
C ----------------------------------------------------
      DOUBLE PRECISION FUNCTION VRKHS(R,N,M,L,NR,BETA,RPT,ALPHA)
C ----------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BETA(1),RPT(NR),ALPHA(NR)

      V = 0.0D0
      DO I=1,NR
        RS = MIN(R,RPT(I))
        RL = MAX(R,RPT(I))
        Q  = 0.0D0
        DO K=1,N
          Q = Q+ BETA(K)*(RS/RL)**(K-1)
        ENDDO
        Q = Q*(1.0D0/RL)**(L*(M+1))
        V = V + Q*ALPHA(I)
      ENDDO
      VRKHS = V

      END
C ----------------------------------------------------
      SUBROUTINE CALCPL(NANGLE,THETA,P)
C ----------------------------------------------------
      IMPLICIT NONE
      INTEGER I, NANGLE, MU
      DOUBLE PRECISION THETA,PI
      DOUBLE PRECISION P(NANGLE)
      DOUBLE PRECISION PN(NANGLE)

      P(1)=1.D0
      P(2)=THETA
      DO I=1,NANGLE-2
        P(I+2)=(P(I+1)*P(2)*DBLE(2*I+1)-P(I)*DBLE(I))/DBLE(I+1)
      ENDDO

      END
C ----------------------------------------------------
      SUBROUTINE GETRAD(R, RR, V, NANG, NRA, RPT, WT,
     O                  ALPHA, BETA, N, M, L, NRMAX)
C ----------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(NRMAX, NANG), NRA(NANG), ALPHA(NRMAX, NANG)
      DIMENSION RPT(NRMAX,NANG), BETA(N), RR(NANG), WT(NANG),P(NANG)
      LOGICAL FLAT

      DO I=1,NANG
        RR(I) = VRKHS(R,N,M,L,NRA(I),BETA,RPT(1,I),ALPHA(1,I))
      ENDDO
      END
