      SUBROUTINE ASSLEG(P,LMAX,X,N)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  GENERATE ASSOCIATED LEGENDRE FUNCTIONS
      DIMENSION P(N)

      P0=1.D0
      P1=X
      P(1)=P0
      P(2)=P1
      IND=2

      DO 100 L=2,LMAX
        TEMP=(DBLE(L+L-1)*X*P1-DBLE(L-1)*P0)/DBLE(L)
        P0=P1
        P1=TEMP
        IND=IND+L
        P(IND)=P1
  100 CONTINUE
C
C  NOW THE ONES FOR K.GT.0
C  NOTE THAT THIS ISN'T STABLE FOR HIGH ORDER AND X NEAR +/-1
C
      SINTH=SQRT(1.D0-X*X)

      DO 200 K=1,LMAX
        IND=K*(K+1)/2
      DO 200 L=K,LMAX
        INDM=IND
        IND=IND+L
        INDP=IND+1
        P(INDP)=-(DBLE(L-K+1)*X*P(IND)-DBLE(L+K-1)*P(INDM))/SINTH
  200 CONTINUE

      RETURN
      END
