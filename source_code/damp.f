      SUBROUTINE DAMP(KMIN,KMAX,BETA,R,P,DP)
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  SUBROUTINE FOR THE EFFICIENT CALCULATION OF THE TANG-TOENNIES
C  DAMPING FUNCTIONS AND THEIR FIRST DERIVATIVES
C  SEE J. CHEM. PHYS. 80, 3726 (1984).
C
      DIMENSION P(KMAX),DP(KMAX)
C
      Y=1.D0
      Z=1.D0
      BR=BETA*R
      DO K=1,KMAX
        Y=Y*BR/DBLE(K)
        Z=Z+Y
        P(K)=Z
      ENDDO
C
      Z=EXP(-BR)
      DO K=KMIN,KMAX
        DP(K)=(P(K)-P(K-1))*BETA*Z
      ENDDO
      DO K=KMIN,KMAX
        P(K)=1.D0-Z*P(K)
      ENDDO
      RETURN
      END
