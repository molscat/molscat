      SUBROUTINE KSYM(AK,N)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS SUBROUTINE SYMMETRISES THE MATRIX AK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AK(N,N)

      DO 10 I=1,N
      DO 10 J=1,I
        TMP=0.5D0*(AK(I,J)+AK(J,I))
        AK(I,J)=TMP
        AK(J,I)=TMP
   10 CONTINUE

      RETURN
      END
