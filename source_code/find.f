      INTEGER FUNCTION FIND(I,J,IG,NG)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  FUNCTION TO FIND A PARTICULAR FOURIER COMPONENT IN A LIST
C  OF COMPONENTS, AND RETURN THE POSITION OF THE REQUIRED
C  COMPONENT.
C
      DIMENSION IG(2,NG)
C
      II=I
      JJ=J
      CALL ORDER(II,JJ)
      FIND=0
      DO N=1,NG
        IF (II.NE.IG(1,N) .OR. JJ.NE.IG(2,N)) CYCLE
        FIND=N
        EXIT
      ENDDO
      RETURN
      END
