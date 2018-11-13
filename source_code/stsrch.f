      LOGICAL FUNCTION STSRCH(LETTER,LSTR,N,I4)
C  Copyright (C) 2018 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  FUNCTION TO COMPARE ELEMENTS OF LSTR WITH LETTER
C  ON ENTRY: LETTER IS SINGLE UPPERCASE LETTER;
C            LSTR IS ARRAY OF LETTERS OF LENGTH (N).
C            ALL UNCHANGED ON EXIT.
C  ON EXIT:  I4 IS INDEX OF ELEMENT OF LSTR THAT MATCHES LETTER
C               OR 0 IF NO MATCH.
C
      CHARACTER(1) LETTER,LSTR(N)
      IF (N.LE.0) GOTO 9000
      DO 1000 I=1,N
        IF (LSTR(I).NE.LETTER) GOTO 1000
        I4=I
        STSRCH=.TRUE.
        RETURN
 1000 CONTINUE

 9000 STSRCH=.FALSE.
      I4=0
      RETURN
      END