C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3

C  This set of subroutines used to print out common information in header
C  and footer messages
C========================================================================
      SUBROUTINE PROGVS(PDATE)
C  PRINT VERSION NUMBER AND COPYRIGHT STATEMENT
      IMPLICIT NONE

      CHARACTER(20) PDATE

      WRITE(6,10) PDATE

   10 FORMAT(' |',76X,'|'/' |',14X,
     4       'Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur',14X,'|'/
     4       ' |',76X,'|'/
     5       ' |',31X,'Version ',A20,17X,'|'/
     6       ' |',76X,'|')
      RETURN
      END
C=========================================================================
      SUBROUTINE TIMEST(CDATE,CTIME)
C  PRINT TIME STAMP FOR BEGINNING OF CALCULATION
      IMPLICIT NONE

      CHARACTER(11) CDATE
      CHARACTER(9)  CTIME

      WRITE(6,10) CDATE,CTIME
   10 FORMAT(' |',23X,'Run on ',A11,2X,' at ',A9,20X,'|'/' |',76X,'|')
      RETURN
      END
C=========================================================================
      SUBROUTINE TIMEMS(TOTIME,NUSED,MXSAVE)
C  PRINT AMOUNT OF CPU TIME AND MEMORY USED
      IMPLICIT NONE

      INTEGER MXSAVE,NUSED
      DOUBLE PRECISION TOTIME

      WRITE(6,10) TOTIME,NUSED,MXSAVE
   10 FORMAT(' |',20X,'This run used',F11.2,' cpu secs and',19X,'|',/,
     8       ' |',9X,I10,' of the allocated',I10,' words of storage',
     9       13X,'|'/' |',76X,'|')
      RETURN
      END
C=========================================================================
      SUBROUTINE CPRMSG(CPROG,PDATE)
C  PRINT COPYRIGHT INFORMATION
      IMPLICIT NONE

      CHARACTER(7)  CPROG
      CHARACTER(20) PDATE

      WRITE(6,10) TRIM(CPROG),PDATE

   10 FORMAT(//2X,'This program is free software: you can ',
     1            'redistribute it and/or modify it under'/
     2         2X,'the terms of the GNU General Public License, ',
     3            'version 3, as published by'/
     4         2X,'the Free Software Foundation.'//
     5         2X,'Publications resulting from the use of this ',
     6            'program should refer to'/
     7         2X,'J. M. Hutson & C. R. Le Sueur, ',A,
     8            ' computer code ',A/
     9         2X,'or to the published paper in Computer Physics ',
     1            'Communications when available.'/)
      RETURN
      END
