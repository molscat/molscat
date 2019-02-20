      SUBROUTINE GDATE(DATE)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THESE ROUTINES ARE MACHINE-DEPENDENT, AND MUST BE SIMULATED.
C  THEY SHOULD RETURN STRINGS CONTAINING THE CURRENT DATE & TIME.
C  THIS VERSION IS FOR ANY COMPILER SUPPORTING THE f90 DATE_AND_TIME FUNCTION.
C
      CHARACTER DATE*11
      CHARACTER CDATE*8, CTIME*10, CZONE*5
      CHARACTER(3) MONTH(12)
      INTEGER VALUES(8)
      DATA MONTH/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     1           'Oct','Nov','Dec'/

      CALL date_and_time(CDATE,CTIME,CZONE,VALUES)
      DATE=CDATE(7:8)//' '//MONTH(VALUES(2))//' '//CDATE(1:4)

      RETURN
      END
C========================================================================
      SUBROUTINE GTIME(TIME)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      CHARACTER TIME*9
      CHARACTER CDATE*8, CTIME*10, CZONE*5
      INTEGER VALUES(8)

      CALL date_and_time(CDATE,CTIME,CZONE,VALUES)
      TIME=CTIME(1:2)//':'//CTIME(3:4)//':'//CTIME(5:6)

      RETURN
      END
