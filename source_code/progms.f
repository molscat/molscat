C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
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
     4       'Copyright (C) 2025 J. M. Hutson & C. R. Le Sueur',14X,'|'/
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

      WRITE(6,10) TOTIME
   10 FORMAT(' |',22X,'This run used',F11.2,' cpu secs',21X,'|',/,
     9       ' |',76X,'|')
      RETURN
      END
C=========================================================================
      SUBROUTINE CPRMS(PDATE)
C  PRINT COPYRIGHT INFORMATION
      IMPLICIT NONE

      COMMON /CNTROL/ CDRIVE
      CHARACTER(1) CDRIVE

      CHARACTER(20) PDATE
      CHARACTER(48) CODE1
      CHARACTER(80) CODE2
      CHARACTER(80) PAPER

      IF (CDRIVE.EQ.'M') THEN
        CODE1='MOLSCAT: a program for non-reactive quantum'
        CODE2='scattering calculation on atomic and molecular '//
     1        'collisions'
        PAPER='J. M. Hutson & C. R. Le Sueur, Comput. Phys. '//
     1        'Commun., 241, 9-18 (2019).'
      ELSEIF (CDRIVE.EQ.'B') THEN
        CODE1='BOUND: a program for bound states of'
        CODE2='interacting pairs of atoms and molecules'
        PAPER='J. M. Hutson & C. R. Le Sueur, Comput. Phys. '//
     1        'Commun., 241, 1-8 (2019).'
      ELSEIF (CDRIVE.EQ.'F') THEN
        CODE1='FIELD: a program for bound states of'
        CODE2='interacting pairs of atoms and molecules as a '//
     1        'function of external field'
        PAPER='J. M. Hutson & C. R. Le Sueur, Comput. Phys. '//
     1        'Commun., 241, 1-8 (2019).'
      ENDIF

      WRITE(6,10) TRIM(CODE1),TRIM(CODE2),TRIM(PDATE),TRIM(PAPER)

   10 FORMAT(//2X,'This program is free software: you can ',
     1            'redistribute it and/or modify it under'/
     2         2X,'the terms of the GNU General Public License, ',
     3            'version 3, as published by'/
     4         2X,'the Free Software Foundation.'//
     5         2X,'Publications resulting from the use of this ',
     6            'program should cite both the'/
     7         2X,'version of the program used:'//
     8         2X,'J. M. Hutson & C. R. Le Sueur, ',A/
     9         2X,A,','/2X,'Version ',A,', ',
     A            'https//github.com/molscat/molscat'//
     B         2X,'and the published paper:'//
     C         2X,A/)
C MOLSCAT:       DOI REFERENCE IS: https://doi.org/10.1016/j.cpc.2019.02.014
C BOUND & FIELD: DOI REFERENCE IS: https://doi.org/10.1016/j.cpc.2019.02.017
      RETURN
      END
C=========================================================================
