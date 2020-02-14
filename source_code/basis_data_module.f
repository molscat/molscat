      module basis_data
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE sizes, ONLY: MXELVL, MXJLVL, MXROTS, MXSYMS
      IMPLICIT NONE

C  DATA FOR BASIS SET-RELATED QUANTITIES
      INTEGER          :: IDENT, JHALF, ISYM(MXSYMS), ISYM2(MXSYMS),
     1                    JMIN, J2MIN, JMAX, J2MAX, JSTEP, J2STEP,
     2                    JLEVEL(MXJLVL), NJLQN9, NLEVEL

      DOUBLE PRECISION :: ELEVEL(MXELVL), EMAX, ROTI(MXROTS), SPNUC,
     1                    WT(2)
C  REPLACES THIS COMMON BLOCK
C     COMMON /CMBASE/ ROTI,ELEVEL,EMAX,
C    1                WT,SPNUC,NLEVEL,JLEVEL,JMIN,JMAX,JSTEP,
C    2                ISYM,J2MIN,J2MAX,J2STEP,ISYM2,JHALF,
C    3                IDENT
C  SPECIFICATIONS FOR MOLSCAT(&BASIS) COMPATIBILITY. . .
      DOUBLE PRECISION ALPHAE(2),BE(2),DE(2),A(2),B(2),C(2),WE(2),
     1                 WEXE(2),DJ,DJK,DK,DT
      INTEGER J1MIN,J1MAX,J1STEP,KMAX,KSET
      EQUIVALENCE (ROTI(1),BE(1),A(1)),
     1            (ROTI(3),ALPHAE(1),B(1)),
     2            (ROTI(5),DE(1),C(1)),
     3            (ROTI(7),WE(1),DJ),
     4            (ROTI(8),DJK),
     5            (ROTI(9),WEXE(1),DK),
     6            (ROTI(10),DT),
     7            (JMIN,J1MIN),(JMAX,J1MAX),(JSTEP,J1STEP),
     8            (J2MAX,KMAX,KSET)
      end module basis_data
