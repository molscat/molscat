      module potential
C  Copyright (C) 2025 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE sizes, ONLY: MXOMEG, MXLMDA
      IMPLICIT NONE

C  DATA FOR POTENTIAL-RELATED QUANTITIES
      INTEGER          :: NVLBLK, NCONST, NRSQ, NEXBLK, NEXTRA,
     1                    NEXTMS(MXOMEG), NDGVL, IREF, LAMBDA(MXLMDA)
c  VL is arranged: MXLAM, NCONST, NRSQ, NEXBLK (total NVLBLK)
c                  |------------------|
c                    normally = NHAM

      DOUBLE PRECISION :: CONLEN, CONFRQ, EP2RU, VCONST(MXOMEG)

      CHARACTER(10)    :: RMNAME, EPNAME

      LOGICAL          :: VRESET
C  REPLACES THIS COMMON BLOCK
C     COMMON /VLUSE / NCONST,NRSQ,VCONST,IREF,
C    1                IFIELD,NNCONST,MCONST
      end module potential
