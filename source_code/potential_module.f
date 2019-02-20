      module potential
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE sizes, ONLY: MXOMEG, MXLMDA
      IMPLICIT NONE

C  DATA FOR POTENTIAL-RELATED QUANTITIES
      INTEGER          :: NDGVL, NCONST, NRSQ, IREF, NVLBLK, NEXTRA,
     1                    NEXBLK, NEXTMS(MXOMEG), LAMBDA(MXLMDA)

      DOUBLE PRECISION :: VCONST(MXOMEG), RSCALE, EPSIL

      CHARACTER(10)    :: RUNAME, RMNAME, EPNAME
c  VL is arranged: MXLAM, NCONST, NRSQ, NEXBLK (total NVLBLK)
c                  |------------------|
c                    normally = NPOTL
C  REPLACES THIS COMMON BLOCK
C     COMMON /VLUSE / NCONST,NRSQ,VCONST,IREF,
C    1                IFIELD,NNCONST,MCONST
      end module potential
