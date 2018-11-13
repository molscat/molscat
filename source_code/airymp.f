      SUBROUTINE  AIRYMP(X, FTHETA, FPHI, XMMOD, XNMOD)
C  This subroutine is part of the MOLSCAT, BOUND and FIELD suite of programs
C
C  AUTHOR:  MILLARD ALEXANDER
C  CURRENT REVISION DATE: 23-SEPT-87
C  SUBROUTINE TO RETURN THE MODULI AND PHASES OF THE AIRY FUNCTIONS AND
C  DERIVATIVES
C ----------------------------------------------------------------------
C  VARIABLES IN CALL LIST:
C    X     ARGUMENT OF AIRY FUNCTIONS
C    FTHETA, XMMOD         ON RETURN: CONTAIN THE (DOUBLE PRECISION)
C                                 PHASE AND MODULUS OF AI(X) AND BI(X)
C                                 (SEE BELOW).
C    FPHI, XNMOD           ON RETURN: CONTAIN THE (DOUBLE PRECISION)
C                                 PHASE AND MODULUS OF AI'(X) AND BI'(X)
C                                 (SEE BELOW).
C ----------------------------------------------------------------------
C  FOR NEGATIVE X
C ----------------------------------------------------------------------
C  THE MODULI AND PHASES ARE DEFINED BY
C      AI(-X) = M(X) COS[THETA(X)]
C      BI(-X) = M(X) SIN[THETA(X)]
C      AI'(-X) = N(X) COS[PHI(X)]
C      BI'(-X) = N(X) SIN[PHI(X)]
C  IN OTHER WORDS
C          2              2        2
C      M(X)  = SQRT[ AI(X)  + BI(X)  ]
C          2               2         2
C      N(X)  = SQRT[ AI'(X)  + BI'(X)  ]
C      THETA(X) = ATAN [ BI(X) / AI(X) ]
C      PHI(X)   = ATAN [ BI'(X) / AI'(X) ]
C  TO DETERMINE THESE MODULI AND PHASES WE USE THE SUBROUTINE
C  SCAIRY, WRITTEN BY D. MANOLOPOULOS (SEPT. 1986)
C  THIS SUBROUTINE RETURNS THE FOLLOWING QUANTITIES:
C     SCAI, SCBI, SCAIP, SCPIB, AND ZETA, WHERE
C     FOR  X .LT. -5.0
C     AI(X) = SCAI * COS(ZETA) + SCBI * SIN(ZETA)
C     BI(X) = SCBI * COS(ZETA) - SCAI * SIN(ZETA)
C     AI'(X) = SCAIP * COS(ZETA) + SCBIP * SIN(ZETA)
C     BI'(X) = SCBIP * COS(ZETA) - SCAIP * SIN(ZETA)
C     WHERE ZETA = (2/3) * (-X) ** (3/2) + PI/4
C
C     FOR  -5.0 .LE. X .LE. 0.0
C
C     AI(X) = SCAI
C     BI(X) = SCBI
C     AI'(X) = SCAIP
C     BI'(X) = SCBIP
C     AND ZETA = 0
C ----------------------------------------------------------------------
C  FOR POSITIVE X
C ----------------------------------------------------------------------
C  THE MODULI AND PHASES ARE DEFINED BY
C      AI(X) = M(X) SINH[THETA(X)]
C      BI(X) = M(X) COSH[THETA(X)]
C      AI'(X) = N(X) SINH[PHI(X)]
C      BI'(X) = N(X) COSH[PHI(X)]
C  IN OTHER WORDS
C          2              2        2
C      M(X)  = SQRT[ BI(X)  - AI(X)  ]
C          2               2         2
C      N(X)  = SQRT[ BI'(X)  - AI'(X)  ]
C      THETA(X) = ATANH [ AI(X) / BI(X) ]
C      PHI(X)   = ATANH [ AI'(X) / BI'(X) ]
C  HERE THE THE EXPONENTIALLY SCALED AIRY FUNCTIONS
C  AI(X), AI'(X), BI(X), BI'(X) ARE:
C      AI(X)  = AI(X)  * EXP[ZETA]
C      AI'(X) = AI'(X) * EXP[ZETA]
C      BI(X)  = BI(X)  * EXP[-ZETA]
C      BI'(X) = BI'(X) * EXP[-ZETA]
C  TO DETERMINE THESE MODULI AND PHASES WE USE THE SUBROUTINE
C  SCAIRY, WRITTEN BY D. MANOLOPOULOS (SEPT. 1986)
C  THIS SUBROUTINE RETURNS THE FOLLOWING QUANTITIES:
C     SCAI, SCBI, SCAIP, SCPIB, AND ZETA
C  IN TERMS OF WHICH THE EXPONENTIALLY SCALED AIRY FUNCTIONS ARE DEFINED
C   AI(X) = SCAI * EXP(-ZETA)
C   BI(X) = SCBI * EXP(+ZETA)
C   AI'(X) = SCAIP * EXP(-ZETA)
C   BI'(X) = SCBIP * EXP(+ZETA)
C   WHERE ZETA = (2/3) * X ** (3/2)
C
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X, FTHETA, FPHI, XMMOD, XNMOD, SCAI,
     :                 SCBI, SCAIP, SCBIP, ZETA, RATIO
      CALL SCAIRY(X, SCAI, SCBI, SCAIP, SCBIP, ZETA)
      IF (X.LE.0.D0) THEN
        XMMOD = SQRT( SCAI ** 2 + SCBI ** 2)
        XNMOD = SQRT( SCAIP ** 2 + SCBIP ** 2)
        FTHETA = ATAN2(SCBI, SCAI)
        FPHI = ATAN2(SCBIP, SCAIP)
        IF (X.LT.(-5.0D0) ) THEN
          FTHETA = FTHETA - ZETA
          FPHI = FPHI - ZETA
        ENDIF
      ELSE
        XMMOD = SQRT( - SCAI ** 2 + SCBI ** 2)
        XNMOD = SQRT( - SCAIP ** 2 + SCBIP ** 2)
        RATIO = SCAI / SCBI
        FTHETA = 0.5D0 * LOG( (1.D0 + RATIO) / (1.D0 - RATIO) )
        RATIO = SCAIP / SCBIP
        FPHI = 0.5D0 * LOG( (1.D0 + RATIO) / (1.D0 - RATIO) )
      ENDIF
      RETURN
      END
