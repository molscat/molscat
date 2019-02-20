      DOUBLE PRECISION FUNCTION EXTPOT(RA,CTH,PHI)
C  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE physical_constants
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  Ar-CH4 systematic potential, version 2.
C  Reads in the following parameters in ATOMIC UNITS
C  C8 RATIO ALPHA SCALE
C  C7 C9 RAT8
C  where RATIO=C10*C6/C8**2, ALPHA is the repulsion softness parameter,
C  SCALE is the repulsion scaling, RAT8=delta(C8)/C8,
C  delta(C8) is the anisotropic part of C8.
C  Defaults:
C  ALPHA=1.956/bohr (overlap model) will be used if input ALPHA.LE.0D0
C  C7=405 (20% more than Fowler) if input C7.LT.0D0
C  C9=0
C  RATIO=1.38 (from H-H exact values) if input RATIO.LT.0D0
C  C8=3600 (40% more than Fowler) if input C8.LT.0D0
C  RAT8=-0.17185 (the Fowler result) if input RAT8.GT.0D0
C  SCALE has no obvious default.
C
      DATA A6,A8,A10/0.3648D0,0.3073D0,0.2514D0/
      DATA B6,B8,B10/0.03360D0,0.02469D0,0.02379D0/
      DATA D6,D8,D10/0.001651D0,0.001227D0,0.0005664D0/
      DATA SCP,GMU,GEX/4.3751D0,41.34D0,0.8588D0/

      DAMP6(X) =(1D0-EXP(-X*(A6 +X*(B6 +D6*X))))**6
      DAMP8(X) =(1D0-EXP(-X*(A8 +X*(B8 +D8*X))))**8
      DAMP10(X)=(1D0-EXP(-X*(A10+X*(B10+D10*X))))**10
      CORR(X)=1D0+GMU*EXP(-GEX*X)
      C2PHI=2D0*COS(PHI)**2-1D0
      C4PHI=2D0*C2PHI*C2PHI-1D0
      CTH2=CTH*CTH
      STH2=1D0-CTH2
      F3=CTH*STH2*C2PHI
      F40=(3D0-CTH2*(30D0-35D0*CTH2))/8D0
      F44=STH2*STH2*C4PHI*SQRT(35D0/64D0)
      F4=SQRT(7D0)*F40-SQRT(5D0)*F44
C  Assumes FACE is -root(1/3),0; VERTEX is -root(1/3),90; EDGE is 1,0
      C6A=C6
      C7A=C7*F3/2D0
      C8A=C8
      C8A=C8A+RAT8*C8*F4*0.4D0/SQRT(7D0)
      C9A=C9*F3/2D0
      C10A=RATIO*C8A*C8A/C6A
      SDAMP=SCP*SQRT(C6A/C8A)
      R=RA/RM
      X=RM/RA
      Y=X*X
      REP=SCALE*EXP(-ALPHA*R)*(1D0+0.67D0*F3*SQRT(13.125D0/PI)
     1                           -0.069D0*F4*SQRT(4.5D0/PI))
      XD=SDAMP*R
      DISP=-Y*Y*Y*(C6A*DAMP6(XD)+X*(C7A*DAMP8(XD)+X*(C8A*DAMP8(XD)
     1     +X*(C9A*DAMP10(XD)+X*C10A*DAMP10(XD)))))
      EXTPOT=REP+DISP*CORR(XD)
      EXTPOT=EXTPOT*EPSIL
      RETURN
C================================================================= EXTINT
      ENTRY EXTINT
      C6=90.94D0
      READ(5,*) C8,RATIO,ALPHA,SCALE
      IF (C8.LT.0D0) C8=3600D0
      IF (RATIO.LT.0D0) RATIO=1.38D0
      IF (ALPHA.LE.0D0) ALPHA=1.956D0
      READ(5,*) C7,C9,RAT8
      IF (C7.LT.0D0) C7=405D0
      IF (RAT8.GT.0D0) RAT8=-0.17185D0
C   Conversion factors
      PI=ACOS(-1.D0)
C     RM=0.529177249D0
C     EPSIL=0.21947463067D6
C
C  20-09-2016: UPDATED TO USE MODULE (physical_constants) THAT CONTAINS
C              CONSISTENT AND UP-TO-DATE VALUES
      RM=bohr_to_Angstrom
      EPSIL=hartree_in_inv_cm
C
      WRITE(6,*) ' C6, C8, RATIO = ',C6,C8,RATIO
      WRITE(6,*) ' C7, RAT8, C9 = ',C7,RAT8,C9
      WRITE(6,*) ' ALPHA, SCALE = ',ALPHA,SCALE
      RETURN
      END
