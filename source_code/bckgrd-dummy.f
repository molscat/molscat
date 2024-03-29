      SUBROUTINE BCKGRD(B,SCLNRE,SCLNIM)
C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE efvs, ONLY: LEFVN, LEFVU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(10)  RUNAM
      CHARACTER(LEFVN) CONNAM
      CHARACTER(LEFVU) CONUNT
C ROUTINE FOR ADJUSTING FOR KNOWN BACKGROUND EFFECTS WHEN CONVERGING ON
C AND CHARACTERISING RESONANCES.
C ON ENTRY: B IS EFV OR ENERGY DEPENDING ON WHICH VARIABLE THE
C             CONVERGENCE IS A FUNCTION OF
C           SCLNRE IS EITHER THE REAL PART OF THE SCATTERING LENGTH, THE
C             EIGENPHASE SUM, OR THE REAL PART OF THE S-MATRIX ELEMENT,
C             DEPENDING ON WHAT IS BEING CONVERGED ON.
C           SCLNIM IS EITHER THE IMAGINARY PART OF THE SCATTERING LENGTH
C             OR THE IMAGINARY PART OF THE S-MATRIX ELEMENT, DEPENDING
C             ON WHAT IS BEING CONVERGED ON. FOR CONVERGENCE ON A
C             RESONANCE IN EIGENPHASE SUM IT IS SET TO 0.D0.
C ON EXIT THESE SHOULD HAVE BEEN CHANGED TO TAKE ACOUNT OF THE RELEVANT
C BACKGROUND EFFECTS
C NOTE THAT NOTHING ELSEWHERE IN THE PROGRAM KNOWS IF/WHAT BACKGROUND IS
C BEING ACCOUNTED FOR, SO ANY MESSAGES TO BE WRITTEN TO THE OUTPUT FILE
C MUST BE DONE IN THIS ROUTINE
      RETURN
      ENTRY BCKINI(IDECAY,IPRINT,RUNAME,CONNAM,CONUNT)
C INITIALISATION CALL CONTAINING INFORMATION ON WHAT TYPE OF CONVERGENCE
C AND SOME DETAILS THAT MAY BE NEEDED FOR PRINTING MESSAGES. IF ANYTHING
C NEEDS TO BE READ IN BY THIS ROUTINE, IT SHOULD BE DONE HERE
      RETURN
      END SUBROUTINE
