      SUBROUTINE VINIT(NV,RUNIT,VUNIT)
C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS VERSION OF VINIT IMPLEMENTS THE POTENTIAL FORM OF
C  TIEMANN AND CO-WORKERS FOR DIATOMIC POTENTIAL CURVES.
C  IT MAKES USE OF A MODULE pot_data_Tiemann THAT MUST CONTAIN THE
C  POTENTIAL PARAMETERS FOR THE SPECIFIC SYSTEM.
C  THIS VERSION IS DESIGNED FOR SYSTEMS LIKE ALKALI-METAL DIATOMICS,
C  WITH TWO CURVES THAT DESCRIBE SINGLET AND TRIPLET STATES.
      USE physical_constants, only: bohr_to_Angstrom
      USE potential, ONLY: RMNAME, EPNAME
      USE pot_data_Tiemann
      IMPLICIT NONE
      LOGICAL :: LSETUP, G2B
      INTEGER :: NV, I
      DOUBLE PRECISION :: RUNIT, VUNIT
      DOUBLE PRECISION :: R, RA, XI, V, VMID, VLR, VEXCH
      DOUBLE PRECISION :: DVLR, DVDR, DXIDR, DXINDR
      DOUBLE PRECISION :: POWER, DPOWER
      DOUBLE PRECISION :: GAMOLD, BETOLD, ANEW, ASRNEW, BSRNEW
      DOUBLE PRECISION :: EPSIL = 1.D0
      DOUBLE PRECISION :: VMAX = 1.D4
      SAVE

  100 FORMAT(1X,A,1PG20.13:,A,G20.13,A,G20.13)
  200 FORMAT(/1X,A,I2,A,1PG20.13,A)
C
C  NV=1:  1SIGMA POTENTIAL
C  NV=2:  3SIGMA POTENTIAL
C
      IF (NV.LT.1 .OR. NV.GT.2) THEN
          WRITE(6,*) "NV OUT OF RANGE -- NV:",NV
          STOP ' PROGRAM HALTED IN POTENTIAL SYS_SS_POT '
      ENDIF
      IF (IPRINT.GE.1) WRITE(6,*) ' JMH routine for Tiemann-style',
     1                            ' alkali dimer potentials'
      IF (IPRINT.GE.1) WRITE(6,*) TRIM(POTNAM)
C
C  THREE CHOICES HERE:
C  CALCULATE EXCHANGE POWER GAMMA FROM EXPONENT BETA,
C  OR BETA FROM GAMMA, OR LEAVE BOTH UNCHANGED
C
      GAMOLD=GAMMA
      BETOLD=BETA
      IF (GAMBET.EQ.2) THEN
        BETA = 7.D0 / (bohr_to_Angstrom * (GAMMA + 1.D0))
      ELSEIF (GAMBET.EQ.1) THEN
        GAMMA = 7.D0 / (BETA * bohr_to_Angstrom) - 1.D0
      ENDIF

      IF (IPRINT.GE.2) THEN
        WRITE(6,*)
        IF (GAMBET.EQ.2) THEN
          WRITE(6,100) ' beta  shifted from ',BETOLD,' A-1'
          WRITE(6,100) '                 to ',BETA,' A-1'
        ELSEIF (GAMBET.EQ.1) THEN
          WRITE(6,100) ' gamma shifted from ',GAMOLD
          WRITE(6,100) '                 to ',GAMMA
        ELSE
          WRITE(6,*) ' input gamma and beta unchanged'
          WRITE(6,100) ' beta  = ',BETA,' A-1'
          WRITE(6,100) ' gamma = ',GAMMA
        ENDIF
      ENDIF
C
C  MATCH POTENTIAL AT LONG-RANGE POINT
C  THE PREFERRED PROCEDURE IS TO CALCULATE THE VALUE OF A(0)
C  TO MATCH THE VALUE OF THE LONG-RANGE POTENTIAL AT RLR.
C  HOWEVER, TO REPRODUCE RESULTS FROM OTHER ROUTINES IT IS SOMETIMES
C  DESIRED TO FIX A(0) AT THE ROUNDED VALUE IN A PUBLISHED PAPER.
C  THIS IS ACHIEVED BY SETTING MATCHL TO .FALSE. IN THE DATA MODULE.
C  IF THIS IS DONE, THE ROUTINE PRINTS OUT THE RESULTING MISMATCH.
C  THERE IS ALWAYS A DEIVATIVE DISCONTINUITY AT RLR, WHICH IS PRINTED.
C
      VEXCH = AEX * RLR(NV)**GAMMA * EXP(-BETA*RLR(NV))
      VLR = EXSIGN(NV)*VEXCH
     X    - C6/RLR(NV)**6 - C8/RLR(NV)**8 - C10/RLR(NV)**10
      IF (NEX.GT.0) VLR = VLR - CEX/RLR(NV)**NEX
      DVLR = EXSIGN(NV)*VEXCH*(1.D0/RLR(NV)-BETA)
     X    + 6.D0*C6/RLR(NV)**7 + 8.D0*C8/RLR(NV)**9
     X    + 10.D0*C10/RLR(NV)**11
      IF (NEX.GT.0) DVLR = DVLR + DBLE(NEX)*CEX/RLR(NV)**(NEX+1)
      XI = (RLR(NV) - RM(NV))/(RLR(NV)+B(NV)*RM(NV))
      VMID = POWER(XI,A(0,NV),NA(NV))
C  DERIVATIVE OF XI
      DXIDR = (B(NV)+1.D0) * RM(NV) / (RLR(NV)+B(NV)*RM(NV))**2
      DVDR = DXIDR*DPOWER(XI,A(0,NV),NA(NV))
      ANEW = A(0,NV)+VLR-VMID

      IF (IPRINT.GE.2) THEN
        WRITE(6,200) ' For potential ',NV,' at RLR =',RLR(NV),' A'
        WRITE(6,100) '                 with C6 =',C6,' and C8 =',C8
        WRITE(6,100) ' V_LR  = ',VLR,'cm-1'
        WRITE(6,100) ' V_mid = ',VMID,'cm-1, dV/dR = ',DVDR,'cm-1/A'
        IF (MATCHL) THEN
          WRITE(6,100) ' A(0) shifted from ',A(0,NV),' cm-1'
          WRITE(6,100) '                to ',ANEW ,' cm-1'//
     1                 ' to match value of V(R) at RLR'
        ELSE
          WRITE(6,*) ' A(0) not shifted to match V(R) at RLR'
          WRITE(6,100) ' Potential discontinuity at RLR is ',
     1                 ANEW-A(0,NV),' cm-1'
          WRITE(6,100) ' Discontinuity / value is ',(ANEW-A(0,NV))/VMID
        ENDIF

        WRITE(6,100) ' Derivative discontinuity at RLR'//
     1               ' (intrinsic to functional form) is ',DVLR-DVDR,
     1               ' cm-1/A'
        WRITE(6,100) ' Discontinuity / value is ',(DVLR-DVDR)/DVLR
      ENDIF

      IF (MATCHL) A(0,NV) = ANEW
C
C  MATCH AT SHORT-RANGE POINT
C  THE PREFERRED PROCEDURE IS TO CALCULATE VALUES OF ASR AND BSR
C  FROM THE POTENTIAL AND ITS DERIVATIVE AT RSR.
C  HOWEVER, TO REPRODUCE RESULTS FROM OTHER ROUTINES IT IS SOMETIMES
C  DESIRED TO FIX ASR AND/OR BSR AT THE ROUNDED VALUES IN A PUBLISHED
C  PAPER. THIS IS ACHIEVED BY SETTING MATCHV AND/OR MATCHD TO .FALSE.
C  IN THE DATA MODULE. IF THIS IS DONE, THE ROUTINE PRINTS OUT THE
C  RESULTING MISMATCH FOR INFORMATION.
C
      XI = (RSR(NV) - RM(NV))/(RSR(NV)+B(NV)*RM(NV))
      VMID = POWER(XI,A(0,NV),NA(NV))
C  DERIVATIVE OF XI
      DXIDR = (B(NV)+1.D0) * RM(NV) / (RSR(NV)+B(NV)*RM(NV))**2
      DVDR = DXIDR*DPOWER(XI,A(0,NV),NA(NV))

      BSRNEW = -DVDR*RSR(NV)**(NSR(NV)+1.D0)/NSR(NV)
      IF (IPRINT.GE.2) THEN
        WRITE(6,200) ' For potential ',NV,' at RSR =',RSR(NV),' A'
        WRITE(6,100) ' V     = ',VMID,'cm-1, dV/dR = ',DVDR,
     1             'cm-1/A and n = ',NSR(NV)
        IF (MATCHD) THEN
          WRITE(6,100) ' B(SR) shifted from ',BSR(NV),' cm-1 A^n'
          WRITE(6,100) '                 to ',BSRNEW ,' cm-1 A^n'
        ELSE
          WRITE(6,*) ' B(SR) not shifted to match dV/dR'
          WRITE(6,100) ' Derivative of V_SR at RSR is ',
     1                 -NSR(NV)*BSR(NV)/RSR(NV)**(NSR(NV)+1.d0),
     2                 ' cm-1/A'
          WRITE(6,100) ' Derivative discontinuity  is ',
     1                 -NSR(NV)*(BSRNEW-BSR(NV))/RSR(NV)
     2                 **(NSR(NV)+1.d0),' cm-1/A'
          WRITE(6,100) ' Discontinuity / value     is ',
     1                 (BSRNEW-BSR(NV))/BSRNEW
        ENDIF
      ENDIF
      IF (MATCHD) BSR(NV) = BSRNEW

      ASRNEW = VMID - BSR(NV)/RSR(NV)**NSR(NV)
      IF (IPRINT.GE.2) THEN
        IF (MATCHV) THEN
          WRITE(6,100) ' A(SR) shifted from ',ASR(NV),' cm-1'
          WRITE(6,100) '                 to ',ASRNEW ,' cm-1'
        ELSE
          WRITE(6,*) ' A(SR) not shifted to match V(R)'
          WRITE(6,100) ' Potential discontinuity at RSR is ',
     1                 ASRNEW-ASR(NV),' cm-1'
          WRITE(6,100) ' Discontinuity / value is ',
     1                 (ASRNEW-ASR(NV))/ASRNEW
        ENDIF
      ENDIF
      IF (MATCHV) ASR(NV) = ASRNEW
C
C  SET MOLSCAT/BOUND ENERGY UNITS TO CM-1
C  SET MOLSCAT/BOUND LENGTH UNITS TO A OR BOHR AS IN DATA MODULE
C  NOTE THAT THIS CHOICE OVERRIDES ANY VALUE OF RM INPUT IN &POTL
C  AND DETERMINES THE UNITS THAT MUST BE USED FOR ANY SPIN-SPIN TERM
C
C  INTERNAL UNITS USED IN THIS VSTAR ARE ALWAYS CM-1 AND A
C  R IS CONVERTED TO R ON ENTRY AND V IS CONVERTED TO EPSIL ON EXIT
C
      RUNIT = RUNITM
      VUNIT = EPSIL
      RMNAME(1:8)=LENUNT
      RETURN
C
      ENTRY VSTAR(NV,R,V)
C
C  CALCULATE POTENTIAL POINT:
C  FIRST CONVERT INPUT R FROM MOLSCAT RM UNITS TO ANGSTROM
C
      RA = R*RUNITM
C
      IF (RA.LT.RSR(NV)) THEN
          V = ASR(NV)+BSR(NV)/RA**NSR(NV)
      ELSEIF (RA.LE.RLR(NV)) THEN
          XI = (RA - RM(NV))/(RA+B(NV)*RM(NV))
          V = POWER(XI,A(0,NV),NA(NV))
      ELSE
          VEXCH = AEX * RA**GAMMA * EXP(-BETA*RA)
          V = EXSIGN(NV)*VEXCH - C6/RA**6 - C8/RA**8 - C10/RA**10
          IF (NEX.GT.0) V = V - CEX/RA**NEX
      ENDIF
C
      V = V / EPSIL
C
C     LIMIT VERY LARGE POSITIVE VALUES TO IMPROVE STEP-SIZE CONVERGENCE
C     PARTICULARLY FOR IPROP = 6
C
      IF (V.GT.VMAX) V = VMAX
      RETURN
C
      ENTRY VSTAR1(NV,R,V)
C
C  CALCULATE DERIVATIVE POINT:
C  FIRST CONVERT INPUT R FROM MOLSCAT RM UNITS TO ANGSTROM
C
      RA = R*RUNITM
C
      IF (RA.LT.RSR(NV)) THEN
        V = -NSR(NV)*BSR(NV)/RA**(NSR(NV)+1.D0)
      ELSEIF (RA.LE.RLR(NV)) THEN
        XI = (RA - RM(NV))/(RA+B(NV)*RM(NV))
        DXIDR = (B(NV)+1.D0) * RM(NV) / (RA+B(NV)*RM(NV))**2
        DVDR = DXIDR*DPOWER(XI,A(0,NV),NA(NV))
        V = DVDR
      ELSE
          VEXCH = AEX * RA**GAMMA * EXP(-BETA*RA)
          V = EXSIGN(NV) * (GAMMA/RA-BETA) * VEXCH
     1      + 6.D0*C6/RA**7 + 8.D0*C8/RA**9 + 10.D0*C10/RA**11
          IF (NEX.GT.0) V = V + DBLE(NEX)*CEX/RA**(NEX+1)
      ENDIF
C
C     CONVERT DERIVATIVE TO EXTERNAL LENGTH UNITS
C
      V = V * RUNITM / EPSIL
C
      RETURN
C
      ENTRY VSTAR2(NV,R,V)
C
C  SECOND DERIVATIVES NOT IMPLEMENTED BUT WOULD BE EASY IF NEEDED
C
      WRITE(6,*) ' CALLED VSTAR2: SECOND DERIVATIVES OF',
     1           ' TIEMANN-STYLE POTENTIAL NOT IMPLEMENTED'
      STOP
      END
C
      DOUBLE PRECISION FUNCTION POWER(X,A,N)
C
C  EVALUATE A POWER SERIES WITH COEFFICIENTS IN A
C
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION X,A(N)

      POWER=0.D0
      DO 10 I=N,2,-1
        POWER=POWER+A(I)
        POWER=POWER*X
   10 CONTINUE
      POWER=POWER+A(1)

      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DPOWER(X,A,N)
C
C  EVALUATE DERIVATIVE OF A POWER SERIES WITH COEFFICIENTS IN A
C
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION X,A(N)

      DPOWER=0.D0
      DO 10 I=N,3,-1
        DPOWER=DPOWER+DBLE(I-1)*A(I)
        DPOWER=DPOWER*X
   10 CONTINUE
      DPOWER=DPOWER+A(2)

      RETURN
      END
