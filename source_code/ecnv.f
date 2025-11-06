      SUBROUTINE ECNV(EUNITS,TOCM,UNAME,IPRINT)
C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE physical_constants, ONLY: GHz_in_inv_cm, Hz_in_inv_cm,
     c                              K_in_inv_cm, MHz_in_inv_cm,
     c                              eV_in_inv_cm, erg_in_inv_cm,
     c                              hartree_in_inv_cm,
     c                              kJ_per_mol_in_inv_cm,
     c                              kcal_per_mol_in_inv_cm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  THIS ROUTINE ACCEPTS A 4 BYTE INPUT - EUNITS - AND DETERMINES
C  A UNITS TYPE AND ACCORDINGLY A CONVERSION FACTOR TO (CM-1).
C
C  IMPLEMENTED UNITS ARE
C  1) CM-1   2) DEG. K  3) MHZ  4) GHZ  5) EV  6) ERG  7) A.U.
C  8) KJ/MOL 9) KCAL/MOL
C
C  16-10-16: THIS ROUTINE NOW PRINTS OUT A STATEMENT ABOUT WHICH SET OF
C            PHYSICAL CONSTANTS ARE BEING USED.
C
C  ON ENTRY: EUNITS SHOULD BE AN INTEGER (1-9) SPECIFYING THE CORRESPONDING
C            UNIT OR A CHARACTER CODE CORRESPONDING TO UNIT.  UNCHANGED.
C            IF EUNITS IS NOT AN INTEGER, ECNVX IS CALLED AND ATTEMPTS TO
C            INTERPRET IT AS A STRING.
C  ON EXIT:  TOCM IS NUMERICAL VALUE CORRESPONDING TO CONVERSION FROM UNIT TO
C            WAVENUMBERS;
C            UNAME IS NAME OF UNIT.
C
C  NOTE THAT 1 ERG = 1E-7 JOULE
C  AND THE CALORIE HERE IS THE THERMOCHEMICAL CALORIE (4.184 J BY
C  DEFINITION), NOT THE "INTERNATIONAL TABLE" CALORIE (4.1868 J)
C
      INTEGER EUNITS
      PARAMETER (MXUNIT=9)
      CHARACTER(8) LTYP(MXUNIT),UNAME
      DIMENSION ECONV(MXUNIT)
      DATA LTYP /'CM-1', 'K', 'MHZ', 'GHZ', 'EV ', 'ERG', 'AU',
     1           'KJ/MOL', 'KCAL/MOL'/
C
      ECONV(1)=1d0
      ECONV(2)=K_in_inv_cm
      ECONV(3)=MHz_in_inv_cm
      ECONV(4)=GHz_in_inv_cm
      ECONV(5)=eV_in_inv_cm
      ECONV(6)=erg_in_inv_cm
      ECONV(7)=hartree_in_inv_cm
      ECONV(8)=kJ_per_mol_in_inv_cm
      ECONV(9)=kcal_per_mol_in_inv_cm

      IVAL=EUNITS
C  ARITHMETIC IF: <0,  =0,  >0
C     IF (IVAL) 2000,1000,1001
      IF (IVAL.LT.0 .OR. IVAL.GT.MXUNIT) THEN
C  IF EUNITS IS STRING RATHER THAN INTEGER VALUE, ECNVX IS CALLED
C  AND ATTEMPTS TO INTERPRET IT
C
C  IF ALPHANUMERICS CANNOT BE SUPPORTED, BELOW SHOULD PRINT ERROR
C  MESSAGE AND TERMINATE.
C
        CALL ECNVX(EUNITS,IVAL)
        TOCM=ECONV(IVAL)
        UNAME=LTYP(IVAL)
      ELSEIF (IVAL.EQ.0) THEN
        TOCM=1.D0
        UNAME=LTYP(1)
C
      ELSE
        TOCM=ECONV(IVAL)
        UNAME=LTYP(IVAL)
      ENDIF
      RETURN
      END
