      SUBROUTINE J6J(J2,J3,L1,L2,L3,IVAL,J1MIN,D6J)
C  This subroutine is part of the MOLSCAT, BOUND and FIELD suite of programs
      IMPLICIT DOUBLE PRECISION (A-H,J-Z)
      DIMENSION D6J(*)
      DATA ZERO/0.D0/,TENTH/0.1D0/,HALF/0.5D0/,ONE/1.D0/,TWO/2.D0/,
     $     CONST/1.D-12/
      E(J1S)=SQRT((J1S-MJ23S)*(J23P1S-J1S)*(J1S-ML23S)*(L23P1S-J1S))
      F(J1,JJP1)=(TWO*J1+ONE)*(JJP1*(FACT-JJP1-TWO*LLP1)+FACT2)
C
C  THIS ROUTINE CALCULATES THE 6-J {J1 J2 J3} COEFFICIENTS FOR ALL PERMISSIBLE
C                                  {L1 L2 L3}
C  VALUES OF J1 FOR FIXED VALUES OF J2, J3, L1, L2, AND L3 USING THE
C  RECURSIVE ALGORITHM OF K. SCHULTEN AND R. G. GORDON, J. MATH. PHYS.
C  VOL. 16, P. 1961, (1975).
C  PROGRAMMED BY D. E. FITZ, 10/22/79
C  MODIFIED BY S. GREEN 20 AUG 93 TO TEST DIMENSION ON D6J
C
      MXDIM=IVAL
      JJP2=J2*(J2+ONE)
      JJP3=J3*(J3+ONE)
      LLP1=L1*(L1+ONE)
      LLP2=L2*(L2+ONE)
      LLP3=L3*(L3+ONE)
      MJ23S=(J2-J3)**2
      ML23S=(L2-L3)**2
      J23P1S=(J2+J3+ONE)**2
      L23P1S=(L2+L3+ONE)**2
      FACT2=(LLP2-LLP3)*(JJP2-JJP3)
      FACT=JJP2+JJP3+LLP2+LLP3
      J1MIN=MAX(ABS(J2-J3),ABS(L2-L3))
      J1MAX=MIN(J2+J3,L2+L3)
      IVAL=INT(J1MAX-J1MIN+ONE+TENTH)
      IF (IVAL.GT.MXDIM) THEN
        WRITE(6,*) ' J6J: ARRAY D6J TOO SMALL. NEEDS ',IVAL,
     1             ' BUT ONLY ',MXDIM,' SUPPLIED'
        STOP
      ENDIF
C
C  TEST FOR OTHER TRIANGULAR INEQUALITES.
C
      IL1=INT(TWO*L1+TENTH)
      IL2=INT(TWO*L2+TENTH)
      IL3=INT(TWO*L3+TENTH)
      IJ2=INT(TWO*J2+TENTH)
      IJ3=INT(TWO*J3+TENTH)
      IF (IJ2.LE.IL1+IL3 .AND. IJ2.GE.ABS(IL1-IL3) .AND.
     $    IJ3.LE.IL1+IL2 .AND. IJ3.GE.ABS(IL1-IL2)) GOTO 11

      DO 12 I=1,IVAL
   12   D6J(I)=ZERO
      RETURN
C
   11 INMID=(IVAL+3)/2
      SGNV=J2+J3+L2+L3
      SGN=ONE
      ISIGN=INT(SGNV+TENTH)
      IF (MOD(ISIGN,2).NE.0) SGN=-ONE
      D6J(1)=HALF
C
C  UPWARD RECURSION.
C
      IF (IVAL.EQ.1) GOTO 40

      JJP1=J1MIN*(J1MIN+ONE)
      F1=F(J1MIN,JJP1)
      J1=J1MIN+ONE
      J1S=J1*J1
      E2=E(J1S)
      IF (J1MIN.LT.TENTH) GOTO 15

      D6J(2)=-F1*D6J(1)/(E2*J1MIN)
      GOTO 16

   15 D6J(2)=-HALF*(LLP2+JJP2-LLP1)*D6J(1)/SQRT(JJP2*LLP2)
   16 SCALED=D6J(2)
      IF (IVAL.EQ.2) GOTO 40

      DO 21 IJ2=3,INMID
        JJP1=J1*(J1+ONE)
        F1=F(J1,JJP1)
        J1=J1+ONE
        E1=E2
        J1S=J1*J1
        E2=E(J1S)
   21   D6J(IJ2)=-(F1*D6J(IJ2-1)+J1*E1*D6J(IJ2-2))/(E2*(J1-ONE))
      SCALED=D6J(INMID)
      IEXC=5
      IF (ABS(SCALED).GT.CONST) GOTO 18

      INMID=INMID-1
      SCALED=D6J(INMID)
      IEXC=3
      GOTO 30

   18 IF (IVAL.EQ.3) GOTO 40
C
C  DOWNWARD RECURSION.
C
   30 D6J(IVAL)=HALF
      J1=J1MAX
      J1S=J1*J1
      JJP1=J1*(J1+ONE)
      F1=F(J1,JJP1)
      E1=E(J1S)
      D6J(IVAL-1)=-F1*D6J(IVAL)/(E1*(J1+ONE))
      IEND=IVAL-INMID
      IF (IVAL.LE.IEXC) GOTO 31

      DO 32 IJ2=2,IEND
        J1=J1-ONE
        E2=E1
        J1S=J1*J1
        JJP1=J1*(J1+ONE)
        E1=E(J1S)
        F1=F(J1,JJP1)
   32   D6J(IVAL-IJ2)=-(J1*E2*D6J(IVAL-IJ2+2)+F1*D6J(IVAL-IJ2+1))/
     $                 (E1*(J1+ONE))
C
C  MATCH UPWARD AND DOWNWARD RECURSIVE RESULTS BY SCALING.
C
   31 SCALED=SCALED/D6J(INMID)
      DO 33 IJ2=INMID,IVAL
   33   D6J(IJ2)=SCALED*D6J(IJ2)
C
C  NORMALIZE RESULTS AND SET PHASE.
C
   40 TOTAL=ZERO
      DO 41 IJ2=1,IVAL
        J1=J1MIN+DBLE(IJ2-1)
   41   TOTAL=TOTAL+(TWO*J1+ONE)*D6J(IJ2)**2
      RNORM=ONE/SQRT(TOTAL*(TWO*L1+ONE))
      IF ((SGN*D6J(IVAL)).LT.ZERO) RNORM=-RNORM
      DO 42 IJ2=1,IVAL
   42   D6J(IJ2)=D6J(IJ2)*RNORM

      RETURN
      END
