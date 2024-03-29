      SUBROUTINE LDPROP(N,MXLAM,NHAM,
     1                  Z,U,VL,IV,EINT,CENT,P,DIAG,
     2                  RSTART,RSTOP,NSTEP,DR,NODES,
     3                  ERED,EP2RU,CM2RU,RSCALE,IPRINT)
C  Copyright (C) 2022 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  SUBROUTINE FOR JOHNSON'S LOG-DERIVATIVE PROPAGATOR
C  SEE B R JOHNSON, J COMP PHYS 13, 445 (1973)
C      B R JOHNSON, J CHEM PHYS 67, 4086 (1977)
C  THE ORIGINS OF THIS ROUTINE ARE LOST IN THE MISTS OF TIME
C
C  COMMENTS AND SOME F77 CONSTRUCTS ADDED BY J M HUTSON 2006
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  COMMON BLOCK FOR CONTROL OF USE OF PROPAGATION SCRATCH FILE
      LOGICAL IREAD,IWRITE,IREADR,IWRITR
      COMMON /PRPSCR/ ESHIFT,ISCRU,ISCRUR,IREAD,IWRITE,IREADR,IWRITR

      DIMENSION U(N,N),Z(N,N),P(MXLAM),VL(2),IV(2),EINT(N),CENT(N),
     1          DIAG(N)

      H=DR/2.D0
      D1 = H*H/3.D0
      D2 = 2.D0*D1
      D4 = -D1/16.D0
C
      R = RSTART
      NODES=0
C
      IF (IREAD) THEN
        READ(ISCRU) U
        DO 130 I=1,N
  130     U(I,I)=U(I,I)-ESHIFT
      ELSE
        CALL HAMMAT(U,N,R,P,VL,IV,ERED,EINT,CENT,EP2RU,CM2RU,
     1              RSCALE,DIAG,MXLAM,NHAM,IPRINT)
        IF (IWRITE) WRITE(ISCRU) U
      ENDIF

      DO 150 J = 1,N
        DO 140 I = J,N
  140     Z(I,J) = H*Z(I,J)+D1*U(I,J)
  150   Z(J,J) = 1.D0+Z(J,J)
C
C  START PROPAGATION LOOP
C
      DO 260 ISTEP = 1,NSTEP
        R = R+H
        IF (IREAD) THEN
          READ(ISCRU) U
          ESH=-D4*ESHIFT
          DO 160 I=1,N
  160       U(I,I)=U(I,I)+ESH
        ELSE
          CALL HAMMAT(U,N,R,P,VL,IV,ERED,EINT,CENT,EP2RU,CM2RU,
     1                RSCALE,DIAG,MXLAM,NHAM,IPRINT)
          DO 180 J = 1,N
            DO 170 I = J,N
  170         U(I,J) = D4*U(I,J)
  180       U(J,J) = 0.125D0+U(J,J)
          IF (IWRITE) WRITE(ISCRU) U
        ENDIF
C
C  INVERT U = (1/8)[1-(H**2/6)*VMAT] FOR EQ 7 OF JOHNSON (1973)
C
        CALL SYMINV(U,N,N,NCU)
        IF (NCU.GT.N) GOTO 900
C
C  INVERT 1+Z FOR EQ 10 OF JOHNSON (1973)
C
        CALL SYMINV(Z,N,N,NCZ)
        IF (NCZ.GT.N) GOTO 900
        NODES=NODES+NCZ
C
C  FIRST HALF OF SECTOR:
C  EQ 10 FROM 1973 WITH W=4 AND U FROM INVERSION
C  USING (I+Z)^{-1} Z    = I - (I+Z)^{-1}
C  AND   (I-H2*V)^{-1} V = (1/H2) [I - (I - H2*V)^{-1}]
C
        DO 210 J = 1,N
          DO 200 I = J,N
  200       Z(I,J) = U(I,J)-Z(I,J)
  210     Z(J,J) = Z(J,J)-6.D0
C
        CALL SYMINV(Z,N,N,NCZ)
        IF (NCZ.GT.N) GOTO 900
        NODES=NODES+NCZ-NCU
C
C  SECOND HALF OF SECTOR WITH W=2, OR W=1 FOR LAST POINT
C
        R = R+H
        IF (ISTEP.EQ.NSTEP) D2=D1
        IF (IREAD) THEN
          READ(ISCRU) U
          ESH=-D2*ESHIFT
          DO 220 I=1,N
  220       U(I,I)=U(I,I)+ESH
        ELSE
          CALL HAMMAT(U,N,R,P,VL,IV,ERED,EINT,CENT,EP2RU,CM2RU,
     1                RSCALE,DIAG,MXLAM,NHAM,IPRINT)
          DO 240 J=1,N
            DO 230 I=J,N
  230         U(I,J)=D2*U(I,J)
  240       U(J,J)=U(J,J)+2.D0
          IF (IWRITE) WRITE(ISCRU) U
        ENDIF
C
C  EQ 10 FROM 1973 WITH U=2+D2*V
C
        DO 250 J = 1,N
        DO 250 I = J,N
  250     Z(I,J) = U(I,J)-Z(I,J)
  260 CONTINUE
C
C  FINISHED PROPAGATING. CONVERT 1+Z BACK TO Y WITH EQ 21 FROM 1977
C
      HI = 1.D0/H
      DO 280 J = 1,N
        DO 270 I = J,N
          Z(I,J) = HI*Z(I,J)
  270     Z(J,I) = Z(I,J)
  280   Z(J,J) = Z(J,J)-HI

      RETURN
C
  900 WRITE(6,901)
  901 FORMAT(/' *** ERROR IN SYMINV CALLED FROM LDPROP - TERMINATING')
      STOP
      END
