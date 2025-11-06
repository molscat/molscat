      SUBROUTINE LDPROP(N,MXLAM,NHAM,
     1                  Z,VL,IV,EINT,CENT,
     2                  RSTART,RSTOP,NSTEP,DR,NODES,
     3                  ERED,EP2RU,CM2RU,RSCALE,IPRINT)
C  Copyright (C) 2025 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE potential, ONLY: NCONST, NRSQ
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

      DIMENSION Z(N,N),VL(*),IV(*),EINT(N),CENT(N)
      ALLOCATABLE :: DIAG(:),P(:),U(:,:)

      H=DR/2.D0
      D1 = H*H/3.D0
      D2 = 2.D0*D1
      D4 = -D1/16.D0
C
      R = RSTART
      NODES=0
C
      ALLOCATE (U(N,N))
      IF (IREAD) THEN
        READ(ISCRU) U
        DO I=1,N
          U(I,I)=U(I,I)-ESHIFT
        ENDDO
      ELSE
        ALLOCATE (DIAG(N),P(MXLAM+NCONST+NRSQ))
        CALL HAMMAT(U,N,R,P,VL,IV,ERED,EINT,CENT,EP2RU,CM2RU,
     1              RSCALE,DIAG,MXLAM,NHAM,IPRINT)
        IF (IWRITE) WRITE(ISCRU) U
      ENDIF

      DO J = 1,N
        DO I = J,N
          Z(I,J) = H*Z(I,J)+D1*U(I,J)
        ENDDO
        Z(J,J) = 1.D0+Z(J,J)
      ENDDO
C
C  START PROPAGATION LOOP
C
!  Start of long DO loop #1
      DO ISTEP = 1,NSTEP
        R = R+H
        IF (IREAD) THEN
          READ(ISCRU) U
          ESH=-D4*ESHIFT
          DO I=1,N
            U(I,I)=U(I,I)+ESH
          ENDDO
        ELSE
          CALL HAMMAT(U,N,R,P,VL,IV,ERED,EINT,CENT,EP2RU,CM2RU,
     1                RSCALE,DIAG,MXLAM,NHAM,IPRINT)
          DO J = 1,N
            DO I = J,N
              U(I,J) = D4*U(I,J)
            ENDDO
            U(J,J) = 0.125D0+U(J,J)
          ENDDO
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
        DO J = 1,N
          DO I = J,N
            Z(I,J) = U(I,J)-Z(I,J)
          ENDDO
          Z(J,J) = Z(J,J)-6.D0
        ENDDO
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
          DO I=1,N
            U(I,I)=U(I,I)+ESH
          ENDDO
        ELSE
          CALL HAMMAT(U,N,R,P,VL,IV,ERED,EINT,CENT,EP2RU,CM2RU,
     1                RSCALE,DIAG,MXLAM,NHAM,IPRINT)
          DO J=1,N
            DO I=J,N
              U(I,J)=D2*U(I,J)
            ENDDO
            U(J,J)=U(J,J)+2.D0
          ENDDO
          IF (IWRITE) WRITE(ISCRU) U
        ENDIF
C
C  EQ 10 FROM 1973 WITH U=2+D2*V
C
        DO J = 1,N
        DO I = J,N
          Z(I,J) = U(I,J)-Z(I,J)
        ENDDO
        ENDDO
      ENDDO
!  End of long DO loop #1

      IF (.NOT.IREAD) DEALLOCATE (DIAG,P)
      DEALLOCATE (U)
C
C  FINISHED PROPAGATING. CONVERT 1+Z BACK TO Y WITH EQ 21 FROM 1977
C
      HI = 1.D0/H
      DO J = 1,N
        DO I = J,N
          Z(I,J) = HI*Z(I,J)
          Z(J,I) = Z(I,J)
        ENDDO
        Z(J,J) = Z(J,J)-HI
      ENDDO

      RETURN
C
  900 WRITE(6,901)
  901 FORMAT(/' *** ERROR IN SYMINV CALLED FROM LDPROP - TERMINATING')
      STOP
      END
