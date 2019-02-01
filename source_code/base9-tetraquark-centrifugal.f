      module BASE9_SUITE
C  Copyright (C) 2018 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  J. M. Hutson, January 2019
C  This plug-in basis set suite is for a system of two heavy quarks
C  and two light quarks, coded with coupled Born-Oppenheimer potentials

      IMPLICIT NONE
c  quantities available to all routines in this basis-set suite
      integer          :: LLTL,LLTU,ISLT,ISHV
c  nqns and nlabvs are used to set variables nqn and nlabv that are
c  returned from set9
c  nqns is one greater than the number of quantum labels for a pair state
c  nlabvs is the number of indices for each term in the potential expansion
c  nqn is passed back into various routines in the basis-set suite
c  nlabvs is used internally in potin9
      integer, parameter        :: nqns=7, nlabvs=3
      end module BASE9_SUITE
C=========================================================================
      SUBROUTINE BAS9IN(PRTP, IBOUND, IPRINT)
      USE potential
      USE efvs  ! only needed if efvs are included in the calculation
c     USE basis_data
      USE base9_suite
      IMPLICIT NONE

      CHARACTER(32), INTENT(OUT)   :: PRTP

      INTEGER,       INTENT(INOUT) :: IBOUND

      INTEGER,       INTENT(IN)    :: IPRINT

      integer :: iconst,iextra,iefv

c  bas9in may read any input variables needed to define the basis set
c  that are not read by base and passed in module basis_data
      NAMELIST / BASIS9 / LLTL,LLTU,ISLT,ISHV

      PRTP = 'tetraquark with only centrifugal'
      IBOUND = 0

C  these quantities are included in the potential module

c  +ve only if H_intl is diagonal and some terms contributing to the
c  internal energy need to be recalculated during the course of a run
      NDGVL = 0 ! number of terms contributing to internal energy

c  +ve only if H_intl is non-diagonal
      NCONST = 0 ! number of terms in H_intl
      if (nconst.gt.0) then
        do iconst=1,nconst
          VCONST(iconst) = 0.D0 ! only set those expansion coefficients
                                ! which do not need to be re-evaluated
                                ! during the course of a run.
        enddo
      endif

c  only +ve if L^2 is non-diagonal
      NRSQ = 1

c  needed (for MOLSCAT) only if extra operators are required to resolve
c  degeneracies in H_intl
      NEXTRA = 0 ! number of extra operators

      do iextra=1,nextra
        nextms(iextra) = 0 ! number of coupling matrices used for each
                           ! extra operator
      enddo

c  needed only if efvs are included in the calculation
      NEFV = 0 ! number of external field variables

      do iefv=1,nefv
        efvnam(iefv) = 'efv name' ! name for efv
        efvunt(iefv) = 'units'    ! units for efv
      enddo
      mapefv = 1 ! position of first efv in vconst

      LLTL=1
      LLTU=1
      ISLT=0
      ISHV=0
      READ(5,BASIS9)

      WRITE(6,601) LLTL,LLTU,ISLT,ISHV
  601 FORMAT('  TETRAQUARK BASIS SET ROUTINE BY J. M. HUTSON, JAN 2019'/
     1 '  BORN-OPPENHEIMER BASIS WITH L_light FROM',I2,' TO',I2/
     2 '  AND ALL ALLOWED BODY-FIXED PROJECTIONS LAMBDA OF L_light'//
     3 '  FIXED SPINS S_light =',I2,' AND S_heavy =',I2//
     4 '  JTOT IS TOTAL ANGULAR MOMENTUM (el IF ALL SPINS ZERO)'/
     5 '  SYM BLOCK 2 IS EVEN COMBINATION OF +|LAMBDA| AND -|LAMBDA|'/
     6 12X,'1 IS  ODD COMBINATION AND EXCLUDES LAMBDA = 0'//
     7 '  BASIS SET [({LTOT L_light LAMBDA} S_light) N S_heavy] JTOT'//
     8 '  THIS INITIAL VERSION INCLUDES ONLY CENTRIFUGAL COUPLING',
     9 ' BETWEEN LAMBDA STATES AND DOES NOT INCLUDE:'/
     A '  - NON-CENTRIFUGAL COUPLING BETWEEN LAMBDA STATES'/
     B '  - COUPLING BETWEEN DIFFERENT VALUES OF L_light'/
     C '  - COUPLING BETWEEN DIFFERENT VALUES OF LTOT OR N'/
     D '  - COUPLING BETWEEN DIFFERENT VALUES OF S_light OR S_heavy')

      RETURN
      END SUBROUTINE BAS9IN
C=========================================================================
      SUBROUTINE SET9(LEVIN, EIN, NSTATE, JSTATE, NQN, QNAME, NBLOCK,
     1                NLABV, IPRINT)
      USE basis_data ! needed only if H_intl is diagonal
      USE potential  ! needed only to check NCONST
      USE base9_suite
      IMPLICIT NONE

      INTEGER, INTENT(OUT)      :: NQN, NBLOCK, NLABV, NSTATE,
     1                             JSTATE(*)
c  JSTATE is dimensioned (NSTATE,NQN), but because NSTATE is not known
c  on entry, it cannot be given those dimensions explicitly.
c
      CHARACTER(8), INTENT(OUT) :: QNAME(10)

      LOGICAL, INTENT(IN)       :: LEVIN, EIN

      INTEGER, INTENT(IN)       :: IPRINT

      integer :: iqn,iqn1,iqn2,iqn3,iqn4,iqn5,iqn6,ilevel,iloop,istate

      NQN = nqns ! one greater than the number of quantum labels for pair states
      qname(1) = 'L_light'
      qname(2) = 'Lambda'
      qname(3) = 'S_light'
      qname(4) = 'LTOT-N'
      qname(5) = 'S_heavy'
      qname(6) = 'N-JTOT'

c  number of symmetry blocks for each JTOT
      NBLOCK = 2

c  number of quantum labels for a potential term
C  For now Lambda, L_light, S_light
      NLABV = nlabvs

      if (levin) then
        write(6,*) 'jlevel set in &basis namelist'
        stop
      endif

c  set jlevel, elevel and nlevel only if H_intl is diagonal
      if (nconst.eq.0) then
c  nqlev is the number of quantum numbers that affect the pair energy.
c  it is internal to set9 and not used elsewhere.
c       nqlev=2
        ilevel=0
        do iqn1=LLTL,LLTU
        do iqn2=0,iqn1
c         jlevel(1+nqlev*ilevel)=iqn1
c         jlevel(2+nqlev*ilevel)=iqn2
c         ilevel=ilevel+1
          if (.not.ein) then
            elevel(ilevel) = 0.D0 ! internal energy as function of quantum numbers
          endif
        enddo
        enddo
        nlevel=ilevel
      endif

c  count the number of states and then assign values to the jstate array

c  this is based on basis functions where LTOT = L_light + L_heavy
c  but there is strong coupling to the inter-heavy axis so both
c  LTOT and L_light have projection Lambda onto the axis
c  LTOT then couples to S_light to give resultant N
c  N couples to S_heavy to give resultant JTOT

c  The allowed values of N and LTOT depend on JTOT:
c  to store them in a JTOT-independent way, the values placed in
c  iqn4 and iqn6 and LTOT-N and N-JTOT, respectively.

c  In this initial base9, ISLT and ISHV are both 0, so LTOT = N = JTOT.
c
      do iloop=1,2
        istate=0
        do iqn1=LLTL,LLTU                  ! L light
        do iqn2=0,iqn1                     ! Lambda
        iqn3=ISLT                          ! S light
        do iqn4=-ISLT,ISLT                 ! LTOT-N
        iqn5=ISHV                          ! S heavy
        do iqn6=-ISHV,ISHV                 ! N-JTOT
          istate=istate+1
          if (iloop.eq.2) then
            jstate(istate)=iqn1
            jstate(istate+nstate)=iqn2
            jstate(istate+nstate*2)=iqn3
            jstate(istate+nstate*3)=iqn4
            jstate(istate+nstate*4)=iqn5
            jstate(istate+nstate*5)=iqn6
          endif
        enddo
        enddo
        enddo
        enddo
        nstate=istate
      enddo

      RETURN
      END SUBROUTINE SET9
C=========================================================================
      SUBROUTINE BASE9(LCOUNT, N, JTOT, IBLOCK, JSTATE, NSTATE, NQN,
     1                 JSINDX, L, IPRINT)
      USE base9_suite
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: N

      INTEGER, INTENT(OUT)   :: JSINDX(N), L(N)

      LOGICAL, INTENT(IN)    :: LCOUNT

      INTEGER, INTENT(IN)    :: JTOT, IBLOCK, NSTATE, NQN,
     1                          JSTATE(NSTATE,NQN), IPRINT

      integer :: ifunc,istate,lambda,ltot

      ifunc=0
      do istate=1,nstate
        LAMBDA=jstate(istate,2)
        IF (LAMBDA.EQ.0 .AND. IBLOCK.EQ.1) CYCLE
        LTOT=JTOT+jstate(istate,4)+jstate(istate,6)
        if (LAMBDA.GT.LTOT) CYCLE
c  now counting only those basis functions included in the current symmetry block
          ifunc=ifunc+1
          if (.not.lcount) then
            jsindx(ifunc)=istate ! pointer to basis function
            l(ifunc)=0           ! not used for NRSQ>0
          endif
      enddo
      n=ifunc

      RETURN
      END SUBROUTINE BASE9
C=========================================================================
      SUBROUTINE POTIN9(ITYPP, LAM, MXLAM, NPTS, NDIM, XPT, XWT, MXPT,
     1                  IVMIN, IVMAX, L1MAX, L2MAX, MXLMB, XFN, MX,
     2                  IXFAC)
      USE base9_suite, ONLY : LLTL, LLTU, ISLT, NLABVS
      IMPLICIT NONE

      INTEGER, INTENT(INOUT)          :: ITYPP, MXLAM

      INTEGER, INTENT(OUT)            :: LAM(*)

c  the quantities below are used only if quadrature is to be used
c  to project out potential expansion coefficients
      INTEGER, INTENT(INOUT)          :: NDIM, NPTS(NDIM), IXFAC, MX

      DOUBLE PRECISION, INTENT(INOUT) :: XFN(*)

      DOUBLE PRECISION, INTENT(OUT)   :: XPT(MXPT,NDIM), XWT(MXPT,NDIM)

      INTEGER, INTENT(IN)             :: MXPT, IVMIN, IVMAX, L1MAX,
     1                                   L2MAX, MXLMB

      integer :: l1,l2,nlabv
      LOGICAL :: LVRTP
      DATA LVRTP/.FALSE./

      itypp = 9

      if (itypp.eq.9) then
        mxlam=0
        do l1=LLTL,LLTU
        do l2=0,l1
          lam(1+mxlam*nlabvs)=l1
          lam(2+mxlam*nlabvs)=l2
          lam(3+mxlam*nlabvs)=ISLT
          mxlam=mxlam+1
        enddo
        enddo
      endif

      RETURN
      END SUBROUTINE POTIN9
C=========================================================================
      SUBROUTINE CPL9(N, IBLOCK, NPOTL, LAM, MXLAM, NSTATE, JSTATE,
     1                JSINDX, L, JTOT, VL, IV, CENT, DGVL, IBOUND,
     2                IEXCH, IPRINT)
      USE base9_suite
      USE potential
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: VL(NVLBLK,N*(N+1)/2), CENT(N),
     1                                 DGVL(N,NDGVL)

      INTEGER, INTENT(OUT)          :: IV(NVLBLK,N*(N+1)/2)

      INTEGER, INTENT(IN)           :: N, IBLOCK, NPOTL, MXLAM,
     1                                 NSTATE, JSTATE(NSTATE,NQNS),
     2                                 JSINDX(N), L(N), JTOT, IBOUND,
     3                                 IEXCH, IPRINT, LAM(*)

      integer :: irc,ipotl,idgvl,lcol,lrow
      integer :: icol,i1col,i2col,i3col,i4col,i5col,i6col
      integer :: irow,i1row,i2row,i3row,i4row,i5row,i6row
      logical :: livuse
      data livuse/.false./

c  Option to use IV array not supported in base9.
c  This code is included to silence compiler warning
      if (livuse) then
        iv(1,1) = 0
      endif

      DO ICOL = 1, N

c  if nconst = 0, there can be contributions to the internal energy that
c  need to be recalculated during the course of the calculation,
c  for example to take account of diagonal terms involving external fields.
c  These are included in the DGVL array and there are NDGVL of them.
        if (nconst.eq.0 .and. ndgvl.gt.0) then
          do idgvl=1,ndgvl
            dgvl(n,idgvl) = 0.D0 ! the idgvl-th contribution to the
                                 ! internal energy
          enddo
        endif

c  if ibound /= 0, the centrifugal operator is not simply L(L+1).
c  if nrsq = 0, its diagonal elements are stored in CENT.
        lcol=l(icol)
        if (ibound.ne.0 .and. nrsq.eq.0) then
          cent(icol)=0.D0 ! the diagonal value of the L^2 operator for
                          ! this basis function
        endif
      enddo

c  calculate the coupling matrices and store them in the VL array
      IRC = 0
      DO ICOL = 1, N
        i1col=jstate(jsindx(icol),1)
        i2col=jstate(jsindx(icol),2)
        i3col=jstate(jsindx(icol),3)
        i4col=jstate(jsindx(icol),4)
        i5col=jstate(jsindx(icol),5)
        i6col=jstate(jsindx(icol),6)
        lcol=jtot+i4col+i6col
        DO IROW = 1, ICOL
          i1row=jstate(jsindx(irow),1)
          i2row=jstate(jsindx(irow),2)
          i3row=jstate(jsindx(irow),3)
          i4row=jstate(jsindx(irow),4)
          i5row=jstate(jsindx(irow),5)
          i6row=jstate(jsindx(irow),6)
          lrow=jtot+i4row+i6row
          IRC = IRC + 1
          DO IPOTL = 1, NVLBLK

c  the following will usually involve calculations that depend on the quantum labels
c  of the row and column and of the potential (or H_intl) expansion term
            if (ipotl.le.mxlam) then
              vl(ipotl, irc) = 0.D0 ! the ipotl-th term in the potential
c  the ipotl-mxlam-nconst-th term in the interaction potential: diagonal if qns match
              if (i1row.eq.i1col .and. i2row.eq.i2col .and.
     1            i3row.eq.i3col .and. i4row.eq.i4col .and.
     2            i5row.eq.i5col .and. i6row.eq.i6col) then
                   if (i1row.eq.lam(1+3*(ipotl-1)) .and.
     1                 i2row.eq.lam(2+3*(ipotl-1)) .and.
     2                 i3row.eq.lam(3+3*(ipotl-1)))
     3                 vl(ipotl, irc) = 1.D0
              endif
            elseif (ipotl.le.mxlam+nconst) then ! not used here because NCONST=0
              vl(ipotl, irc) = 0.D0 ! the ipotl-mxlam-th term in H_intl
            elseif (ipotl.le.mxlam+nconst+nrsq) then
              vl(ipotl, irc) = 0.D0
c  the ipotl-mxlam-nconst-th term in L^2 ! zero unless diagonal in all but i2 (Lambda)
              if (i1row.eq.i1col .and.
     1            i3row.eq.i3col .and. i4row.eq.i4col .and.
     2            i5row.eq.i5col .and. i6row.eq.i6col) then
                if (i2row.eq.i2col) then
                  vl(ipotl, irc) =
     1              dble(lrow*(lrow+1) + i1row*(i1row+1) - 2*i2row**2) ! diagonal
                else
                  if (abs(i2row-i2col).eq.1) vl(ipotl, irc) =
     1              sqrt(dble((i1row*(i1row+1)-i2row*i2col)
     2                       *(lrow*(lrow+1)-i2row*i2col)))          ! Coriolis
                  if (i2row.eq.0 .or. i2col.eq.0)
     1              vl(ipotl, irc) = sqrt(2.D0) * vl(ipotl, irc)     ! Sym factor
                endif
              endif
            else                                ! not used here because NEXTRA=0
              vl(ipotl, irc) = 0.D0 ! the ipotl-mxlam-nconst-nrsq-th
                                    ! term used for extra operators
            endif
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CPL9
C=========================================================================
c     SUBROUTINE DEGEN9(JJ1, JJ2, DEGFAC)
c     USE base9_suite
c     IMPLICIT NONE

c     DOUBLE PRECISION, INTENT(OUT) :: DEGFAC

c     INTEGER, INTENT(IN)           :: JJ1, JJ2

c     RETURN
c     END SUBROUTINE DEGEN9
C=========================================================================
c     SUBROUTINE THRSH9(IREF, MONQN, NQN, EREF, IPRINT)
c     USE base9_suite
c     IMPLICIT NONE

c     DOUBLE PRECISION, INTENT(OUT) :: EREF

c     INTEGER, INTENT(IN)           :: IREF, MONQN(NQN), NQN, IPRINT

c     RETURN
c     END SUBROUTINE THRSH9
C=========================================================================
c     SUBROUTINE EFV9(IFVARY)
c     USE efvs
c     USE potential
c     USE base9_suite
c     IMPLICIT NONE

c     INTEGER, INTENT(IN) :: IFVARY

c     RETURN
c     END SUBROUTINE EFV9
