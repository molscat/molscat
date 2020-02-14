      module BASE9_SUITE
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  CR Le Sueur Oct 2018
C  This file provides a skeleton version of a plug-in basis-set suite
C  with some indications of how to adapt it for your own use

C  It must be read in conjunction with the full documentation,
C  which specifies the values that must be given to the variables
C  and arrays that are set in these routines.

      IMPLICIT NONE

c  quantities available to all routines in this basis-set suite
      double precision :: xxxxx
      integer          :: ixxxx
c  nqns and nlabvs are used to set variables nqn and nlabv that are
c  returned from set9
c  nqns is one greater than the number of quantum labels for a pair state
c  nlabvs is the number of indices for each term in the potential expansion
c  nqn is passed back into various routines in the basis-set suite
c  nlabvs is used internally in potin9
      integer, parameter        :: nqns=2, nlabvs=1
      end module BASE9_SUITE
C=========================================================================
      SUBROUTINE BAS9IN(PRTP, IBOUND, IPRINT)
      USE potential
      USE efvs  ! only needed if efvs are included in the calculation
      USE base9_suite
c     USE basis_data
      IMPLICIT NONE

      CHARACTER(32), INTENT(OUT)   :: PRTP

      INTEGER,       INTENT(INOUT) :: IBOUND

      INTEGER,       INTENT(IN)    :: IPRINT

      integer :: iconst,iextra,iefv

c  bas9in may make use of input variables passed in module basis_data
c  by uncommenting the USE line above and may also read any additional
c  variables required, for example in a namelist block such as
c     NAMELIST / BASIS9 / extra_variables

      PRTP = 'collision type'
      IBOUND = 0

C  NDGVL, NCONST, NRSQ and NEXTRA are in module potential

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
      NRSQ = 0

c  needed only for MOLSCAT and only if extra operators are required
c  to resolve degeneracies in H_intl
      NEXTRA = 0 ! number of extra operators

      if (nextra.gt.0) then
        do iextra=1,nextra
          nextms(iextra) = 0 ! number of coupling matrices used for each
                             ! extra operator
        enddo
      endif

c  needed only if efvs are included in the calculation
      NEFV = 0 ! number of external field variables

      if (nefv.gt.0) then
        do iefv=1,nefv
          efvnam(iefv) = 'efv name' ! name for efv
          efvunt(iefv) = 'units'    ! units for efv
        enddo
        mapefv = 1 ! position of first efv in vconst
      endif

c  read any additional quantities required to define the basis set
c     READ(5,&BASIS9)

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

      integer :: iqn,iqn1,iqn1mn,iqn1mx,iqn1st,ilevel,iloop,istate,
     1           nqlev

      NQN = nqns ! one greater than the number of quantum labels
      do iqn=1,nqn-1
        qname(iqn) = 'qn label' ! name of quantum label
      enddo

c  number of symmetry blocks for each JTOT
      NBLOCK = 1

c  number of indices for each term in the potential expansion
      NLABV = nlabvs

      if (levin) then
        write(6,*) 'jlevel set in &basis namelist'
        stop
      endif

c  set jlevel, elevel and nlevel only if H_intl is diagonal

c  nqlev is the number of quantum numbers that affect the pair energy.
c  it is internal to set9 and not used elsewhere.
      nqlev=2

      if (nconst.eq.0) then
        iqn1mn=jmin
        iqn1mx=jmax
        iqn1st=jstep
        ilevel=0
        do iqn1=iqn1mn,iqn1mx,iqn1st
c       do iqn2=iqn2mn,iqn2mx,iqn2st ! other quantum labels to cycle over
          jlevel(1+nqlev*ilevel)=iqn1
c         jlevel(2+nqlev*ilevel)=iqn2
          ilevel=ilevel+1
          if (.not.ein) then
            elevel(ilevel) = 0.D0 ! internal energy as function of quantum numbers
          endif
        enddo
        nlevel=ilevel
      endif

c  count the number of states and then assign values to the jstate array.
c  this example (with loops over iqn2 and iqncpl uncommented) is for
c  two angular momenta iqn1 and iqn2 that couple to give a resultant iqncpl
      do iloop=1,2
        istate=0
        do iqn1=iqn1mn,iqn1mx,iqn1st
c       do iqn2=iqn2mn,iqn2mx,iqn2st       ! other quantum labels to cycle over
c       do iqncpl=abs(iqn1-iqn2),iqn1+iqn2 ! and maybe some coupling between them
          istate=istate+1
          if (iloop.eq.2) then
            jstate(istate)=iqn1
c           jstate(istate+nstate)=iqn2     ! other quantum labels stored
c           jstate(istate+nstate*2)=iqncpl ! and coupling label also stored
          endif
c       enddo
c       enddo
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

      integer :: ifunc,istate,iqn1,lmin,lmax,ll

      ifunc=0
      do istate=1,nstate
        iqn1=jstate(istate,1)
c       iqn2=jstate(istate,2)
c       ... etc
        lmin=0
        do ll=lmin,lmax
c         if (.not.(some conditions on jtot, iblock, ll, iqn1 (etc))) cycle

c  now counting only those basis functions included in the current symmetry block
          ifunc=ifunc+1
          if (.not.lcount) then
            jsindx(ifunc)=istate ! pointer to basis function
            l(ifunc)=ll          ! value of L for basis function
          endif
        enddo
      enddo
      n=ifunc

      RETURN
      END SUBROUTINE BASE9
C=========================================================================
      SUBROUTINE POTIN9(ITYPP, LAM, MXLAM, NPTS, NDIM, XPT, XWT, MXPT,
     1                  IVMIN, IVMAX, L1MAX, L2MAX, MXLMB, XFN, MX,
     2                  IXFAC)
      USE base9_suite
      IMPLICIT NONE

      INTEGER, INTENT(INOUT)          :: ITYPP, MXLAM

      INTEGER, INTENT(OUT)            :: LAM(*)

      LOGICAL :: LVRTP
      DATA LVRTP/.FALSE./

      INTEGER, INTENT(IN)             :: IVMIN, IVMAX, L1MAX, L2MAX

      INTEGER :: L1

c  the remaining quantities are used only if quadrature is to be used
c  to project out potential expansion coefficients

      INTEGER, INTENT(INOUT)          :: NDIM, NPTS(NDIM), IXFAC, MX

      DOUBLE PRECISION, INTENT(INOUT) :: XFN(*)

      DOUBLE PRECISION, INTENT(OUT)   :: XPT(MXPT,NDIM), XWT(MXPT,NDIM)

      INTEGER, INTENT(IN)             :: MXLMB, MXPT

      double precision, allocatable :: fn1(:,:) !,fn2(:,:) etc
      integer :: il,npttot,nfun,ipt,ipt1tm,ipt1,i,ix
      double precision pt1,wt1

c  optionally set itypp to one of the precoded values (1 to 8) and do nothing else
      itypp = 1
c  alternatively define a special-purpose set of potential indices for itypp=9
      if (itypp.eq.9) then
        mxlam=0
        do l1=0,l1max
c       do l2=0,l2max
c       do l3=0,l2max
c       do lcpl=abs(11-l2),l1+l2 ! and maybe some coupling between them
          lam(1+mxlam*nlabvs)=l1
c         lam(2+mxlam*nlabvs)=l2
c         lam(3+mxlam*nlabvs)=l3
          mxlam=mxlam+1
c       enddo
c       enddo
        enddo
      endif

c  if potential coefficients are to be obtained by quadrature and itypp=9,
c  the points and weights must be defined here.
c  This is complicated and seldom necessary, but the code might look something like this:

      if (itypp.eq.9 .and. lvrtp) then

c  set values for npts if not set in namelist
c       npts(1)=

c  zero everything to start with
        xpt(1:npts(1),1)=0.D0
        xwt(1:npts(1),1)=0.D0
c  generate set of quadrature points and weights
c       call getqpt(1,xpt(1,1),xwt(1,1))

c  set up an array to hold the values of the functions for this quadrature
        allocate (fn1(0:l1max,npts(1)))

c  generate sets of functions at each quadrature point
c       call quad1(xpt(1,1),fn1,l1max)

c  do the same for all remaining quadratures
        if (ndim.gt.1) then
c         npts(2)=
c         call getqpt(2,xpt(1,2),xwt(1,2))
c         allocate (fn2(0:l2max,npts(2)))
c         call quad2(xpt(1,2),fn2,l2max)
        endif

c  calculate the total number of points in the quadrature
        npttot=product(npts(1:ndim))

c  calculate the total number of evaluated functions
        nfun=npttot*mxlam

c  set where the evaluated functions are going to be stored
c  (at the end of the x array)
        ixfac=mx-nfun

c  set the start point for storing evaluated functions
        ix=ixfac

c  cycle over the total number of points in the quadrature
        do ipt=1,npttot

c  work out the point number for the first quadrature
          ipt1tm=ipt
          ipt1=mod(ipt,npts(1))
c  and get the point and weight corresponding to that value
          pt1=xpt(ipt1,1)
          wt1=xwt(ipt1,1)

c  do the same for all remaining quadratures
c         ipt2tm=(ipt1tm-1)/npts(1)+1
c         ipt2=mod(ipt2tm,npts(2))
c         pt2=xpt(ipt2,2)
c         wt2=xwt(ipt2,2)

c  cycle over the potential expansion
          do i=1,mxlam

c  get the index for the first label
            l1=lam(i)

c  get indices for all the other labels
c           l2=lam(i+mxlam)

c  set the index for this evaluated function
            ix=ix+1

c  calculate the evaluated function and store it
            xfn(ix)=fn1(l1,ipt1)*wt1 !*fn2(l2,ipt2)*wt2 etc...
          enddo
        enddo

c  reset the size of the x array so these evaluated functions can't be overwritten
        mx=ixfac

c  these arrays are no longer needed
        deallocate (fn1)
c       deallocate (fn2) ! etc
      endif

      RETURN
      END SUBROUTINE POTIN9
C=========================================================================
      SUBROUTINE CPL9(N, IBLOCK, NHAM, LAM, MXLAM, NSTATE, JSTATE,
     1                JSINDX, L, JTOT, VL, IV, CENT, DGVL, IBOUND,
     2                IEXCH, IPRINT)
      USE base9_suite
      USE potential
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: VL(NVLBLK,N*(N+1)/2), CENT(N),
     1                                 DGVL(N,NDGVL)

      INTEGER, INTENT(OUT)          :: IV(NVLBLK,N*(N+1)/2)

      INTEGER, INTENT(IN)           :: N, IBLOCK, NHAM, MXLAM,
     1                                 NSTATE, JSTATE(NSTATE,NQNS),
     2                                 JSINDX(N), L(N), JTOT, IBOUND,
     3                                 IEXCH, IPRINT, LAM(*)

      integer :: irc,icol,irow,i1col,i1row,lcol,lrow,iham,idgvl
      logical livuse
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
        i1col=jstate(jsindx(icol),1) ! other quantum labels may also be needed
        lcol=l(icol)
        DO IROW = 1, ICOL
          i1row=jstate(jsindx(irow),1) ! other quantum labels may also be needed
          lrow=l(irow)
          IRC = IRC + 1
          DO IHAM = 1, NVLBLK

c  the following will usually involve calculations that depend on the quantum labels
c  of the row and column and of the potential (or H_intl) expansion term
            if (iham.le.mxlam) then
              VL(IHAM, IRC) = 0.D0 ! the iham-th term in the potential
            elseif (iham.le.mxlam+nconst) then
              vl(iham, irc) = 0.D0 ! the iham-mxlam-th term in H_intl
            elseif (iham.le.mxlam+nconst+nrsq) then
              vl(iham, irc) = 0.D0 ! the iham-mxlam-nconst-th term in L^2
            else
              vl(iham, irc) = 0.D0 ! the iham-mxlam-nconst-nrsq-th
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
