      module BASE9_SUITE
C  Copyright (C) 2018 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  CR Le Sueur Oct 2018
C  This file provides a skeleton version of a plug-in basis set suite
c  with some indications of how to adapt it for your own use
      IMPLICIT NONE

c  names of quantities used by other routines in this basis set suite
      double precision :: xxxxx
      integer          :: ixxxx
      integer, parameter        :: nqns=2
      end module BASE9_SUITE
C=========================================================================
      SUBROUTINE BAS9IN(PRTP, IBOUND, IPRINT)
      USE potential
      USE efvs  ! only needed if efvs are included in the calculation
      USE base9_suite
      IMPLICIT NONE

      CHARACTER(32), INTENT(OUT)   :: PRTP

      INTEGER,       INTENT(INOUT) :: IBOUND

      INTEGER,       INTENT(IN)    :: IPRINT

      integer :: iconst,iextra,iefv

      PRTP = 'collision type'
      IBOUND = 0

C  these quantities are included in the potential module

c  only +ve if H_intl is diagonal and some terms contributing to the
c  internal energy need to be recalculated during the course of a run
      NDGVL = 0 ! number of terms contributing to internal energy

c  only +ve if H_intl is non-diagonal
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

c  only needed (for MOLSCAT) if extra operators are required to resolve
c  degeneracies in H_intl
      NEXTRA = 0 ! number of extra operators

      if (nextra.gt.0) then
        do iextra=1,nextra
          nextms(iextra) = 0 ! number of coupling matrices used for each
                             ! extra operator
        enddo
      endif

c  only needed if efvs are included in the calculation
      NEFV = 0 ! number of external field variables

      if (nefv.gt.0) then
        do iefv=1,nefv
          efvnam(iefv) = 'efv name' ! name for efv
          efvunt(iefv) = 'units'    ! units for efv
        enddo
        mapefv = 1 ! position of first efv in vconst
      endif

      RETURN
      END SUBROUTINE BAS9IN
C=========================================================================
      SUBROUTINE SET9(LEVIN, EIN, NSTATE, JSTATE, NQN, QNAME, NBLOCK,
     1                NLABV, IPRINT)
      USE basis_data ! only needed if H_intl is diagonal
      USE potential  ! only needed to check NCONST
      USE base9_suite
      IMPLICIT NONE

      INTEGER, INTENT(OUT)      :: NQN, NBLOCK, NLABV, NSTATE,
     1                             JSTATE(*)
c  JSTATE is dimensioned (NSTATE,NQN), but because NSTATE is not known
c  on entry, we can't give it those dimensions explicitly.
c
      CHARACTER(8), INTENT(OUT) :: QNAME(10)

      LOGICAL, INTENT(IN)       :: LEVIN, EIN

      INTEGER, INTENT(IN)       :: IPRINT

      integer :: iqn,iqn1,iqn1mn,iqn1mx,iqn1st,ilevel,iloop,istate

      NQN = nqns ! one more than the number of quantum labels
      do iqn=1,nqn-1
        qname(iqn) = 'qn label' ! name of quantum label
      enddo

      NBLOCK = 1

      NLABV = 1

      if (levin) then
        write(6,*) 'jlevel set in &basis namelist'
        stop
      endif

c  set jlevel, elevel and nlevel only if H_intl is diagonal
      if (nconst.eq.0) then
        iqn1mn=jmin
        iqn1mx=jmax
        iqn1st=jstep
        ilevel=0
        do iqn1=iqn1mn,iqn1mx,iqn1st
c       do iqn2=iqn2mn,iqn2mx,iqn2st ! other quantum labels to cycle over
          ilevel=ilevel+1
          jlevel(ilevel)=iqn1
c  not sure about this if statement...
          if (.not.ein) then
            elevel(ilevel)=roti(1)*iqn1*(iqn1+1)
          endif
        enddo
        nlevel=ilevel
      endif

c  count the number of states and then assign values to the jstate array
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
c       iqn2=jstate(istate,2) ... etc
        lmin=0
        do ll=lmin,lmax
c         if (.not.(some conditions on jtot, iblock, ll, iqn1 (etc))) cycle

c  now counting only basis functions included in the current symmetry block
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

C the quantities below are only utilised if quadrature is to be used
      INTEGER, INTENT(INOUT)          :: NDIM, NPTS(NDIM), IXFAC, MX

      DOUBLE PRECISION, INTENT(INOUT) :: XFN(*)

      DOUBLE PRECISION, INTENT(OUT)   :: XPT(MXPT,NDIM), XWT(MXPT,NDIM)

      INTEGER, INTENT(IN)             :: MXPT, IVMIN, IVMAX, L1MAX,
     1                                   L2MAX, MXLMB

      double precision, allocatable :: fn1(:,:) !,fn2(:,:) etc
      integer :: il,npttot,nfun,ipt,ipt1tm,ipt1,i,ix,l1,iloop
      double precision pt1,wt1
      LOGICAL :: LVRTP
      DATA LVRTP/.FALSE./

      itypp = 1

      if (itypp.eq.9) then
        il=0
        do iloop=1,2
          do l1=0,l1max
c         do l2=0,l2max
            il=il+1
            if (iloop.eq.2) then
              lam(il)=l1
c             lam(il+mxlam)=l2
            endif
c         enddo
          enddo
          mxlam=il
        enddo
      endif

c  using quadrature is complicated.  Code should look something like this:
      if (lvrtp .and. itypp.eq.9) then

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
     1                                 NSTATE, JSTATE(NSTATE,NQNs),
     2                                 JSINDX(N), L(N), JTOT, IBOUND,
     3                                 IEXCH, IPRINT, LAM(*)

      integer :: irc,icol,irow,i1col,i1row,lcol,lrow,ipotl,idgvl
      logical livuse
      data livuse/.false./

c  Option to use IV array not supported.
c  This code included to silence compiler warning
      if (livuse) then
        iv(1,1) = 0
      endif

      DO ICOL = 1, N

c  if nconst = 0, there can be contributions to the internal energy which
c  may need to be recalculated during the course of the calculation.
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
          DO IPOTL = 1, NVLBLK

c  these will usually involve calculations that depend on the quantum labels
c  of the row and column and of the potential (or H_intl) expansion term
            if (ipotl.le.npotl) then
              VL(IPOTL, IRC) = 0.D0 ! the ipotl-th term in the potential
            elseif (ipotl.le.npotl+nconst) then
              vl(ipotl, irc) = 0.D0 ! the ipotl-npotl-th term in H_intl
            elseif (ipotl.le.npotl+nconst+nrsq) then
              vl(ipotl, irc) = 0.D0 ! the ipotl-npotl-nconst-th term in L^2
            else
              vl(ipotl, irc) = 0.D0 ! the ipotl-npotl-nconst-nrsq-th
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
