SUBROUTINE YTRANS(Y,EVEC,EINT,WVEC,                      &
                  JSINDX,L,N,P,VL,IV,                    &
                  MXLAM,NHAM,ERED,EP2RU,CM2RU,DEGTOL,    &
                  NOPEN,IBOUND,CENT,IPRINT,LQUIET)
!  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
!  Distributed under the GNU General Public License, version 3
USE potential
!  Routine to transform the log derivative matrix (Y) with primitive
!  basis into that with asymptotic basis.
!
!  ON ENTRY,  Y HOLDS THE LOG DERIVATIVE MATRIX
!             (Represented in the primitive basis)
!
!  ON EXIT,   Y HOLDS THE TRANSFORMED LOG DERIVATIVE MATRIX
!             (Represented in the asymptotic basis)
!
!             EINT HOLDS THE THRESHOLD ENERIES
!             WVEC HOLDS THE WAVE VECTORS CORRESPONDING TO EINT
!
!             EVEC is used as temporary space, though (if necessary)
!             the L values corresponding to the eigenvalues of the
!             centrifugal operator are contained in the first column
!             of EVEC, for passing into YTOK.
!
!  Several different YTRANS amalgamated by CRLS 24-04-15
!
!  altered by CRLS 17-07-2018 to use new structure for extra operators
!
!  This latest version has the capability to find degenerate sets of
!  eigenvalues and then to form linear (orthonormal) combinations of
!  them based on diagonalisation of futher operators (which may or may
!  not be included in the Hamiltonian).
!  It assumes that if NRSQ/=0, the terms that would be multiplied by 1/r^2
!  are the centrifugal operator terms, and so extracts eigenvalues of this
!  operator.  The resulting eigenvalues of the L^2 operator are stored
!  in CENT, ready for passing into YTOK.
!
!  This version is capable of dealing with primitive bases that are
!  body-fixed (see, for example T.V. Tscherbul and A. Dalgarno, JCP,
!  133 (2010) 184104).
!
implicit none

integer, intent(in)            ::N,IV(1),JSINDX(N),IBOUND,IPRINT,MXLAM,NHAM
integer, intent(inout)         ::L(N)
integer, intent(out)           ::NOPEN
double precision, intent(in)   ::VL(1),DEGTOL,EP2RU,CM2RU,ERED
double precision, intent(inout)::Y(N,N),CENT(N)
double precision, intent(out)  ::EVEC(N,N),EINT(N),WVEC(N)
logical                          LQUIET
!  array used as workspace
double precision P(1)

!-----------------------------------------------------------

integer          :: MX,IXNEXT,NIPR,IDUMMY
double precision :: X
COMMON /MEMORY/ MX,IXNEXT,NIPR,IDUMMY,X(1)

character(1)     :: CDRIVE
COMMON /CNTROL/ CDRIVE

integer          :: IVLFL, IVLU
COMMON /VLFLAG/ IVLFL
COMMON /VLSAVE/ IVLU

!  internal variables
double precision, allocatable  :: Wsub(:,:,:),eigmag(:)
double precision, allocatable  :: eval(:,:),wks(:,:,:),r_L(:)
integer, allocatable           :: nchan(:),mconst(:),maskCorL(:)
integer                           jprint,Lmin,Lmax,LL,nn,icol,irow,i,j,imax, &
                                  it,i_op,n_ops,irsq,i_dim,iextra,ibeg,      &
                                  NEXTRA_local,index_CENT
integer, external              :: IDAMAX
double precision                  realL,tol_L,wmax,CENTcur,CENTmax
logical                           no37,nowarn,zdegen
character(6)                   :: l_one
character(19)                  :: string
data no37,nowarn/.True.,.False./
data tol_L/1d-10/

jprint=IPRINT
if (LQUIET) jprint=0

if (NCONST.eq.0 .and. NRSQ.eq.0 .and. NEXTRA.eq.0) then
  return
endif

NEXTRA_local=NEXTRA
if (cdrive.ne.'M') then
  NEXTRA_local=0 ! can ignore extra operators for BOUND/FIELD
  string='thresholds'
else
  string='asymptotic channels'
endif

n_ops=min(1,NCONST)+min(1,NRSQ)+NEXTRA_local
allocate (mconst(n_ops))
i=0
if (NCONST.ge.1) then
  i=1
  mconst(1)=MXLAM+NCONST
endif
if (NRSQ.ge.1) then
  i=i+1
  mconst(i)=NRSQ
endif
do iextra=1,NEXTRA_local
  mconst(iextra+i)=NEXTMS(iextra)
enddo

if (IVLU.ne.0 .or. IVLFL.gt.0) then
  write(6,*) "I haven't coded this case: more work needed"
  stop
endif

if ((jprint.ge.6 .and. n_ops.gt.1 .and. CDRIVE.eq.'M') .or. jprint.ge.10) then
  write(6,890) string
endif
890  format(/'  Coefficients of operators to be diagonalised to find ',a)

allocate(wks(N,N,n_ops))
ibeg=1
do i_op=1,n_ops
  i=ibeg
  i_dim=mconst(i_op)
  if (nconst.gt.0 .and. i_op.eq.1) then
    P(1:mxlam)=0.d0
    P(mxlam+1:mxlam+nconst)=VCONST(1:nconst)*cm2ru/ep2ru
    call PERTRB(1.d30,P,NHAM,0)
    P(1:i_dim)=P(1:i_dim)*ep2ru
  else
    P(1:i_dim)=1.d0
  endif
  if ((jprint.ge.6 .and. n_ops.gt.1 .and. CDRIVE.eq.'M') .or. jprint.ge.10) then
    write(6,900) 'operator #',i_op,P(1:i_dim)/CM2RU
  endif
900  format(3x,a,i2,1p,10(g17.10,1x))
! next 7 lines are analogous to code in WAVVEC
  do j=1,N
    call dgemv('T',i_dim,j,1.0d0,VL(i),NVLBLK,P,1,0.0d0,wks(1,j,i_op),1)
    i=i+j*NVLBLK
  enddo
  ibeg=ibeg+mconst(i_op)
  call dsyfil('L',N,wks(1,1,i_op),N)
enddo

Lmin=0
Lmax=0

allocate(maskCorL(N))
if (NRSQ.eq.0 .and. IBOUND.eq.0) then
  Lmin=minval(L)
  Lmax=maxval(L)
  do i=1,N
    maskCorL(i)=L(i)
  enddo

elseif (NRSQ.eq.0 .and. IBOUND.ne.0) then
! in this case we want to block the asymptotic Hamiltonian by values in
! the CENT array rather than by values of L
  CENTcur=minval(CENT)
  CENTmax=maxval(CENT)
  maskCorL=0
  index_CENT=0
  do while (CENTcur.le.CENTmax)
    index_CENT=index_CENT+1
    do i=1,N
      if (maskCorL(i).ne.0) cycle
      if (abs(CENT(i)-CENTcur).le.DEGTOL) then
        maskCorL(i)=index_CENT
      endif
    enddo
    CENTcur=minval(CENT,maskCorL.eq.0)
  enddo
  Lmin=1
  Lmax=index_CENT
endif

if (.not.no37 .and. NRSQ.eq.0) then
  write(37) N,Lmin,Lmax
  write(37) JSINDX,L
endif

allocate (r_L(N))

EVEC=0.d0
do LL=Lmin,Lmax

  nn=N
  if (NRSQ.eq.0) then
    nn=0
    do i=1,N
      if (maskCorL(i).ne.LL) cycle
      if (IBOUND.ne.0) CENTcur=CENT(i)
      nn=nn+1
    enddo
  endif

  if (nn.eq.0) cycle

! can't be sure this will fit into any array passed into ytrans, so
! allocate it here
  allocate (Wsub(nn,nn,n_ops),eval(nn,n_ops))

! copy the currently relevant elements of the W matrices into Wsub
  icol=0
  do i=1,N
    if (NRSQ.eq.0 .and. maskCorL(i).ne.LL) cycle
    icol=icol+1
    irow=0
    do j=1,N
      if (NRSQ.eq.0 .and. maskCorL(j).ne.LL) cycle
      irow=irow+1
      do i_op=1,n_ops
        Wsub(irow,icol,i_op)=Wks(j,i,i_op)
      enddo
    enddo
  enddo

  if (jprint.ge.25) then
    if (NRSQ.eq.0 .and. IBOUND.eq.0) then
      write(6,'(/,a,I3)') '  Asymptotic Hamiltonian for L = ',LL
    elseif (NRSQ.eq.0 .and. IBOUND.ne.0) then
      write(6,'(/,a,G12.5)') '  Asymptotic Hamiltonian for CENT = ',CENTcur
    else
      write(6,'(/,a)') '  Asymptotic Hamiltonian'
    endif
    write(6,'(a,I3,a,I3)') '  submatrix',nn,'*',nn
    call MATPRN(6,Wsub,nn,nn,nn,2,Wsub, 'submatrix is:',1)
  endif

!  recursively (if necessary) diagonalise submatrices of asymptotic
!  hamiltonian/extra operators
  if ((jprint.ge.6 .and. n_ops.gt.1 .and. CDRIVE.eq.'M') .or. jprint.ge.10) then
    if (NRSQ.eq.0 .and. IBOUND.eq.0) then
      write(6,910) 'Calling multop for L =',LL
    elseif (NRSQ.eq.0 .and. IBOUND.ne.0) then
      write(6,911) 'Calling multop for CENT =',CENTcur
    else
      if (jprint.ge.10) write(6,*) ' Calling multop'
    endif
  endif
910   format(/2x,a,1x,i2)
911   format(/2x,a,1x,G12.5)
  call multop(1,n_ops,eval,Wsub,DEGTOL,nn,nn,1,jprint)

  if ((jprint.ge.6 .and. n_ops.gt.1 .and. CDRIVE.eq.'M') .or. jprint.ge.10) then
    write(6,*)
    if (NRSQ.eq.0 .and. IBOUND.eq.0) then
      write(6,910) 'Eigenvalues of all operators for L =',LL
    elseif (NRSQ.eq.0 .and. IBOUND.ne.0) then
      write(6,911) 'Eigenvalues of all operators for CENT =',CENTcur
    else
      if (NCONST.eq.0) then
        l_one='first'
      else
        l_one='second'
      endif
      write(6,912) ' Eigenvalues of all operators (',l_one,' one is L^2)'
912   format(a,a,a)
    endif
    i=0
    do icol=1,N
      if (NRSQ.eq.0 .and. maskCorL(icol).ne.LL) cycle
      i=i+1
      write(6,920) i,icol,(eval(i,i_op),i_op=1,n_ops)
    enddo
    write(6,*)
  endif
920  format(i5,1x,i5,1x,10F17.10)

  if (NRSQ.gt.0) then
! copy eigenvalues of L^2 operator into r_L
    if (NCONST.gt.0) then
      irsq=2
    else
      irsq=1
    endif
    r_L(:)=eval(:,irsq)
  endif


!  Choose sign of each eigenvector to make the element of largest
!  magnitude positive. The intent is to smooth out discontinuities
!  in off-diagonal elements of the transformed matrix.
  do i = 1, nn
!  IDAMAX is BLAS function
!  this function returns (minimum) index of the element
!  having maximum absolute value. Returned index is from 1 to nn.
    imax = IDAMAX(nn,Wsub(:,i,1),1)

!  When the largest magnitude of amplitude is negative...
    if (Wsub(imax,i,1).lt.0.0D0) then

!  DSCAL is BLAS routine. Change the sign (multiply vector by -1).
      call DSCAL(nn,-1.0D0,Wsub(:,i,1),1)
    endif
  enddo

!  check for degeneracies
  if (.not.nowarn) then
    allocate(eigmag(n_ops))
    do i_op=1,n_ops
      eigmag(i_op)=maxval(eval(:,i_op))-minval(eval(:,i_op))
    enddo
    do j=1,nn-1
      do i=j+1,nn
        if (abs(eval(i,1)-eval(j,1)).lt.degtol*eigmag(1)) then
          zdegen=.true.
          do i_op=2,n_ops
            if (abs(eval(i,i_op)-eval(j,i_op)).gt. &
                degtol*eigmag(i_op)) zdegen=.false.
          enddo
          if (zdegen) then
            write(6,*)
            write(6,'(a,i5,a,i5,a)') '  Warning in YTRANS. Eigenvalues ', &
                                     j,' and ',i, &
                                     ' near-degenerate for all operators'
            do i_op=1,n_ops
              write(6,'(a,i3)') '  i_op =',i_op
              write(6,"(10x,'Eval(',i2,')',18x,'Eval(',i2,')',18x,'diff')") &
                    j,i
              write(6,*) eval(j,i_op),eval(i,i_op), &
                         eval(j,i_op)-eval(i,i_op)
            enddo
            write(6,*) ' Channels may be mixed'
          endif
        endif
      enddo
    enddo
    deallocate(eigmag)
  endif

!  Copy eigenvector submatrix into position in evec
  icol=0
  do i=1,N
    if (NRSQ.eq.0 .and. maskCorL(i).ne.LL) cycle
    icol=icol+1
    irow=0
    do j=1,N
      if (NRSQ.eq.0 .and. maskCorL(j).ne.LL) cycle
      irow=irow+1
      EVEC(j,i)=Wsub(irow,icol,1)
    enddo
    if (.not.no37 .and. NRSQ.eq.0) then
      write(37) LL,eval(i,1)
      write(37) (EVEC(j,i),j=1,N)
    endif
  enddo
  deallocate (Wsub)

!  Set up EINT array (threshold energies) from eigenvalues
  if (NCONST.gt.0) then
    icol=0
    do i=1,N
      if (NRSQ.eq.0 .and. maskCorL(i).ne.LL) cycle
      icol=icol+1
      EINT(i)=eval(icol,1)
    enddo
  endif

  deallocate(eval)
enddo

deallocate(wks)

!---------------------------------------------------------------
!Print sorted by the energy

if (jprint.ge.15 .and. NRSQ.gt.0 .and. NCONST.gt.0) then
! print energy-sorted list of thresholds, together with nearest integer
! value for L
  allocate(nchan(N))
  call chnsrt(nchan,EINT,N)

  write(6,*)
  write(6,*) ' ERED = ',ERED
  write(6,*) ' Index      Asymp channel no.     L       Energy'
  do i=1,N
    write(6,'(i5,i16,i15,E28.15)') &
               i,nchan(i),nint(sqrt(r_L(i)+0.25d0)-0.5d0),EINT(nchan(i))/CM2RU
  enddo
  write(6,*)
  deallocate(nchan)
endif

!  Count the number of open channels
call wvcalc(WVEC,wmax,ERED,EINT,NOPEN,N)

! Transform to asymptotic basis
if (lquiet) then
! next 5 lines to get round bug? in matmul which causes it to crash
! sometimes if large matrix multiplies are allocated back to same space
  allocate (eval(N,N),wks(N,N,2))
  eval=transpose(EVEC)
  wks(:,:,1)=matmul(eval,Y)
  Y=matmul(wks(:,:,1),evec)
! this is original code
! Y=matmul(matmul(transpose(EVEC),Y),EVEC)
  if (iprint.ge.20) then
    call MATPRN(6,EVEC,N,N,N,3,EVEC,' Eigenvectors (last time):',1)
    call MATPRN(6,Y,N,N,N,3,Y,' transformed Y:',1)
    do ll=1,NHAM
      it=ll
      do j=1,N
      do i=1,j
        wks(i,j,1)=vl(it)
        wks(j,i,1)=vl(it)
        it=it+NHAM
      enddo
      enddo
      write (6,*) ' For potential term ',ll,':'
      if (iprint.ge.23) call MATPRN(6,wks,N,N,N,2,wks,'original VL',1)
      call trnsfm(eval,wks(:,:,1),wks(:,:,2),N,.false.,.false.)
      call MATPRN(6,wks,N,N,N,2,wks,'transformed VL',1)
    enddo
  endif
  deallocate (eval,wks)
endif

! Place eigenvalues of L^2 operator in array CENT ready to be used
! by YTOK
if (NRSQ.ne.0) then
930 format(2x,a,i4,a,f8.4)
  do i=1,N
    realL=r_L(i)
    if (realL.lt.-tol_L) then
      write(6,*) ' Eigenvalue of centrifugal operator is not physical'
      write(6,930) 'Solution #',i,' is ',realL
      stop
    endif
    realL=sqrt(realL+0.25d0)-0.5d0
    if (IBOUND.eq.0) then
      if (abs(realL-dble(nint(realL))).gt.tol_L) then
        write(6,*) ' Eigenvalues of centrifugal operator are not integer'
        write(6,930) 'Solution #',i,' is ',realL
      endif
      L(i)=nint(realL)
    else
      if (abs(realL-dble(nint(realL))).le.tol_L) L(i)=nint(realL)
    endif
    CENT(i)=r_L(i)
  enddo
elseif (IBOUND.ne.0) then
  do i=1,N
    realL=CENT(i)
    if (realL.lt.-tol_L) then
      write(6,*) ' Eigenvalue of centrifugal operator is not physical'
      write(6,930) 'Solution #',i,' is ',realL
      stop
    endif
  enddo
else
  CENT(:)=dble(L*(L+1))
endif

deallocate (r_L)

return
end
