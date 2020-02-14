recursive subroutine multop(i_op,n_ops,eval,evec,degtol,nn,ncurr,ibeg,iprint)
!  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
!  Distributed under the GNU General Public License, version 3
USE potential, only: nconst,nrsq
! This subroutine diagonalises the (long-range) Hamiltonian matrix in order to
! find the threshold energies and eigenstates.  It takes each operator in turn
! (which may or may not be included in the Hamiltonian itself) and diagonalises
! the matrix over the basis functions.  It then transforms the subsequent
! operator matrices into that basis.  For each set of degenerate eigenvalues,
! it then selects that sub-block for the next operator and diagonalises for that
! next operator within that sub-block.  It uses recursion to achieve this.
!
! written by CRLS 12-11-13
!
! altered by CRLS 10-03-15 to include operators that are not part of the
! Hamiltonian.
!
! altered by CRLS 24-03-15 to reorder eigenvectors into an order that
! matches the original Hamiltonian ordering as far as is possible.
!
! altered by CRLS 17-07-2018 to use new structure for extra operators
implicit none

integer,intent(in):: i_op,nn,ncurr,ibeg,iprint,n_ops
! on entry: i_op  is index of recursion
!           nn    is size of complete hamiltonian
!           ncurr is size of current hamiltonian sub-block
!           ibeg  is start index of current sub-block within whole hamiltonian
!           iprint controls output level

double precision,intent(inout):: eval(nn,n_ops), evec(nn,nn,n_ops)
! on entry: evec contains W matrices for all hamiltonian parts, in basis
!                appropriate to H(i_op-1)
! on exit:  eval(:,i_op) contains eigenvalues of H(i_op)
!           evec(:,i_op) contains eigenvectors of H(i_op) in basis
!                          appropriate to H(i_op-1)
! on final exit: eval      contains eigenvalues for complete Hamiltonian
!                evec(:,1) contains eigenvectors for complete Hamiltonian,
!                          in original basis
double precision,intent(in)::degtol
! on entry: degtol is a constant against which the criterion of degeneracy is
! tested

! internal variables
integer ifail,info,jbeg,jend,nnext,nsets,nopinH
integer i,j,ij,iset,j_op,ival,icol,irow,jval,ieval,jvec,jpoint,nbefore
integer icount,iend,ivec
logical zinc,zorder
double precision query(1),avgval,maxel

integer, allocatable::setbeg(:),setend(:),point(:)
double precision, allocatable::tot_evec(:,:),sum_eval(:), &
                               wks(:,:),etmp(:)
logical, allocatable::zused(:)

iend=ibeg+ncurr-1

! set up the pointer array for rearranging the ordering of eigenvectors
if (i_op .eq. 1) call evbord(nn,n_ops)

if (iprint.ge.10) then
  if (i_op .eq. 1) write(6,899) 'W matrix has size',nn
  write(6,900) 'On this call to multop, operator number (i_op) =',i_op, &
               ', IBEG = ',ibeg,', and NCURR = ',ncurr
endif
899 format(/2x,a,1x,i6)
900 format(/2x,a,i2,a,i2,a,i2)

if (iprint.ge.25 .and. i_op.gt.1) then
  call matprn(6,evec(ibeg,ibeg,i_op),nn,ncurr,ncurr, &
              2,evec(ibeg,ibeg,i_op),' Submatrix of this operator',1)
endif

! diagonalise current W matrix
! on return evec contains current eigenvectors
allocate (wks(nn,ncurr))
wks(1:ncurr,1:ncurr)=evec(ibeg:iend,ibeg:iend,i_op)
call diagvc(wks,nn,ncurr,eval(ibeg,i_op),evec(ibeg,ibeg,i_op))
deallocate (wks)

! check whether current operator is in Hamiltonian in order to know whether
! reordering might be necessary
nopinH=min(1,nconst)+min(1,nrsq)
zorder=i_op.gt.nopinH

! if current operator is not in Hamiltonian, construct pointer array for
! eigenvectors/values to match ordering within Hamiltonian
! (but only apply it just before returning control to ytrans - see call to
! eveord below)
if (zorder .and. ncurr.gt.1) then
  call ev_ord(evec,ibeg,ncurr,nn,n_ops,i_op,iprint)
endif

if (iprint.ge.15) then
  write(6,*) ' Eigenvalues of this submatrix are'
  write(6,905)(eval(ibeg+ival,i_op),ival=0,ncurr-1)
905 format(1x,(5(1x,1pg22.15)))
  if (iprint.ge.25) call matprn(6,evec(ibeg,ibeg,i_op),nn,ncurr,ncurr, &
                                3,evec(ibeg,ibeg,i_op),' with eigenvectors',1)
endif

if (i_op .eq. n_ops) return ! end of recursion

! transform subsequent W matrices into this basis
if (iprint.ge.10) then
  write(6,*) ' Submatrices of subsequent operators ', &
             'are transformed into this basis'
endif

do j_op=i_op+1,n_ops
  evec(ibeg:iend,ibeg:iend,j_op)=                           &
       matmul(transpose(evec(ibeg:iend,ibeg:iend,i_op)),    &
                 matmul(evec(ibeg:iend,ibeg:iend,j_op),     &
                        evec(ibeg:iend,ibeg:iend,i_op)))
  if (iprint.ge.30) then
    write(6,940) j_op
940 format(/'  For operator #',i2,':')
    call matprn(6,evec(ibeg,ibeg,j_op),nn,ncurr,ncurr,2, &
                evec(ibeg,ibeg,j_op),' transformed submatrix',1)
  endif
enddo
if (iprint.ge.25) write(6,*) ' Finished transformation'

! search for sets of degenerate eigenvalues
allocate (setbeg(ncurr),setend(ncurr))
iset=1
setbeg(1)=1
setend(1)=1

do ival=2,ncurr
  jval=ibeg-1+ival
  if (abs(eval(jval,i_op)-eval(jval-1,i_op)).gt.degtol) then
    iset=iset+1
    setbeg(iset)=ival
  endif
  setend(iset)=ival
enddo

nsets=iset

! loop over degenerate sets of eigenvalues
do iset=1,nsets
  jbeg=setbeg(iset)+ibeg-1
  jend=setend(iset)+ibeg-1
  nnext=setend(iset)-setbeg(iset)+1
  if (iprint.ge.10 .and. jend.gt.jbeg) &
                   write(6,910) setbeg(iset),setend(iset), &
                   eval(jend,i_op)-eval(jbeg,i_op),degtol, &
                   (eval(ival,i_op),ival=jbeg,jend)
910 format(/'  Eigenvalues',i4,' to',i4,' differ by only',e12.5, &
           ' < DEGTOL =',e12.5,' (in internal units)'/9(1x,1pg22.15))

! set (almost) degenerate eigenvalues to be exactly degenerate
  if (nnext.gt.1) then
    if (iprint.ge.10) write(6,930) 'Called multop for operator', &
              i_op+1,' to construct a suitable linear combination ', &
              'of these eigenvectors'
930 format(/2x,a,i2,a,a)

! and call multop for the next operator with this subset...
    call multop(i_op+1,n_ops,eval,evec,degtol,nn,nnext,jbeg,iprint)
    do irow=1,nn
    do icol=jbeg,jend
      if (irow.ge.jbeg .and. irow.le.jend) cycle
      evec(irow,icol,i_op+1)=0d0
      evec(icol,irow,i_op+1)=0d0
    enddo
    enddo

! (but don't bother calling multop if the current sub-block has size 1)
  else
    do j_op=i_op+1,n_ops
      eval(jbeg,j_op)=evec(jbeg,jbeg,j_op)
      evec(jbeg,:,j_op)=0d0
      evec(:,jbeg,j_op)=0d0
      evec(jbeg,jbeg,j_op)=1d0
    enddo
  endif

enddo

deallocate (setbeg,setend)

if (iprint.ge.10 .and. i_op .ne. 1) then
  write(6,*)
  write(6,950) 'Returning to operator ',i_op-1
950 format(2x,a,i2)
endif
if (i_op .ne. 1) return
if (nconst .eq. 0) return

! the following is only executed just before return to YTRANS

allocate (tot_evec(nn,nn),sum_eval(nn))
if (n_ops.le.2) sum_eval=eval(:,1)

! multiply all the eigenvectors together...
tot_evec=evec(:,:,n_ops)
do j_op=n_ops-1,1,-1
  tot_evec=matmul(evec(:,:,j_op),tot_evec)
  if (j_op.ne.2) cycle
  do j=1,nn
    sum_eval(j)=dot_product(eval(:,1),tot_evec(j,:)**2)
  enddo
enddo

! replace the eigenvalues of the Hamiltonian with expectation values
eval(:,1)=sum_eval


if (iprint.ge.15) then
  write(6,*)
  write(6,*) ' Complete EIGENVALUES before reordering:'
  write(6,45) sum_eval(:)
  if (iprint.ge.25) call matprn(6,tot_evec,nn,nn,nn, &
                                3,tot_evec,' Complete EIGENVECTORS:',1)
endif
deallocate (tot_evec,sum_eval)

! then rearrange the ordering of the eigenvalues and eigenvectors
if (iprint.gt.24) &
  write(6,*) ' Reordering current set of eigenvalues and vectors'

call eveord(nn,evec,eval,n_ops)

if (iprint.ge.15) then
  write(6,*)
  write(6,*) ' Complete EIGENVALUES after reordering:'
  write(6,45) eval(:,1)
  if (iprint.ge.25) call matprn(6,evec,nn,nn,nn, &
                                3,evec,' Complete EIGENVECTORS:',1)
endif

45 format((3X,7(F22.15)))

return
end
