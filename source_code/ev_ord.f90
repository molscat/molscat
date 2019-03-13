subroutine ev_ord(evec,ibeg,ncurr,nn,nnctot,iconst,iprint)
!  Copyright (C) 2019 J. M. Hutson & C. R. Le Sueur
!  Distributed under the GNU General Public License, version 3
!
!  Written by CRLS 26-03-16
! This subroutine reorders eigenvectors and values for operators subsequent to
! those included in the Hamiltonian to match the original order.
implicit none

double precision,intent(inout)::evec(nn,nn,nnctot),eval(nn,nnctot)
! on entry: evec contains the eigenvectors
! on entry to eveord:  evec contains the eigenvectors
!                      eval contains the eigenvalues
! on exit from eveord: evec contains the reordered eigenvectors
!                      eval contains the reordered eigenvalues
integer,intent(in)::ibeg,ncurr,nn,nnctot,iconst,iprint
! on entry: ibeg is start index of current sub-block
!           ncurr is size of current sub-block
!           nn is size of complete basis
!           nnctot=nnconst+min(1,nrsq)
!           iconst is index of current operator
!           iprint controls output level

! local variables
integer iend,i,ivec,icmp,mcmp,mvec,posrep,j,jcmp,jmin
double precision cmpbig,cmpgt,summin
logical, allocatable::zmaskv(:),zmaskc(:)
double precision, allocatable::cmpmax(:),sumrep(:)
integer, allocatable::maxcmp(:,:),irep(:),jpoint(:,:)
save jpoint


! additional variables used by eveord entry point
integer ival,jconst,iv,jvec,jval,kconst
integer, allocatable::kvec(:),kval(:),ktemp(:)
double precision, allocatable::tot_evec(:,:),tr_evec(:,:),tr_eval(:)

iend=ibeg+ncurr-1
if (nnctot.eq.1) return
allocate (zmaskv(ibeg:iend),zmaskc(ibeg:iend),cmpmax(ibeg:iend), &
          maxcmp(ibeg:iend,ncurr), &
          irep(ibeg:iend),sumrep(ncurr))

! these govern which eigenvectors and which components are still left in the set
! to be searched through
zmaskv=.false.
zmaskc=.false.

! cycle round all the eigenvectors, eliminating them one-by-one
outer: do i=1,ncurr-1

! for the remaining eigenvectors, find the largest unused component
         cmpgt=0.d0
         do ivec=ibeg,iend
           if (zmaskv(ivec)) cycle
           cmpbig=0.d0
           do icmp=ibeg,iend
             if (zmaskc(icmp)) cycle
             if (abs(evec(icmp,ivec,iconst)).gt.cmpbig) then
               cmpbig=abs(evec(icmp,ivec,iconst))
               mcmp=icmp
             endif
           enddo
           cmpmax(ivec)=cmpbig
           maxcmp(ivec,1)=mcmp

! check to see whether the largest unused component is unique in that vector
           irep(ivec)=1
           do icmp=ibeg,iend
             if (zmaskc(icmp) .or. icmp.eq.mcmp) cycle
             if (abs(abs(evec(icmp,ivec,iconst))-cmpbig).lt.1.d-10) then
               irep(ivec)=irep(ivec)+1
               maxcmp(ivec,irep(ivec))=icmp
             endif
           enddo

! finally, find the eigenvector which has the largest unused component overall
           if (cmpbig.gt.cmpgt) then
             cmpgt=cmpbig
             mvec=ivec
           endif
         enddo

! looking at the largest unused component, check whether other vectors have
! their largest component in the same direction
         do j=1,irep(mvec)
           jcmp=maxcmp(mvec,j)
           posrep=0
           sumrep(j)=0.d0
           do ivec=ibeg,iend
             if (zmaskv(ivec) .or. ivec.eq.mvec) cycle
             do icmp=1,irep(ivec)
               if (maxcmp(ivec,icmp).eq.jcmp) then
                 posrep=posrep+1
                 sumrep(j)=sumrep(j)+cmpmax(ivec)
               endif
             enddo
           enddo

! if they don't, then choose this vector to be associated with this component
! position
           if (posrep.eq.0) then
             zmaskv(mvec)=.true.
             zmaskc(jcmp)=.true.
             jpoint(mvec,iconst)=jcmp
             cycle outer
           endif

! (check the same over all repeated greatest components)
         enddo

! if the largest unused component is unique within that vector, but other
! vectors have their largest component in the same place, assign this component
! position to this vector anyway
         if (irep(mvec).eq.1) then
           zmaskv(mvec)=.true.
           zmaskc(jcmp)=.true.
           jpoint(mvec,iconst)=jcmp
           write(6,*) ' *** Warning ***'
           write(6,100) '  Eigenvector ',mvec,' assigned to position',jcmp
100        format(a,i4,a,i4)
           write(6,*) ' but other eigenvectors also have their maximum',&
                      ' component in that position'
           write(6,*)
           cycle outer
         endif

! finally, if the largest unused component is not unique, and other vectors have
! their largest component in the same place(s) too, choose the position of the
! largest unused component that has the smallest overlap with other vectors
         jmin=0
         summin=1.d30
         do j=1,irep(mvec)
           jcmp=maxcmp(mvec,j)
           if (sumrep(j).lt.summin) then
             summin=sumrep(j)
             jmin=jcmp
           endif
         enddo

         zmaskv(mvec)=.true.
         zmaskc(jmin)=.true.
         jpoint(mvec,iconst)=jmin
         write(6,*) ' *** Warning ***'
         write(6,*) ' Eigenvectors are thoroughly mixed.  Somewhat arbitrary', &
                    ' choices have been made'
       enddo outer

! by now, just one vector and one component are not masked.
do ivec=ibeg,iend
  if (zmaskv(ivec)) cycle
  do icmp=ibeg,iend
    if (zmaskc(icmp)) cycle
    jpoint(ivec,iconst)=icmp
  enddo
enddo

deallocate (zmaskv,zmaskc,cmpmax,maxcmp,irep,sumrep)
if (iprint.ge.25) then
  write(6,110) 'For operator #',iconst
  write(6,111) 'Vector         ',(ivec,ivec=ibeg,iend)
  write(6,111) 'chosen ordering',(jpoint(ivec,iconst),ivec=ibeg,iend)
110 format(2x,a,i3)
111 format(2x,a,10(1x,i3))
endif

return

entry evbord(nn,nnctot)
! set up the pointer array used for reordering the eigenvectors and eigenvalues

if (nnctot.eq.1) return
allocate(jpoint(nn,nnctot))
do ival=1,nn
do jconst=1,nnctot
  jpoint(ival,jconst)=ival
enddo
enddo

return

entry eveord(nn,evec,eval,nnctot)
! reorder the eigenvectors and eigenvalues
! also multiply the eigenvectors together in reordered form to obtain the
! reordered total eigenvectors

!evec(:,ivec,2)->evec(:,jpoint(ivec,2),2)
! which means that evec(ival,ivec,3)->evec(jpoint(ival,2),
!                                          jpoint(jpoint(ivec,3),2),3)
! and
! evec(ival,ivec,4)->evec(jpoint(jpoint(ival,3),2),
!                         jpoint(jpoint(jpoint(ivec,4),3),2),4)
! etc

if (nnctot.eq.1) return
allocate (kvec(nn),kval(nn),ktemp(nn))

kval=jpoint(:,1)

allocate (tr_evec(nn,nn),tr_eval(nn),tot_evec(nn,nn))

tot_evec=evec(:,:,1)
do jconst=2,nnctot
  ktemp=jpoint(:,jconst)
  kvec=ktemp
  do kconst=jconst-1,2,-1
    kvec=jpoint(ktemp(:),kconst)
    ktemp=kvec
  enddo
  do ivec=1,nn
    jvec=kvec(ivec)
    do ival=1,nn
      jval=kval(ival)
      tr_evec(jval,jvec)=evec(ival,ivec,jconst)
    enddo
    tr_eval(jvec)=eval(ivec,jconst)
  enddo

  evec(:,:,jconst)=tr_evec
  eval(:,jconst)=tr_eval
  tot_evec=matmul(tot_evec,tr_evec)

  if (jconst.eq.nnctot) exit
  kval=kvec
enddo

evec(:,:,1)=tot_evec

deallocate(tot_evec,tr_evec,tr_eval,jpoint,kvec,kval,ktemp)

return
end
