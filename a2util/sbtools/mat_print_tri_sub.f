
      subroutine mat_print_tri_sub(which,diag,mat,nrow,title,collab,
     &                             rowlab,numwid,ndec,
     &                             r0_in,r1_in,c0_in,c1_in)
      implicit none

      integer nrow,numwid,ndec,r0_in,r1_in,c0_in,c1_in
      character*1 which,diag
      double precision mat(nrow,nrow)
      character*(*) title,collab,rowlab

      integer pagewid
      parameter (pagewid=80)

      character*80 f,fn,fb,frh,frb,fch,fc,fr
      character*10 tmpn,tmpm,tmpl,tmpw,tmpd,tmp
      integer rownumlen,numcol,col0,col1,n,irow,icol,cols(pagewid),off,
     &        nblank,c,r0,r1,c0,c1,cc,r,rr
      integer  strlen
      external strlen
      logical top

      c0=c0_in
      c1=c1_in
      r0=r0_in
      r1=r1_in

c Check some error conditions...
      if (nrow.lt.1 .or. ndec.lt.0 .or. numwid.lt.1 .or.
     &    numwid.gt.pagewid/2 .or. ndec.gt.numwid-2 .or.
     &    r0.lt.1 .or. r0.gt.nrow .or. c0.lt.1 .or. c0.gt.nrow .or.
     &    r0.gt.r1 .or. c0.gt.c1 .or. (which.ne.'u' .and. which.ne.'U'
     &    .and. which.ne.'l' .and. which.ne.'L') .or. (diag.ne.'d' .and.
     &    diag.ne.'D' .and. diag.ne.'n' .and. diag.ne.'N')) then
         write(*,*) '@MAT_PRINT_TRI_SUB: invalid parameter'
         stop
      end if

      top=(which.eq.'u' .or. which.eq.'U')
      off=1
      if (diag.eq.'d' .or. diag.eq.'D') off=0

c Some columns/rows end up totally blank.  Remove them.
      if (top) then
         if (c0.lt.r0+off) c0=r0+off
         if (r1.gt.c0-off) r1=c1-off
      else
         if (c1.gt.r1-off) c1=r1-off
         if (r0.lt.c0+off) r0=c0+off
      end if

c Some combinations don't include anything to print
      if ((top .and. r0.gt.c0+off) .or. (.not.top .and. r1.lt.c1-off))
     &   return

c The number of digits in the row
      rownumlen=1+(log10(nrow*1.d0))
c The maximum number of columns to print
      numcol=(pagewid-(strlen(rowlab)+1+rownumlen))/(numwid+1)

c Formats are:
c   fn = '(x,fW.D)'            : a single number
c   fb = '(W+1x)'              : a blank number field
c   frh= '(aN,x,iM,'           : a row header
c   frb= '(N+M+1x,'            : a blank row header
c   f  = '(aN,x,iM,L(x,fW.D))' : a left justified set of numbers
c   fch= '(aN+M+1,L(x,iW))'    : a column header
c   fc = '(N+M+1x,L(x,W('-')))': the 2nd line of the column header
c where
c   N : strlen(rowlab)
c   M : rownumlen
c   L : numcol
c   W : numwid
c   D : ndec

      write(tmpw,'(i10)') numwid
      write(tmpd,'(i10)') ndec
      call delspc(tmpw,tmpw)
      call delspc(tmpd,tmpd)
      fn='(x,f'//tmpw(1:strlen(tmpw))//'.'//tmpd(1:strlen(tmpd))//')'
      write(tmp,'(i10)') numwid+1
      call delspc(tmp,tmp)
      fb='('//tmp(1:strlen(tmp))//'x)'

      write(tmpn,'(i10)') strlen(rowlab)
      write(tmpm,'(i10)') rownumlen
      call delspc(tmpn,tmpn)
      call delspc(tmpm,tmpm)
      frh='(a'//tmpn(1:strlen(tmpn))//',x,i'//tmpm(1:strlen(tmpm))//','
      write(tmp,'(i10)') strlen(rowlab)+rownumlen+1
      call delspc(tmp,tmp)
      frb='('//tmp(1:strlen(tmp))//'x,'

      write(tmpl,'(i10)') numcol
      call delspc(tmpl,tmpl)
      f=frh(1:strlen(frh))//tmpl(1:strlen(tmpl))//fn(1:strlen(fn))//')'

      fch='(a'//tmp(1:strlen(tmp))//','//tmpl(1:strlen(tmpl))//'(x,i'
     &    //tmpw(1:strlen(tmpw))//'))'

      col0=c0
      col1=c0+numcol-1
      if (col1.gt.c1) col1=c1
      write(*,*) title
      write(*,*)

   10 continue

c The actual number of columns to print
      n=col1-col0+1
      write(tmpn,'(i10)') n
      call delspc(tmpn,tmpn)
      fc=frb(1:strlen(frb))//tmpn(1:strlen(tmpn))//'(x,'//
     &   tmpw(1:strlen(tmpw))//'(''-'')))'

      do icol=col0,col1
         cols(icol-col0+1)=icol
      end do
      write(*,fch) collab,(cols(icol-col0+1),icol=col0,col1)
      write(*,fc)

      r=r0
      rr=r1
      if (top .and. rr.gt.col1-off) rr=col1-off
      if (.not.top .and. r.lt.col0+off) r=col0+off

      do irow=r,rr
         c=col0
         cc=col1
         nblank=0
         if (top .and. col0.lt.irow+off) then
            nblank=irow+off-col0
            c=irow+off
         endif
         if (.not.top .and. col1.gt.irow-off) cc=irow-off
         if (nblank.eq.0) then
            write(*,f) rowlab,irow,(mat(irow,icol),icol=c,cc)
         else
            write(tmp,'(i10)') nblank
            call delspc(tmp,tmp)
            fr=frh(1:strlen(frh))//tmp(1:strlen(tmp))//fb(1:strlen(fb))
     &         //','//tmpl(1:strlen(tmpl))//fn(1:strlen(fn))//')'
            write(*,fr) rowlab,irow,(mat(irow,icol),icol=c,cc)
         end if
      end do
      write(*,*)

      if (col1.lt.c1) then
         col0=col0+numcol
         col1=col1+numcol
         if (col1.gt.c1) col1=c1
         goto 10
      end if

      return
      end

