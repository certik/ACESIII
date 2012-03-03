
      subroutine mat_print_sub(mat,nrow,ncol,title,collab,rowlab,
     &                         numwid,ndec,r0,r1,c0,c1)
      implicit none

      integer nrow,ncol,numwid,ndec,r0,r1,c0,c1
      double precision mat(nrow,ncol)
      character*(*) title,collab,rowlab

      integer pagewid
      parameter (pagewid=80)

      character*80 f,fn,fb,frh,frb,fch,fc
      character*10 tmpn,tmpm,tmpl,tmpw,tmpd,tmp
      integer rownumlen,numcol,col0,col1,n,irow,icol,cols(pagewid)
      integer  strlen
      external strlen

c Check some error conditions...
      if (nrow.lt.1 .or. ncol.lt.1 .or. ndec.lt.0 .or. numwid.lt.1
     &    .or. numwid.gt.pagewid/2 .or. ndec.gt.numwid-2 .or.
     &    r0.lt.1 .or. r0.gt.nrow .or. c0.lt.1 .or. c0.gt.ncol .or.
     &    r0.gt.r1 .or. c0.gt.c1) then
         write(*,*) '@MAT_PRINT_SUB: invalid parameter'
         stop
      end if

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

      fch='(a'//tmp(1:strlen(tmp))//','
     &        //tmpl(1:strlen(tmpl))//'(x,i'
     &        //tmpw(1:strlen(tmpw))//'))'

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

      do irow=r0,r1
         write(*,f) rowlab,irow,(mat(irow,icol),icol=col0,col1)
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

