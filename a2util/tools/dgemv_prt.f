      subroutine dgemv_prt(trans,m,n,alpha,a,lda,x,incx,beta,y,incy,
     &                     label,labely,start,step)
      implicit none

      double precision alpha, beta
      integer incx, incy, lda, m, n
      character*1 trans
      double precision a(lda,*), x(*), y(*)
      character*(*) label, labely
      integer start, step

      character*(80) header
      integer i, info
c
      info = 0
      if (.not.((trans.eq."N").or.(trans.eq."n").or.
     &          (trans.eq."T").or.(trans.eq."t").or.
     &          (trans.eq."C").or.(trans.eq."c")    )) then
         info = 1
      else if (m.lt.0) then
         info = 2
      else if (n.lt.0) then
         info = 3
      else if (lda.lt.max(1,m)) then
         info = 6
      else if (incx.eq.0) then
         info = 8
      else if (incy.eq.0) then
         info = 11
      end if
      if (info.ne.0) then
         call xerbla("DGEMV_PRT",info)
         return
      end if

      if ((m.eq.0).or.(n.eq.0).or.
     $    ((alpha.eq.(0.0D0)).and.(beta.eq.(1.0D0)))) return

c This was in the original code. I wonder why?
      if (alpha.eq.(0.0D0)) return

      if ((trans.eq."N").or.(trans.eq."n")) then
         do i = 0,(m-1)
            write(header,'(4A,I3,3A,I3,A,F10.6,A)')
     &         label," ",
     &         labely,"(",start+(i*step),") = ",
     &         labely,"(",start+(i*step),") [* ",beta,"] +"
            call ddot_prt(n,a((1+i),1),lda,x,incx,header,alpha)
         end do
      else
         do i = 0,(n-1)
            write(header,'(4A,I3,3A,I3,A,F10.6,A)')
     &         label," ",
     &         labely,"(",start+(i*step),") = ",
     &         labely,"(",start+(i*step),") [* ",beta,"] +"
            call ddot_prt(m,a(1,(1+(i*lda))),1,x,incx,header,alpha)
         end do
      end if
      return
      end
