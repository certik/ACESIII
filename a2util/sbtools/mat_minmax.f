
c This routine takes one or two arrays of a given length, and
c returns a value depending on the values of op and neg.
c
c   op  neg   val
c    0    0   max( abs(a(i)-b(i)),           i=1,length )
c    0    1   max( a(i)-b(i),                i=1,length )
c    0   -1   max( b(i)-a(i),                i=1,length )
c    0    2   max( abs(abs(a(i))-abs(b(i))), i=1,length )
c    1    0   max( abs(a(i)),                i=1,length )
c    1    1   max( a(i),                     i=1,length )
c   -1    0   min( abs(a(i)),                i=1,length )
c   -1    1   min( a(i),                     i=1,length )
c
c When op!=0, the matrix b is ignored.

      subroutine mat_minmax(op,neg,val,length,a,b)
      implicit none

      integer op,neg,length
      double precision a(*),b(*),val

      integer i

      if (op.eq.0) then
         if (neg.eq.0) then
            val=abs(a(1)-b(1))
            do i=2,length
               val=max(val,abs(a(i)-b(i)))
            end do
         else if (neg.eq.1) then
            val=a(1)-b(1)
            do i=2,length
               val=max(val,a(i)-b(i))
            end do
         else if (neg.eq.-1) then
            val=b(1)-a(1)
            do i=2,length
               val=max(val,b(i)-a(i))
            end do
         else if (neg.eq.2) then
            val=abs(abs(a(1))-abs(b(1)))
            do i=2,length
               val=max(val,abs(abs(a(i))-abs(b(i))))
            end do
         else
            write(*,*) '@MAT_MINMAX: invalid value of neg'
            stop
         end if

      else if (op.eq.1) then
         if (neg.eq.0) then
            val=abs(a(1))
            do i=2,length
               val=max(val,abs(a(i)))
            end do
         else if (neg.eq.1) then
            val=a(1)
            do i=2,length
               val=max(val,a(i))
            end do
         else
            write(*,*) '@MAT_MINMAX: invalid value of neg'
            stop
         end if

      else if (op.eq.-1) then
         if (neg.eq.0) then
            val=abs(a(1))
            do i=2,length
               val=min(val,abs(a(i)))
            end do
         else if (neg.eq.1) then
            val=a(1)
            do i=2,length
               val=min(val,a(i))
            end do
         else
            write(*,*) '@MAT_MINMAX: invalid value of neg'
            stop
         end if

      else
         write(*,*) '@MAT_MINMAX: invalid value of op'
         stop
      end if

      return
      end

