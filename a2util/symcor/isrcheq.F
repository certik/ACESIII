
C EMULATOR OF CRAY SCILIB ROUTINE. LOCATES FIRST ELEMENT OF AN INTEGER
C VECTOR (IA) EQUAL TO A PARTICULAR VALUE (TARGET).

c The UNICOS man page for isrcheq is appended to this file.
c NOTE: Since Cray computers always use 64-bit integers, the vendor's
c       isrcheq accepts both integers and floats in the array argument.
c       This version only accepts integer arrays.

      integer function isrcheq(n,ia,inc,target)
      implicit none
      integer n, ia(*), inc, target, i, ndx
#ifdef _ASSERT
      if (n.lt.0) then
         print *, '@ISRCHEQ: Assertion failed.'
         print *, '          n = ',n
         call errex
      end if
#endif
      if (n.lt.1) then
         isrcheq=0
         return
      end if
      if (inc.ge.0) then
         if (inc.eq.1) then
            do i=1,n
               if (ia(i).eq.target) go to 2
            end do
         else
            ndx=1
            do i=1,n
               if (ia(ndx).eq.target) go to 2
               ndx=ndx+inc
            end do
         end if
      else
         if (inc.eq.-1) then
            ndx=n
            do i=1,n
               if (ia(ndx).eq.target) go to 2
               ndx=ndx-1
            end do
         else
            ndx=1-(n-1)*inc
            do i=1,n
               if (ia(ndx).eq.target) go to 2
               ndx=ndx+inc
            end do
         end if
      end if
 2    isrcheq=i
      return
      end

c NAME
c      ISRCHEQ, ISRCHNE - Searches a vector for the first element equal or
c      not equal to a target
c
c SYNOPSIS
c      index = ISRCHEQ(n,x,incx,target)
c      index = ISRCHNE(n,x,incx,target)
c
c DESCRIPTION
c      ISRCHEQ searches a real or integer vector for the first element that
c      is equal to a real or integer target.
c
c      ISRCHNE searches a real or integer vector for the first element that
c      is not equal to a real or integer target.
c
c      When using the Cray Fortran compiler UNICOS systems, all arguments
c      must be of default kind unless documented otherwise.  On UNICOS
c      systems, default kind is KIND=8 for integer, real, complex, and
c      logical arguments.
c
c      These functions have the following arguments:
c
c      index  Integer.  (output)
c             Index of the first element equal or not equal to target.  If
c             target is not found, n+1 is returned.  If n <= 0, 0 is
c             returned.
c
c      n      Integer.  (input)
c             Number of elements to be searched.
c
c      x      Real or integer array of dimension (n-1)* |incx| + 1.  (input)
c             Array x contains the vector to be searched.
c
c      incx   Integer.  (input)
c             Increment between elements of the searched array.
c
c      target Real or integer.  (input)
c             Value for which to search in the array.
c
c      The Fortran equivalent code for ISRCHEQ is as follows:
c
c                FUNCTION ISRCHEQ(N,X,INCX,TARGET)
c                INTEGER X(*), TARGET
c                J=1
c                ISRCHEQ=0
c                IF(N.LE.0) RETURN
c                IF(INCX.LT.0) J=1-(N-1)*INCX
c                DO 100 I=1,N
c                   IF(X(J).EQ.TARGET) GO TO 200
c                   J=J+INCX
c            100 CONTINUE
c            200 ISRCHEQ=I
c                RETURN
c                END
c
c      Although used as integers internally, you can use real values of x and
c      target, because ISRCHEQ and ISRCHEQ are matching bit patterns.
c
c NOTES
c      ISRCHEQ replaces the ISEARCH routine, but it has an entry point named
c      ISEARCH as well as ISRCHEQ.
c
c      When scanning backward (incx < 0), each routine starts at the end of
c      the vector and moves backward, as follows:
c
c      x(1-incx*(n-1)), x(1-incx*(n-2)), ..., x(1)
c
c      The desired value is at:
c
c      x(1+(index-1)*incx) when incx > 0
c      x(1+(index-n)*incx) when incx < 0

