
c The original ator (now a more complete str2dbl_old) forgave all kinds
c of spaces and tabs in the input string, but it did not allow fractional
c exponents.  str2dbl now ignores all whitespace and simply uses the
c internal unit numbers of Fortran.  It might be slower than str2dbl_old,
c but the parsing logic for all valid reals with embedded whitespace is
c just crazy.

      subroutine str2dbl(str,val,err)
      implicit none
      character*(*) str
      double precision val
      integer err
      integer i, j
      character*132 sz
      character*1 achar, spc, tab
      spc = achar(32)
      tab = achar(9)
      sz = ' '
      j = 1
      do i = 1, len(str)
         if (str(i:i).ne.spc.and.str(i:i).ne.tab.and.j.le.132) then
            sz(j:j) = str(i:i)
            j = j+1
         end if
      end do
      read(sz,fmt=*,err=10,iostat=err) val
 10   continue
c      print *, '@STR2DBL: "',str,'" = ',val
      return
      end

c This routine takes a string, which should contain a single valid float,
c and returns the value of that float. If an error occurs, err is set to:
c     1   if the string is empty
c     2   string contains a sign but no number
c     3   integer part contains non-digit characters (masked by err=5)
c     4   decimal part contains non-digit characters (masked by err=5)
c     5   exponent part doesn't start with 'E' or 'D'
c     6   exponent part empty
c     7   exponent part contains non-digit characters
c     8   a decimal point with no integer or decimal part entered

c For those keen to regular expressions, valid floats are of the form:
c   '^ *[+-]? *([0-9]+[.]?|[0-9]*\.[0-9]+) *([DdEe] *[0-9]+)? *$'

      subroutine str2dbl_old(str,val,err)
      implicit none
      character*(*) str
      double precision val
      integer err

      integer iachar, strlen, i, zref, ipm, j, k
      integer ii, id, il, ie
      character*1 achar, spc, tab, c

      zref = iachar('0')
      spc = achar(32)
      tab = achar(9)

c   o Skip trailing whitespace.
      i = len(str)
      do while ((str(i:i).eq.spc.or.str(i:i).eq.tab)
     &          .and.i.gt.0)
         i = i-1
      end do
      strlen = i

c   o Skip leading whitespace.
      i = 1
      do while ((str(i:i).eq.spc.or.str(i:i).eq.tab)
     &          .and.i.le.strlen)
         i = i+1
      end do
      if (i.gt.strlen) then
         err = 1
         return
      end if
      c = str(strlen:strlen)
      if (c.eq.'D'.or.c.eq.'d'.or.c.eq.'E'.or.c.eq.'e') then
         err = 6
         return
      end if

c   o Register sign and skip whitespace.
      ipm = 1
      if (str(i:i).eq.'-'.or.str(i:i).eq.'+') then
         if (str(i:i).eq.'-') ipm = -1
         i = i+1
         do while ((str(i:i).eq.spc.or.str(i:i).eq.tab)
     &             .and.i.le.strlen)
            i = i+1
         end do
         if (i.gt.strlen) then
            err = 2
            return
         end if
      end if
      c = str(i:i)
      if (c.eq.'D'.or.c.eq.'d'.or.c.eq.'E'.or.c.eq.'e') then
         err = 8
         return
      end if

c   o Read the integer.
      ii = 0
      j = i
      k = iachar(str(j:j))-zref
      do while (0.le.k.and.k.le.9.and.j.le.strlen)
         j = j+1
         k = iachar(str(j:j))-zref
      end do
      if (j.ne.i) call str2int(str(i:j-1),ii,err)
      i = j

c   o Read the decimal up to the exponent.
      id = 0
      il = 0
      if (str(i:i).eq.'.') then
         i = i+1
         j = i
         k = iachar(str(j:j))-zref
         do while (0.le.k.and.k.le.9.and.j.le.strlen)
            j = j+1
            k = iachar(str(j:j))-zref
         end do
         if (j.ne.i) then
            call str2int(str(i:j-1),id,err)
            il = j-i
         end if
         do while ((str(j:j).eq.spc.or.str(j:j).eq.tab)
     &             .and.j.le.strlen)
            j = j+1
         end do
         i = j
      end if

c   o Read the exponent.
      ie = 0
      if (str(i:i).eq.'D'.or.str(i:i).eq.'d'.or.
     &    str(i:i).eq.'E'.or.str(i:i).eq.'e'    ) then
         i = i+1
         if (i.le.strlen) then
            call str2int(str(i:strlen),ie,err)
            if (err.ne.0) then
               err = 7
               return
            end if
         end if
         i = strlen+1
      end if

      if (i.le.strlen) then
         err = 5
         return
      end if

      val = ii*1.d0 + id*10.d0**(-il)
      val = sign(val*10.d0**(ie),dble(ipm))

c      print *, '@STR2DBL: "',str(1:strlen),'" = ',val

      return
      end

