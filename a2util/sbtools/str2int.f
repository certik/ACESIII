
c This routine takes a string, which should contain a single valid integer,
c and returns the value of that integer. If an error occurs, err is set to:
c     1   if the string is empty
c     2   string contains a sign but no integer
c     3   integer part contains non-digit characters

c For those keen to regular expressions, valid integers are of the form:
c   '^ *[+-]? *[0-9]+ *$'

      subroutine str2int(str,val,err)
      implicit none
      character*(*) str
      integer val, err

      integer iachar, strlen, i, zref, ipm, j
      character*1 achar, spc, tab

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

c   o Parse the integer.
      val = 0
      err = 0
      do while (i.le.strlen)
         j = iachar(str(i:i))-zref
         if (j.lt.0.or.j.gt.9) err = 3
         val = 10*val + j
         i = i+1
      end do
      val = ipm*val

c      print *, '@STR2INT: "',str(1:strlen),'" = ',val

      return
      end

c      program unittest
c      implicit none
c      integer i, ierr
c      call str2int("   ",i,ierr)
c      if (ierr.ne.1) stop "ERROR: passes empty string"
c      call str2int(" + ",i,ierr)
c      if (ierr.ne.2) stop "ERROR: passes lone sign"
c      call str2int(" - ",i,ierr)
c      if (ierr.ne.2) stop "ERROR: passes lone sign"
c      call str2int(" 1 0 ",i,ierr)
c      if (ierr.ne.3) stop "ERROR: passes interstitial whitespace"
c      call str2int(" 10 ",i,ierr)
c      if (i.ne.10) stop "ERROR: fails to recognize ' 10 '"
c      call str2int("-10",i,ierr)
c      if (i.ne.-10) stop "ERROR: fails to recognize '-10'"
c      call str2int("+ 10",i,ierr)
c      if (i.ne.10) stop "ERROR: fails to recognize '+ 10'"
c      call str2int(" 10",i,ierr)
c      print *, 'STR2INT unit test: all tests pass'
c      end

