
c This routine reads in integers from iUnit and compares them to ints from
c szRecName (stored in iScr). It returns the next record line in iUnit.

c INPUT
c int iUnit
c char*(*) szRecName
c int iDimScr

c OUTPUT
c char*80 sz
c logical bCont
c int iScr(iDimScr)

#include "iachar.h"

      subroutine test_int(iUnit,sz,bCont,szRecName,iScr,iDimScr)
      implicit none

c ARGUMENT LIST
      integer iUnit
      character*80 sz
      logical bCont
      character*8 szRecName
      integer iScr(*), iDimScr

c INTERNAL VARIABLES
      integer nInts, ndx, iFirst, iLast, iVal, i
      integer iStat
      logical bVal
      character*80 szVal
      character*1 czSpace, czTab

c EXTERNAL FUNCTIONS
      integer c_atol
      external c_atol

c COMMON BLOCKS
      integer            iErrExit
      common /exit_stat/ iErrExit
      save   /exit_stat/

c ----------------------------------------------------------------------

c   o get the record length
      call getrec(0,'JOBARC',szRecName,nInts,iScr)

c   o make sure the record exists
      if (nInts.lt.1) then
         print *, '@TEST_INT: "',szRecName,'" is empty'
         call aces_exit(1)
      end if

c   o read it in
      if (nInts.gt.iDimScr) then
         print *, '@TEST_INT: There is not enough memory to hold the ',
     &            'record.'
         call aces_exit(1)
      end if
      call getrec(1,'JOBARC',szRecName,nInts,iScr)

c   o constants
      czSpace = achar(_IACHAR_SPACE)
      czTab   = achar(_IACHAR_TAB)

c   o keep reading good data from iUnit until a new record or EOF
      ndx = 0
      bVal = .true.
      do while (bVal)
         read(unit=iUnit,fmt='(a)',end=100,err=9999,iostat=iStat) sz
c      o find the first token
         iFirst = 1
         do while ((sz(iFirst:iFirst).eq.czSpace.or.
     &              sz(iFirst:iFirst).eq.czTab).and.
     &             iFirst.le.80)
            iFirst = iFirst + 1
         end do
         iLast = iFirst
         do while ((sz(iLast:iLast).ne.czSpace.and.
     &              sz(iLast:iLast).ne.czTab).and.
     &             iLast.lt.80)
            iLast = iLast + 1
         end do
c      o process the token (skip blank lines)
         if (iFirst.le.80) then
            if (sz(iLast:iLast).eq.czSpace.or.
     &          sz(iLast:iLast).eq.czTab) iLast = iLast - 1
            szVal = sz(iFirst:iLast)//achar(0)
            iVal = c_atol(szVal)
c         o zeroes must be entered exactly as '0'
            bVal = (iVal.ne.0.or.sz(iFirst:iLast).eq.'0')
            if (bVal) then
               ndx = ndx + 1
               if (ndx.gt.nInts) then
                  print *, '@TEST_INT: There is less test data than ',
     &                     'good data'
                  call aces_exit(1)
               end if
               if (iVal.ne.iScr(ndx)) then
                  print *, 'Integer ',ndx,' of record "',szRecName,'"',
     &                     ' does not match.'
                  print *, '   good value = ',iVal
                  print *, '   test value = ',iScr(ndx)
                  iErrExit = iErrExit + 1
               end if
c           end if (bVal)
            end if
c        end if (iFirst.le.80)
         end if
c     end do while (bVal)
      end do

      if (ndx.ne.nInts) then
         print *, '@TEST_INT: WARNING: more test data than good data'
         if (nInts.eq.1) then
            print *, '           "',szRecName,'" has 1 integer'
         else
            print *, '           "',szRecName,'" has ',nInts,' integers'
         end if
         if (ndx.eq.1) then
            print *, '           Only 1 integer was defined.'
         else
            print *, '           Only ',ndx,' integers were defined.'
         end if
      end if

      return

 100  bCont = .false.
      return

 9999 call aces_io_error('TEST_INT',iUnit,iStat)

c     end subroutine test_int
      end

