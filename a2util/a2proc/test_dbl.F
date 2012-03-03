
c This routine reads in doubles from iUnit and compares them to doubles from
c szRecName (stored in dScr). It returns the next record line in iUnit.

c INPUT
c int iUnit
c char*(*) szRecName
c int iDimScr

c OUTPUT
c char*80 sz
c logical bCont
c double dScr(iDimScr)

#include "iachar.h"

      subroutine test_dbl(iUnit,sz,bCont,szRecName,dScr,iDimScr)
      implicit none

c ARGUMENT LIST
      integer iUnit
      character*80 sz
      logical bCont
      character*8 szRecName
      double precision dScr(*)
      integer iDimScr

c INTERNAL VARIABLES
      double precision dVal, dTol, dTmp
      integer iTmp, nDbls, ndx
      integer iStat, i, iFirst, iLast, length
      logical bVal, bAbs
      character*80 szVal
      character*1 czSpace, czTab, czNull

c EXTERNAL FUNCTIONS
      double precision c_atof
      external c_atof

c COMMON BLOCKS
      integer            iErrExit
      common /exit_stat/ iErrExit
      save   /exit_stat/
#include "machsp.com"

c ----------------------------------------------------------------------

c   o get the record length in integers
      call getrec(0,'JOBARC',szRecName,iTmp,dScr)

c   o make sure the record exists as a double
      if (iTmp.lt.iintfp) then
         print *, '@TEST_DBL: "',szRecName,'" is empty'
         call aces_exit(1)
      end if

c   o read it in
      nDbls = iTmp/iintfp
      if (nDbls.gt.iDimScr) then
         print *, '@TEST_DBL: There is not enough memory to hold the ',
     &            'record.'
         call aces_exit(1)
      end if
      call getrec(1,'JOBARC',szRecName,iTmp,dScr)

c   o constants
      czSpace = achar(_IACHAR_SPACE)
      czTab   = achar(_IACHAR_TAB)
      czNull  = achar(0)

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
            length = iLast+1-iFirst
c         o check for absolute value
            if (length.gt.3) then
               bAbs = (sz(iFirst:iFirst+2).eq.'+/-')
               if (bAbs) then
                  iFirst = iFirst + 3
                  length = length - 3
               end if
            else
               bAbs = .false.
            end if
c         o normalize the exponential prefix to 'e'
            szVal = sz(iFirst:iLast)
            do i = 1, length
               if (szVal(i:i).eq.'d'.or.
     &             szVal(i:i).eq.'D'.or.
     &             szVal(i:i).eq.'E') szVal(i:i) = 'e'
            end do
            if (length.lt.80) then
               szVal(length+1:length+1) = achar(0)
            else
               szVal(80:80) = achar(0)
            end if
            dVal = c_atof(szVal)
c            print *,'DEBUG: "',szVal(1:length),'" ?= ',dVal
            if (dVal.eq.0.d0) then
c            o zeroes must be valid
               if (szVal(1:1).ne.'0'.and.szVal(1:1).ne.'.'.and.
     &             szVal(1:1).ne.'+'.and.szVal(1:1).ne.'-') bVal=.false.
               i = 2
               do while (szVal(i:i).ne.czNull.and.bVal.and.
     &                   szVal(i:i).ne.'d'.and.szVal(i:i).ne.'D'.and.
     &                   szVal(i:i).ne.'e'.and.szVal(i:i).ne.'E'
     &                  )
                  bVal = (szVal(1:1).eq.'0'.or.szVal(1:1).ne.'.')
                  i = i + 1
               end do
            end if
            if (bVal) then
               ndx = ndx + 1
               if (ndx.gt.nDbls) then
                  print *, '@TEST_DBL: There is less test data than ',
     &                     'good data'
                  call aces_exit(1)
               end if
c            o check for a second token and return the tolerance
               if (iLast.lt.80) then
                  call test_dbl_tol(szVal,dVal,sz(iLast+1:80),dTol)
               else
                  call test_dbl_tol(szVal,dVal,' ',dTol)
               end if
c            o check for agreement
               if (bAbs) dScr(ndx) = dabs(dScr(ndx))
               dTmp = dabs(dVal-dScr(ndx))
               if (dTmp.gt.dTol) then
                  print *, 'Double ',ndx,' of record "',szRecName,'"',
     &                     ' is not within tolerance.'
                  print '(a,e20.12,a,e11.2)',
     &                  '   good value = ',dVal,     '   +/-  ',dTol
                  print '(a,e20.12,a,e11.2)',
     &                  '   test value = ',dScr(ndx),' off by ',
     &                                     dScr(ndx)-dVal
                  iErrExit = iErrExit + 1
               end if
c           end if (bVal)
            end if
c        end if (iFirst.le.80)
         end if
c     end do while (bVal)
      end do

      if (ndx.ne.nDbls) then
         print *, '@TEST_DBL: WARNING: more test data than good data'
         if (nDbls.eq.1) then
            print *, '           "',szRecName,'" has 1 double'
         else
            print *, '           "',szRecName,'" has ',nDbls,' doubles'
         end if
         if (ndx.eq.1) then
            print *, '           Only 1 double was defined.'
         else
            print *, '           Only ',ndx,' doubles were defined.'
         end if
      end if

      return

 100  bCont = .false.
      return

 9999 call aces_io_error('TEST_DBL',iUnit,iStat)

c     end subroutine test_dbl
      end

c ----------------------------------------------------------------------

      subroutine test_dbl_tol(szVal,dVal,szTol,dTol)
      implicit none

c ARGUMENT LIST
      character*(*) szVal, szTol
      double precision dVal, dTol

c INTERNAL VARIABLES
      integer length, iFirst, iLast, i, iTmp
      character*81 szTmp
      character*1 czSpace, czTab, czNull

c EXTERNAL FUNCTIONS
      double precision c_atof
      external c_atof

c   o constants
      czSpace = achar(_IACHAR_SPACE)
      czTab   = achar(_IACHAR_TAB)
      czNull  = achar(0)
      length  = len(szTol)

c   o find the tolerance token
      iFirst = 1
      do while ((szTol(iFirst:iFirst).eq.czSpace.or.
     &           szTol(iFirst:iFirst).eq.czTab).and.
     &          iFirst.le.length)
         iFirst = iFirst + 1
      end do
      iLast = iFirst
      do while ((szTol(iLast:iLast).ne.czSpace.and.
     &           szTol(iLast:iLast).ne.czTab).and.
     &          iLast.lt.length)
         iLast = iLast + 1
      end do
c   o process the token
      if (iFirst.le.length) then
         if (szTol(iLast:iLast).eq.czSpace.or.
     &       szTol(iLast:iLast).eq.czTab) iLast = iLast - 1
c      o normalize the exponential prefix to 'e'
         szTmp = szTol(iFirst:iLast)//achar(0)
         do i = 1, iLast+1-iFirst
            if (szTmp(i:i).eq.'d'.or.
     &          szTmp(i:i).eq.'D'.or.
     &          szTmp(i:i).eq.'E') szTmp(i:i) = 'e'
         end do
         dTol = c_atof(szTmp)
         if (szTol(iLast:iLast).eq.'%') dTol = dVal*dTol*1.d-2
c     else if (no tolerance token)
      else
c      o base the tolerance on the rounding error of the good value
         i = 1
         do while (szVal(i:i).ne.'d'.and.szVal(i:i).ne.'D'.and.
     &             szVal(i:i).ne.'e'.and.szVal(i:i).ne.'E'.and.
     &             szVal(i:i).ne.czNull)
            iTmp = iachar(szVal(i:i))
            if (48.le.iTmp.and.iTmp.le.57) then
               szTmp(i:i) = '0'
            else
               szTmp(i:i) = szVal(i:i)
            end if
            i = i + 1
         end do
         szTmp(i:i) = '5'
         szTmp(i+1:81) = szVal(i:80)
         dTol = dabs(c_atof(szTmp))
c     end if (iFirst.le.length)
      end if
c      print '(a,e20.12,a,e20.12)', 'DEBUG: ',dVal,' +/- ',dTol

      return
c     end subroutine test_dbl_tol
      end

