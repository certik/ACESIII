
c This routine retrieves a record from the JOBARC file.

c INPUT
c int IFLAG : a behavior flag
c             > 0; extant record is retrieved or process is killed
c             = 0; IRECLEN is modified as output
c                  >=  0; extant record integer-length
c                   = -1; record does not exist
c             < 0; extant record is retrieved or IDEST(1:IRECLEN) is zeroed out
c char*(*) SZARCHIVE : the internal filename of the record archive
c                      NOTE: Currently this is unused because all the archive
c                            statistics relate to JOBARC alone; however, we
c                            could expand the functionality to name any
c                            arbitrary archive.
c char*(*) SZRECNAME : the name of the record to retrieve

c INPUT/OUTPUT
c int IRECLEN : on input,  this is the integer-length of the record to retrieve
c               on output, if IFLAG=0, the int-length of the record is returned

c OUTPUT
c int IDEST : the destination array
c             NOTE: The record length is in units of integers regardless of
c                   the type of data actually stored in the record.

#include "jobarc.h" /* for _JA_RECWD */

      subroutine getrec(iFlag,szArchive,szRecName,iRecLen,iDest)
      implicit none

c ARGUMENTS
      integer iFlag, iRecLen, iDest(*)
      character*(*) szArchive, szRecName

c EXTERNAL FUNCTIONS
      integer iszeq

c INTERNAL VARIABLES
      character*80 szJOBARC
      integer       iJOBARC
      integer iRecNdx
      integer iRec, iBufNdx, iTmp
      integer nLeft, nGet, iOff
      integer iStat, iBuf(_JA_RECWD), i

c COMMON BLOCKS
#include "jobarc.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert job archive subsystem is up
      if (.not.bJAUp) then
         print *, '@GETREC: Assertion failed.'
         print *, '   bJAUp = ',bJAUp
         iTmp = 1
      end if
c   o assert record length is >= 0 if it's supposed to exist
      if ((iFlag.ne.0).and.(iRecLen.lt.0)) then
         print *, '@GETREC: Assertion failed.'
         print *, '   iRecLen = ',iRecLen
         iTmp = 1
      end if
c   o the record name must be between 1 and 8 characters
      if ((len(szRecName).lt.1).or.(len(szRecName).gt.8)) then
         print *, '@GETREC: Assertion failed.'
         print *, '   szRecName = ',szRecName
         iTmp = 1
      end if
c   o the record cannot be named 'OPENSLOT'
      if (len(szRecName).eq.8) then
      if (szRecName.eq.'OPENSLOT') then
         print *, '@GETREC: Assertion failed.'
         print *, '   szRecName = ',szRecName
         iTmp = 1
      end if
      end if
      if (iTmp.ne.0) then
         print *, '   record name = "',szRecName,'"'
         call aces_exit(iTmp)
      end if
#endif /* _ASSERT */

      if ((iFlag.ne.0).and.(iRecLen.lt.0)) return

c ----------------------------------------------------------------------

c   o get the index of the requested record
cYAU - This confused me at first (because I am stupid). Apparently,
c      iszeq declares marker as char*(*) and knows it is an array of
c      strings. Therefore, it steps through marker in 8-Byte increments,
c      which is the string length of one marker element.
      iRecNdx = iszeq(_MAX_JA_RECS,marker,1,szRecName)

c   o the record is not found
      if (iRecNdx.eq.0) then
         if (iFlag.lt.0) then
c         o return zeroes in iDest
            call izero(iDest,iRecLen)
         else
c         o the record has no length
            if (iFlag.ne.0) then
c            o complain and die
               print *, '@GETREC: The record "',szRecName,
     &                  '" is not found.'
               call aces_exit(1)
            end if
            iRecLen = -1
         end if
         return
      end if

c   o the record is found, but the caller doesn't care
      if (iFlag.eq.0) then
         iRecLen = rsize(iRecNdx)
         return
      end if

c   o the record is found, but the caller doesn't know the actual length
      if (iRecLen.ne.rsize(iRecNdx)) then
c      o only bomb for invalid data
         if (iRecLen.gt.rsize(iRecNdx)) then
            print *, '@GETREC: "',szRecName,'" record length mismatch'
            print *, '         requested length = ',iRecLen
            print *, '         stored    length = ',rsize(iRecNdx)
            call aces_exit(1)
         end if
      end if
      if (iRecLen.eq.0) return

c   o find the first physical record and integer index that point to
c     the first element
      iBufNdx = rloc(iRecNdx)
      iTmp    = (iBufNdx-1)/irecwd
      iRec    = 1       + iTmp
      iBufNdx = iBufNdx - iTmp*irecwd

#ifdef OPEN_JOBARC_AS_NEEDED
c   o open the JOBARC file
      call gfname('JOBARC',szJOBARC,iJOBARC)
      open(unit=_JA_UNIT,file=szJOBARC(1:iJOBARC),
     &     form='UNFORMATTED',access='DIRECT',
     &     status='OLD',recl=irecln,err=666,iostat=iStat)
#endif

c   o read/copy the record
      read(unit=_JA_UNIT,rec=iRec,err=666,iostat=iStat) iBuf
      if (iRecLen.eq.1) then
         iDest(1) = iBuf(iBufNdx)
      else
         nLeft = iRecLen
         nGet  = min(nLeft,irecwd+1-iBufNdx)
         call icopy(nGet,iBuf(iBufNdx),1,iDest,1)
         iOff  = 1     + nGet
         nLeft = nLeft - nGet
         do while (nLeft.ne.0)
            nGet = min(nLeft,irecwd)
            iRec = iRec + 1
            read(unit=_JA_UNIT,rec=iRec,err=666,iostat=iStat) iBuf
            call icopy(nGet,iBuf,1,iDest(iOff),1)
            iOff  = iOff  + nGet
            nLeft = nLeft - nGet
         end do
      end if

#ifdef OPEN_JOBARC_AS_NEEDED
c   o close the JOBARC file
      close(unit=_JA_UNIT,status='KEEP',err=666,iostat=iStat)
#endif

      return

c   o JOBARC read error
 666  print *, '@GETREC: read error on JOBARC'
      print *, '         record name = "',szRecName,'"'
      print '(/)'
      call aces_io_error('GETREC',_JA_UNIT,iStat)

c     end subroutine getrec
      end

