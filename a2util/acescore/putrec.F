
c This routine stores a record to the JOBARC file.

c INPUT
c int IFLAG : a behavior flag
c             == 0; nInts are flushed with zeroes (ignoring iSrc)
c             != 0; nInts are copied from iSrc
c char*(*) SZARCHIVE : the internal filename of the record archive
c                      NOTE: Currently this is unused because all the archive
c                            statistics relate to JOBARC alone; however, we
c                            could expand the functionality to name any
c                            arbitrary archive.
c char*(*) SZRECNAME : the name of the record to store
c int IRECLEN : the integer-length of the record to store
c int ISRC : the source array
c            NOTE: The record length is in units of integers regardless of
c                  the type of data actually stored in the record.

#include "jobarc.h" /* for _JA_RECWD */
c#define TRAP_PUTREC_LEAKS

      subroutine putrec(iFlag,szArchive,szRecName,iRecLen,iSrc)
      implicit none

c ARGUMENTS
      integer iFlag, iRecLen, iSrc(*)
      character*(*) szArchive, szRecName

c EXTERNAL FUNCTIONS
      integer iszeq

c INTERNAL VARIABLES
      character*80 szJOBARC
      integer       iJOBARC
      integer iRecNdx
      integer iRec, iBufNdx, iTmp
      integer iRec0, iRec1
      integer nLeft, nPut, iOff
      integer iStat, iBuf(_JA_RECWD), i

c COMMON BLOCKS
#include "jobarc.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert job archive subsystem is up
      if (.not.bJAUp) then
         print *, '@PUTREC: Assertion failed.'
         print *, '   bJAUp = ',bJAUp
         iTmp = 1
      end if
c   o assert record length is >= 0
      if (iRecLen.lt.0) then
         print *, '@PUTREC: Assertion failed.'
         print *, '   iRecLen = ',iRecLen
         iTmp = 1
      end if
c   o the record name must be between 1 and 8 characters
      if ((len(szRecName).lt.1).or.(len(szRecName).gt.8)) then
         print *, '@PUTREC: Assertion failed.'
         print *, '   szRecName = ',szRecName
         iTmp = 1
      end if
c   o the record cannot be named 'OPENSLOT' or '. LOST .'
      if (len(szRecName).eq.8) then
      if (szRecName.eq.'OPENSLOT'.or.szRecName.eq.'. LOST .') then
         print *, '@PUTREC: Assertion failed.'
         print *, '   szRecName = ',szRecName
         iTmp = 1
      end if
      end if
      if (iTmp.ne.0) then
         print *, '   record name = "',szRecName,'"'
         call aces_exit(iTmp)
      end if
#endif /* _ASSERT */

      if (iRecLen.lt.0) return

c ----------------------------------------------------------------------

c   o see if the record exists already
      iRecNdx = iszeq(_MAX_JA_RECS,marker,1,szRecName)

c   o accomodate unsafe programming that increases record lengths
      if (iRecNdx.gt.0) then
         if (iRecLen.gt.rsize(iRecNdx)) then
            print *, '@PUTREC: "',szRecName,'" record length ',
     &               'mismatch'
            print *, '         requested length = ',iRecLen
            print *, '         stored    length = ',rsize(iRecNdx)
#ifdef TRAP_PUTREC_LEAKS
            call aces_exit(1)
#endif
            print *, '!WARNING! LOSING DISK SPACE WHILE WRITING ',
     &               'LARGER RECORD'
            marker(iRecNdx) = '. LOST .'
            iRecNdx = 0
         end if
      end if

c   o It is possible to recover some lost disk space by first searching
c     for any ". LOST ." record that is larger than the new record. The
c     two concepts to keep in mind are 1) stepping through marker(:)
c     since multiple records could be lost and 2) not updating rsize().

      if (iRecNdx.eq.0) then
c      o the record is not found, so find an open slot
         iRecNdx = iszeq(_MAX_JA_RECS,marker,1,'OPENSLOT')
         if (iRecNdx.eq.0) then
c         o complain and die since there are no open slots left
            print *, '@PUTREC: maximum number of JOBARC records ',
     &               'exceeded while attempting to write'
            print *, '         record "',szRecName,'"'
            call aces_exit(1)
         end if
c      o update the statistics for this new record
         marker(iRecNdx) = szRecName
         if (iRecNdx.ne.1) rloc(iRecNdx) =   rloc(iRecNdx-1)
     &                                     + rsize(iRecNdx-1)
         if (2.le.iRecLen.and.iRecLen.le.irecwd) then
c         o make sure the logical record is on one physical record
            iRec0 = 1 + (rloc(iRecNdx)-1)/irecwd
            iTmp  = rloc(iRecNdx)+iRecLen-1
            iRec1 = 1 + (iTmp-1)/irecwd
            if (iRec0.ne.iRec1) then
               iTmp = 1 + (iRec1-1)*irecwd
c               print *, '@PUTREC: Fragmenting JOBARC'
c               print *, '         ',rloc(iRecNdx),
c     &                  ' -> ',iTmp,
c     &                  ' for ',iRecLen,' ints (wasting',
c     &                  1.d2*(iTmp-rloc(iRecNdx))/(iTmp+iRecLen-1),'%)'
               rloc(iRecNdx) = iTmp
            end if
         end if
         rsize(iRecNdx) = iRecLen
      end if
      if (iRecLen.eq.0) then
         bJAMod = .true.
         return
      end if

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

c   o read/copy/write the record
      if (iRec.le.nrecs) then
         read(unit=_JA_UNIT,rec=iRec,err=666,iostat=iStat) iBuf
      end if
      if (iRecLen.eq.1) then
         if (iFlag.ne.0) then
            iBuf(iBufNdx) = iSrc(1)
         else
            iBuf(iBufNdx) = 0
         end if
         write(unit=_JA_UNIT,rec=iRec,err=666,iostat=iStat) iBuf
         nrecs = max(nrecs,iRec)
      else
         nLeft = iRecLen
         nPut  = min(nLeft,irecwd+1-iBufNdx)
         if (iFlag.ne.0) then
            call icopy(nPut,iSrc,1,iBuf(iBufNdx),1)
         else
            do i = 0, nPut-1
               iBuf(iBufNdx+i) = 0
            end do
         end if
         write(unit=_JA_UNIT,rec=iRec,err=666,iostat=iStat) iBuf
         nrecs = max(nrecs,iRec)
         iOff  = 1     + nPut
         nLeft = nLeft - nPut
         do while (nLeft.ne.0)
            nPut = min(nLeft,irecwd)
            iRec = iRec + 1
            if (iRec.le.nrecs) then
               read(unit=_JA_UNIT,rec=iRec,err=666,iostat=iStat) iBuf
            end if
            if (iFlag.ne.0) then
               call icopy(nPut,iSrc(iOff),1,iBuf,1)
            else
               do i = 1, nPut
                  iBuf(i) = 0
               end do
            end if
            write(unit=_JA_UNIT,rec=iRec,err=666,iostat=iStat) iBuf
            nrecs = max(nrecs,iRec)
            iOff  = iOff  + nPut
            nLeft = nLeft - nPut
         end do
      end if

#ifdef OPEN_JOBARC_AS_NEEDED
c   o close the JOBARC file
      close(unit=_JA_UNIT,status='KEEP',err=666,iostat=iStat)
#endif

c   o mark JOBARC as modified
      bJAMod = .true.

      return

c   o JOBARC I/O error
 666  print *, '@PUTREC: I/O error on JOBARC'
      print *, '         record name = "',szRecName,'"'
      print '(/)'
      call aces_io_error('PUTREC',_JA_UNIT,iStat)

c     end subroutine putrec
      end

