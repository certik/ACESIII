
c This routine writes a physical record to a direct access file.

c INPUT
c int IUNIT   : the external unit number of the direct access file
c int IREC    : the index of the physical record to write
c int ISRC    : the (integer) array to write
c int ILENGTH : the integer-length of the record and source array

      subroutine aces_io_write(iUnit,iRec,iSrc,iLength)
      implicit none

c ARGUMENTS
      integer iUnit, iRec, iLength
      integer iSrc(iLength)

c INTERNAL VARIABLES
      integer iStat

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iStat = 0
c   o assert iUnit and iRec are natural
      if ((iUnit.lt.1).or.(iRec.lt.1)) then
         print *, '@ACES_IO_WRITE: Assertion failed.'
         print *, '   iUnit = ',iUnit
         print *, '   iRec  = ',iRec
         iStat = 1
      end if
c   o assert iLength is whole
      if (iLength.lt.0) then
         print *, '@ACES_IO_WRITE: Assertion failed.'
         print *, '   iLength = ',iLength
         iStat = 1
      end if
      if (iStat.ne.0) call aces_exit(iStat)
#endif /* _AASERT */

c ----------------------------------------------------------------------

      write(unit=iUnit,rec=iRec,err=666,iostat=iStat) iSrc
      return

 666  call aces_io_error('ACES_IO_WRITE',iUnit,iStat)

c     end subroutine aces_io_write
      end

