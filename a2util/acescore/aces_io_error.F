
c This routine is the error handler for I/O routines.

c INPUT:
c char*(*) SZROUTINE : the routine name from which aces_io_error is called
c int      IUNIT     : the external unit on which the I/O error occured
c int      ISTAT     : the I/O error number (from IOSTAT=iStat)

      subroutine aces_io_error(szRoutine,iUnit,iStat)
      implicit none
      character*(*) szRoutine
      integer iUnit, iStat
      print '(/)'
      print *, '@',szRoutine,': An ACES I/O Error has occurred.'
#ifdef _DEBUG
      print *, '   UNIT=',iUnit,'; IOSTAT=',iStat
#endif
      call c_strerror
      call aces_exit(iStat)
      end

