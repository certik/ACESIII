
c This routine finalizes the ACES environment by writing out information which
c allows later member executables to run successfully. This must be the last
c ACES routine called by the executable.

      subroutine aces_fin
      implicit none

c TIMING VARIABLES
      integer stime_sec, stime_usec
      integer utime_sec, utime_usec
      integer rtime_sec, rtime_usec
      double precision         ss, st,         us, ut,         rs, rt
      integer          sh, sm,         uh, um,         rh, rm

c COMMON BLOCKS
#include "aces_time.com" /* for timing stats */
#include "parallel_aces.com" /* for parallel environment */

c ----------------------------------------------------------------------

c   o finalize the I/O subsystem
      call aces_io_fin

c   o finalize the job archive subsystem
      call aces_ja_fin

c   o print the total execution time since aces_init
      call c_rutimes(utime_sec,utime_usec,stime_sec,stime_usec)
      call c_gtod   (rtime_sec,rtime_usec)
      ame_stime_out = (1.d-6 * stime_usec) + stime_sec
      ame_utime_out = (1.d-6 * utime_usec) + utime_sec
      ame_rtime_out = (1.d-6 * rtime_usec) + rtime_sec
      if (ame_timed) then
         st = ame_stime_out - ame_stime_in
         ut = ame_utime_out - ame_utime_in
         rt = ame_rtime_out - ame_rtime_in
         call ds2hms(st,sh,sm,ss)
         call ds2hms(ut,uh,um,us)
         call ds2hms(rt,rh,rm,rs)
         if (nprocs.eq.1) then
            print '(/)'
            print '(3a)', "ACES ---SYSTEM--- (total sec)",
     &                        " ----USER---- (total sec)",
     &                         " -WALLCLOCK- (total sec)"
 10         format (a,2(1x,i3,":",i2,":",f5.2," (",f9.1,")"),
     &                 (1x,i3,":",i2,":",f4.1," (",f9.1,")") )
            print 10, "TIME",sh,sm,ss,st,uh,um,us,ut,rh,rm,rs,rt
         else
            if (irank.eq.0) then
 11            format (a,(1x,i3,":",i2,":",f4.1," (",f9.1,")") )
               print 11,
     &       "Total wallclock time for the root process was",rh,rm,rs,rt
            end if
         end if
      end if

      return
c     end subroutine aces_fin
      end

