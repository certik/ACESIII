C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
c---------------------------------------------------------------------------
c   Code to implement a timer subsystem.
c
c   Timers must first be registered by calling subroutine register_timer.
c   This routine accepts a character string describing the timer event, and
c   returns a key to use when referencing the timer.
c
c   As timer data is collected, it is updated by calls to subroutine
c   timer_start, followed by update_timer.  A call to timer_start starts
c   the timer, the corresponding call to update_timer turns it off and 
c   adds the result to the accumulator.
c
c   At the end of the program, the master process accumulates the timer data 
c   and prints a report of each timer that has been registered.
c---------------------------------------------------------------------------

      subroutine init_timers()
c---------------------------------------------------------------------------
c   Initialize the timer data structures.
c---------------------------------------------------------------------------

      implicit none
      include 'timerz.h'
      integer i

      do i = 1, max_timers
         tdesc(i) = ' '
         timers(i) = 0.
         tmark(i) = -1.
         timer_type(i) = 2
      enddo
      return
      end

      subroutine register_timer(desc, type, key)
c---------------------------------------------------------------------------
c   Registers a timer with the timer subsystem.  
c
c   Arguments:
c      desc (char)     	Character description of the timer being measured.
c      type             Type of timer.
c      key		Returned key to use in updating the timer.
c
c   This routine aborts if too many timers are registered.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'timerz.h'
      include 'saved_data.h'
      include 'parallel_info.h'

      character*(*) desc
      integer key, type
      integer i, j, desc_len, str_len
      integer ierr

      if (first_timer) then
         first_timer = .false.

         do i = 1, max_timers
            tdesc(i) = ' '
            timers(i) = 0.
            tmark(i) = -1.
            timer_type(i) = 2
         enddo
      endif

c------------------------------------------------------------------------
c   Validate timer type.
c------------------------------------------------------------------------

      if (type .ne. cpu_timer .and.
     *    type .ne. elapsed_time_timer) then
         print *,'Error: register_timer called with invalid timer type',
     *      ' for ',desc,'.  Type = ',type
         call abort_job()
      endif
c------------------------------------------------------------------------
c   Determine string length as either (1) len(desc) or (2) the
c   null-terminated string length (as in C).
c------------------------------------------------------------------------

      str_len  = max_timer_desc_len
      desc_len = str_len
      do j = 1, str_len
         if (ichar(desc(j:j)) .eq. 0) then
            desc_len = j-1
            go to 50
         endif
      enddo

   50 continue    
      if (desc_len .eq. 0) then
         print *,'Error: register_timer called with desc length 0'
         print *,'desc = ',desc,' type = ',type
         call abort_job()
      endif

c------------------------------------------------------------------------
c   Check to see if timer has already been registered.
c------------------------------------------------------------------------

      do i = 1, max_timers
         if (tdesc(i)(1:desc_len) .eq. desc(1:desc_len) .and.
     *       type .eq. timer_type(i)) then
            key = i   ! duplicate of a previously registered timer.
            return
         endif 
      enddo

      do i = 1, max_timers
         if (tdesc(i)(1:1) .eq. ' ') then
            tdesc(i)  = desc(1:desc_len)
            timers(i) = 0.
            tmark(i)  = -1.
	    timer_type(i) = type
	    key = i 
            return
         endif
      enddo

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      print *,'Error: Too many timers are already in service.'
      print *,'       There is a limit of ',max_timers,' timers.'
      print *,'       Attempted to register ',desc,' on processor ',
     *                   me
      call mpi_abort(mpi_comm_world, 1, ierr)
      return
      end

      subroutine timer_start(key)
c---------------------------------------------------------------------------
c   Starts the timer referenced by "key".
c---------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'timerz.h'
      include 'parallel_info.h'
      include 'trace.h' 
      integer key, ierr
      integer my_type
      integer utime_sec, stime_sec, utime_usec, stime_usec

      if (.not. do_timer) return

      ierr = 0
      if (key .eq. 0) return
      if (key .lt. 1 .or. key .gt. max_timers) then
         print *,'Error: timer_start was called with key ',key,
     *             ' out of range.'
         ierr = 1
      else
         if (tdesc(key) .eq. ' ') then
            print *,'Error: Timer ',key,' is being started without ',
     *              'being registered.'
            ierr = 2
         endif
      endif

      if (ierr .eq. 0 .and.
     *    tmark(key) .ne. -1.) then
         print *,'Task ',me,' Error: Timer ',tdesc(key),
     *             ' is being started ',
     *            'without a corresponding update. key = ',
     *            key
         ierr = 3
      endif

      if (ierr .eq. 0) then
         my_type = timer_type(key)
         if (my_type .eq. elapsed_time_timer) then
            tmark(key) = mpi_wtime()   ! turn on the timer
         else if (my_type .eq. cpu_timer) then
c            call c_rutimes(utime_sec, utime_usec, stime_sec, 
c     *                     stime_usec)
c            tmark(key) = 1.d-6*(utime_usec + stime_usec) + 
c     *                         (utime_sec + stime_sec)
         else 
            print *,'Error: Invalid timer type: key = ',key,
     *             ' type = ',timer_type(key)
            call abort_job()
         endif
      else
         print *,'Task ',me,' Aborting due to timer errors at line ',
     *        current_line
         call mpi_abort(mpi_comm_world, 1, ierr)   ! abort the job.
      endif

      return
      end

      subroutine update_timer(key)
c----------------------------------------------------------------------------
c   Accumulates "time_val" into the timer referenced by "key".
c----------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'timerz.h'
      include 'parallel_info.h'
      include 'trace.h'

      double precision time_val
      integer key
      integer ierr
      integer utime_sec, stime_sec, utime_usec, stime_usec

      if (.not. do_timer) return

      if (key .eq. 0) return

      if (key .lt. 1 .or. key .gt. max_timers) then
         call mpi_comm_rank(mpi_comm_world, me, ierr)
         print *,'Error: Process ',me,' attempted to update timer ',
     *           key,'.'
         call mpi_abort(mpi_comm_world, 1, ierr)
      else

c------------------------------------------------------------------------------
c  Valid timer key.  Determine if the timer has been started.
c-----------------------------------------------------------------------------

         if (timer_type(key) .eq. elapsed_time_timer) then
            time_val       = mpi_wtime() - tmark(key)
         else if (timer_type(key) .eq. cpu_timer) then
c             call c_rutimes(utime_sec,utime_usec,stime_sec,
c     *                      stime_usec)
c             time_val = 1.d-6*(utime_usec+stime_usec) + 
c     *                        (utime_sec+stime_sec) - tmark(key)
         endif

         timers(key)    = timers(key) + time_val
         tmark(key)  = -1.   ! turn off the timer.
      endif

      return
      end  

      subroutine stop_timers()
c-------------------------------------------------------------------------
c   Stops all timers on a process wihtout updating them.  This allows any
c   leftover timers that are still running to be turned off and reused at
c   a later point.
c-------------------------------------------------------------------------
      implicit none
      include 'timerz.h'

      integer i

      if (.not. do_timer) return

      do i = 1, max_timers
         tmark(i) = -1.
      enddo
      return
      end

      subroutine reset_timer_info()
c----------------------------------------------------------------------------
c   Resets all existing timers to 0.
c----------------------------------------------------------------------------
      implicit none
      include 'timerz.h'

      integer i

      if (.not. do_timer) return
      call stop_timers()    ! make sure they are not running.
      
      do i = 1, max_timers
         timers(i) = 0.
      enddo

      return
      end

      subroutine timer_report(siofile, descs, timer_data, contrib,
     *                        sumsq, nprocs)
c-------------------------------------------------------------------------
c   Called by the master to produce a report of the accumulated timer
c   data collected from each integral worker.
c-------------------------------------------------------------------------

      implicit none
      include 'timerz.h'

      integer ranks_per_line
      parameter (ranks_per_line = 4)

      integer max_calls
      parameter (max_calls = 1000)

      integer i, j, k, l, ndx, lndx, nprocs
      integer nline, nworker
      integer ntimers
      double precision timer_data(max_timers,3)
      double precision contrib(max_timers)
      double precision sumsq(max_timers)
      
      character*(max_timer_desc_len) descs(*)
      character*(*) siofile
      character*32 call_table(max_calls)
      character*256 aces_source_dir
      character*120 srcline, token
      character*256 srcfile, sialfile
      character*40  pardo_desc, blkwait_desc, timerdesc
      character*60 line_fmt

      integer n
      integer lineno, lenval, ierr
      integer str_trimlen
      integer nprepare, npreparesum, nrequest, nprequest
      integer ndesc, nxtcall, ncalls
      integer call_line(max_calls)
      double precision call_nproc(max_calls)
 
      integer ranks(ranks_per_line)
      double precision data(ranks_per_line)
      double precision sum, avg, sdev
      double precision temp
      double precision calc_sdev
      double precision call_sum(max_calls), call_sumsq(max_calls)
      real pardo_avg, blkwait_avg, efficiency
      real call_tmin(max_calls), call_tmax(max_calls)
      
      logical source_level_analysis 
      logical proc_timer, pardo_timer, line_timer, timer_exist

c---------------------------------------------------------------------------
c   Build a table of all the unique timer descriptions.
c---------------------------------------------------------------------------

      if (.not. do_timer) return   ! nothing to do.

c---------------------------------------------------------------------------
c   Check for the ACES_SOURCE_DIR env. variable.
c---------------------------------------------------------------------------

      aces_source_dir = ' '
      call c_getenv('ACES_SOURCE_DIR'//char(0),
     *            aces_source_dir, lenval, ierr)
      if (ierr .eq. 0) then
         source_level_analysis = .true.
         lenval = str_trimlen(aces_source_dir)
      else
         source_level_analysis = .false.
      endif

      if (source_level_analysis) go to 1000 

c---------------------------------------------------------------------------
c   Print individual worker's statistics.
c---------------------------------------------------------------------------

      print *,'---------- Summary of Timer Statistics ----------'
      print *,' '

c--------------------------------------------------------------------------
c   Calculate and print average and standard dev. of all timers.
c--------------------------------------------------------------------------

       print 400
       print 500
       call c_flush_stdout()
      do i = 1, max_timers
         avg  = timer_data(i,1) 
         sdev = 0.

         if (descs(i)(1:1) .ne. ' ' .and.
     *       descs(i)(1:1) .ne. char(0)) then
            if (contrib(i) .le. 1.) then
               avg = timer_data(i,1)
               sdev = 0.d0
            else  
               avg = timer_data(i,1)/contrib(i)
               sdev = calc_sdev(sumsq(i), timer_data(i,1), contrib(i))
            endif

            if (timer_data(i,3) .eq. 1.d10)    ! remove init value for min
     *          timer_data(i,3) = 0.d0
       
            if (avg .gt. 0.0005) 
     *         print 200,descs(i)(1:25), avg, sdev, 
     *             timer_data(i,3), timer_data(i,2) 
         endif
      enddo
      return

 1000 continue
 
c-----------------------------------------------------------------------
c   SOURCE LEVEL ANALYSIS: Find the name of the source file.
c-----------------------------------------------------------------------

      n = str_trimlen(siofile)
      if (siofile(n-3:n) .eq. '.sio') then
         srcfile = siofile(1:n-4) // '.sial'
      else
         print *,'Error: Cannot determine SIAL filename'
         print *,'siofile is ',siofile
         call abort_job() 
      endif

      n = str_trimlen(srcfile)
      print *,'------ Timer statistics for SIAL file ',srcfile(1:n),
     *        '------'

c-------------------------------------------------------------------------
c   Open the source file 
c-------------------------------------------------------------------------

      sialfile = aces_source_dir(1:lenval)// '/' // 
     *           srcfile(1:n) 
      open(unit=95, file=sialfile,
     *     form='formatted', status='OLD',err=1100,iostat=ierr)
     
      lineno  = 1
      nxtcall = 0
      go to 1200

 1100 continue
      print *,'Error opening src file: ',
     *     aces_source_dir(1:lenval)//srcfile(1:n) ,
     *     ' iostat = ',ierr
      call abort_job()

 1200 continue
      read(95,'(a)',end=2000) srcline

c------------------------------------------------------------------------
c   Omit any embedded nulls.
c------------------------------------------------------------------------

      nline = str_trimlen(srcline)
      do i = 1, nline
         if (srcline(i:i) .eq. char(0)) srcline(i:i) = ' '
      enddo

      line_fmt = ' '
      if (nline .ge. 100) then
         write (line_fmt, 9000) nline
      else if (nline .ge. 10 ) then
         write (line_fmt, 9050) nline
      else
         write (line_fmt, 9060) nline
      endif
 
c--------------------------------------------------------------------------
c   Do we have timer data to print for this line?
c--------------------------------------------------------------------------

      pardo_timer = .false.
      proc_timer = .false.

c-------------------------------------------------------------------------
c   Pick up the first token off the source line.
c-------------------------------------------------------------------------
       
      token = ' '
      call f_strtoken(srcline, ' ', token)

      if (token(1:4) .eq. 'call' .or.
     *         token(1:4) .eq. 'CALL') then
         nxtcall = nxtcall + 1
         token = ' '
         call f_strtoken(char(0), ' ', token)  ! get proc name 
         if (token(1:1) .eq. char(0)) then
            print *,'Error: Cannot find proc name for CALL at line ',
     *         lineno
            print *,'srcline = ',srcline
            print *,'token ',token
            call abort_job()
         else
            call_table(nxtcall) = token
            call_line(nxtcall)  = lineno
         endif
      else if (token(1:5) .eq. 'pardo' .or.
     *         token(1:5) .eq. 'PARDO') then
         pardo_timer = .true.

         pardo_desc = ' '
         write (pardo_desc,9300) lineno
         blkwait_desc = ' '
         write (blkwait_desc,9400) lineno
      else
         timerdesc = ' ' 
         write (timerdesc,9500) lineno
         ndesc = str_trimlen(timerdesc)
      endif 
      
      if (pardo_timer) then

c--------------------------------------------------------------------------
c   Calculate the average of all pardo and block wait timers for this loop.
c---------------------------------------------------------------------------

         pardo_avg = 0.
         blkwait_avg = 0.
         do i = 1, max_timers
            if (descs(i) .eq. pardo_desc) then
               if (contrib(i) .eq. 0.) then
                  pardo_avg = 0.
               else
                  pardo_avg = timer_data(i,1) / contrib(i)
               endif 
               go to 2200
            endif
         enddo

 2200    continue 
         do i = 1, max_timers
            if (descs(i) .eq. blkwait_desc) then
               if (contrib(i) .eq. 0.) then
                  blkwait_avg = 0.
               else
                  blkwait_avg = timer_data(i,1) / contrib(i)
               endif
               go to 2300 
            endif
         enddo

 2300    continue
         if (pardo_avg .gt. 0.) then
            efficiency = 100.0 * (pardo_avg - blkwait_avg) / pardo_avg
         else
            efficiency = 0.
         endif

c---------------------------------------------------------------------------
c   Now get the server-side data pertaining to the loop.
c---------------------------------------------------------------------------

         call get_server_stats(lineno, nprepare, npreparesum, nrequest,
     *                         nprequest ) 

         print 9200,pardo_avg, blkwait_avg, n, efficiency
         print *,'   Num. prepares ',nprepare,' Num. preparesum ',
     *              npreparesum,' Num. requests ', nrequest,
     *              ' Num. prequest ',nprequest
         print line_fmt, lineno,srcline(1:nline)
      else if (proc_timer) then
      else

c----------------------------------------------------------------------------
c   Check for a timer with this line number encoded into it.
c----------------------------------------------------------------------------

         timer_exist = .false.
      
         do i = 1, max_timers
            call f_strtoken(descs(i), ':', token)
            if (token(1:1) .ne. char(0)) then
               if (token(1:ndesc) .eq. timerdesc) then
                  if (contrib(i) .le. 1.) then
                     avg = timer_data(i,1)
                     sdev = 0.d0 
                  else
                     avg = timer_data(i,1) / contrib(i)
                     sdev = calc_sdev(sumsq(i), timer_data(i,1),
     *                                contrib(i))
                  endif
         
                  timer_exist = .true.
                  go to 1500
	       endif
            endif
         enddo

 1500    continue
         if (timer_exist) then
            line_fmt = ' '
            if (nline .ge. 100) then
               write(line_fmt, 9600) nline
            else if (nline .ge. 10) then
               write(line_fmt, 9650) nline
            else
               write(line_fmt, 9660) nline
            endif 

            print line_fmt, real(avg), real(timer_data(i,3)), 
     *             real(timer_data(i,2)), 
     *             real(sdev), lineno, srcline(1:nline) 
         else 
            print line_fmt, lineno, srcline(1:nline) ! no timer for this line
         endif 
      endif

      lineno = lineno + 1

      go to 1200
 
 2000 continue

c---------------------------------------------------------------------------
c   Make a pass through the proc calls to generate the call data.
c---------------------------------------------------------------------------

      ncalls = nxtcall
      do i = 1, ncalls
         timerdesc = ' '
         write (timerdesc, 9700) call_line(i)
         
         do j = 1, max_timers
            if (descs(j) .eq. timerdesc) then
               call_sum(i)   = timer_data(j,1)
               call_sumsq(i) = sumsq(j)
               call_tmin(i)  = timer_data(j,3)
               call_tmax(i)  = timer_data(j,2)
               call_nproc(i) = contrib(j)
               go to 2400
            endif
         enddo
 2400    continue
 
      enddo

c---------------------------------------------------------------------------
c   Print the call data.
c---------------------------------------------------------------------------

      print *,'---------- Table of proc calls -----------------------'
      print 9750

      do i = 1, ncalls
         if (call_nproc(i) .gt. 0) then
            avg = call_sum(i) / call_nproc(i)
            if (call_nproc(i) .le. 1.) then
               sdev = 0.d0
            else
               sdev = calc_sdev(call_sumsq(i), call_sum(i), 
     *                       call_nproc(i))
            endif
            n = str_trimlen(call_table(i)) 
	    print 9800,call_table(i)(1:n),call_line(i), real(avg), 
     *            call_tmin(i), call_tmax(i), real(sdev)
         endif
      enddo

c---------------------------------------------------------------------------
c   Make one last pass to print block wait time.
c---------------------------------------------------------------------------

      n = 0
      sum = 0.
      timerdesc = 'Block wait time'
      ndesc = str_trimlen(timerdesc)
      do i = 1, max_timers
         if (descs(i)(1:ndesc) .eq. timerdesc) then
            if (contrib(i) .le. 1.d0) then 
               avg = timer_data(i,1)
            else
               avg = timer_data(i,1) / contrib(i)
            endif
            go to 2500
         endif
      enddo
 2500 continue

      print *,'Average block wait time: ',avg  
      print *,'------------------------------------------------------'
      return

  100 format(1x,a25,5(1x,f12.3))
  200 format(1x,a25,1x,f12.3,1x,f12.3,7x,f8.3,3x,f8.3)
  300 format(1x,'Rank               ',9(9x,i4))
  400 format(32x,'Average    Standard deviation',
     *       '  Min. time   Max. time')
  500 format(25x,'-------------   ------------------',
     *       '   --------   --------')

      
 9000 format('(48x,''| '',i6,'': '',a',i3,')')
 9050 format('(48x,''| '',i6,'': '',a',i2,')')
 9060 format('(48x,''| '',i6,'': '',a',i1,')')
 9200 format('*** PARDO LOOP: time ',f12.5,' Block wait time ',f12.5,
     *       ' Num. procs ',i5,' Loop efficiency ',f7.3,'% ***')
 9300 format('Pardo at line ',i6)
 9400 format('Blkwait for pardo ',i6)
 9500 format('Line',i6)
 9600 format('(4(1x,f11.5),''| '',i6,'': '',a',i3,')')
 9650 format('(4(1x,f11.5),''| '',i6,'': '',a',i2,')')
 9660 format('(4(1x,f11.5),''| '',i6,'': '',a',i1,')')
 9700 format('Proc call at line ',i6) 
 9750 format(57x,'Average',6x,'Min. time',5x,'Max. time',4x,
     *         'Std. Dev.')
 9800 format('CALL ',a32,' at line ',i6,': ',4(f12.5,1x))
      end

      subroutine print_timers()
c--------------------------------------------------------------------------
c   Print the contents of all timers for debugging purposes.
c--------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'timerz.h'

      integer i, me, ierr
      
      call mpi_comm_rank(mpi_comm_world, me, ierr)
      print *,'Timer info on process ',me

      do i = 1, max_timers
         if (tdesc(i)(1:1) .ne. ' ') then
            print *,tdesc(i),' ',timers(i)
         endif
      enddo
      return
      end
     
      subroutine print_timer(str, key)
c------------------------------------------------------------------------
c   Prints an individual timer for debugging purposes.
c------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'timerz.h' 
      integer key, me, ierr
      character*(3) str

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      print *,str(1:3), ' Timer ',tdesc(key)(1:20),' on ',me,
     *                  ': ',timers(key)
      return
      end

      double precision function calc_sdev(sumsq, sumdata, contrib)
      double precision sumsq, sumdata, contrib

      double precision temp, dnom

      if (dabs(contrib-1.d0) .le. 0.d-6) then
         calc_sdev = 0.d0
         return
      endif

      temp = sumsq * contrib - sumdata * sumdata
      if (temp .lt. 0.d0) then
         calc_sdev = 0.d0
         return
      endif

      dnom = contrib * (contrib - 1.d0)
      calc_sdev = dsqrt(temp / dnom ) 
      return 
      end
