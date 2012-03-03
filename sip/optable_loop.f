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
      subroutine optable_loop(optable, noptable, array_table, 
     *                   narray_table,
     *                   array_labels, index_table, nindex_table, 
     *                   segment_table, nsegment_table, block_map_table,
     *                   nblock_map_table, 
     *                   scalar_table, nscalar_table, proctab,
     *                   address_table, 
     *                   debug, validate, comm, comm_timer)
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_barrier_data.h'
      include 'scratchpad.h'
      include 'dbugcom.h'
      include 'timerz.h'
      include 'where_table.h'
      include 'context.h'
      include 'checkpoint_data.h'

      common /load_balance/load_balance
      logical load_balance

      include 'int_gen_parms.h'
       
      integer noptable, narray_table, nindex_table, nsegment_table
      integer nblock_map_table, nscalar_table
      integer comm
      integer optable(loptable_entry,noptable)
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer block_map_table(lblock_map_entry,nblock_map_table)
      integer proctab(2,*)
      character*10 array_labels(narray_table)
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)
      
      integer ierr
      integer iopsave, iblk
      integer index, result_array, result_type
      integer block_map_entry(lblock_map_entry)
      integer op1_block_map_entry(lblock_map_entry)
      integer op2_block_map_entry(lblock_map_entry)
      integer opcode
      integer i, j, k
      integer instruction_timer
      integer icount, timer_count 

      integer iseg, jseg

      double precision blk_send_time
      double precision timeslice, t1, t2
      double precision t1_heartbeat, heartbeat_time

      integer flopcount, blk_send_count
      integer comm_timer

      logical debug, validate

      character*40 hbmsg

      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      call mpi_comm_rank(mpi_comm_world, me, ierr)
      call mpi_comm_rank(comm, my_company_rank, ierr)
      call mpi_comm_size(comm, my_company_size, ierr)

      trace = .false.
      load_balance = .true.

      call determine_timeslice(timeslice, timer_count)

      iwhere = 0
      nwhere = 0

c--------------------------------------------------------------------------
c   Initialize the server barrier data requests.
c--------------------------------------------------------------------------

      do i = 1, mx_msg
         server_requests(i) = MPI_REQUEST_NULL
      enddo

c--------------------------------------------------------------------------
c   Initialize the array_table status values.
c--------------------------------------------------------------------------

      do i = 1, narray_table
         array_table(c_array_status,i) = 0
      enddo
c--------------------------------------------------------------------------
c   Set the current segments of each index in the index_table to an 
c   undefined value.  If an array is later referenced without having all
c   its indices defined, an error will be triggered.
c--------------------------------------------------------------------------

      do i = 1, nindex_table
         index_table(c_current_seg,i) = undefined_segment
      enddo
      
c--------------------------------------------------------------------------
c   Clear the scratchpad areas.
c--------------------------------------------------------------------------

      do i = 1, mx_scratchpad
         scratchpad(i) = 0
      enddo

      do i = 1, mx_fscratchpad
         fscratchpad(i) = 0.
      enddo

c--------------------------------------------------------------------------
c   Main processing loop.
c--------------------------------------------------------------------------

         call set_program_context(1, 0, noptable)

         t1 = mpi_wtime()
         t1_heartbeat = t1
         heartbeat_time = 360.0  ! print every 6 minutes
         icount = 0

c---------------------------------------------------------------------------
c   Restart from a previous job if necessary.
c---------------------------------------------------------------------------

         call restart(array_table, narray_table,
     *                index_table,
     *                nindex_table, segment_table, nsegment_table,
     *                block_map_table, nblock_map_table,
     *                scalar_table, nscalar_table,
     *                address_table, optable(1,iop))

 1000    continue
            opcode = optable(c_opcode, iop)
            iopsave = iop
            current_op = iop  ! save operation pointer for debugging  purposes.
            current_line = optable(c_lineno,iop)
c            if (trace .and. 
c     *          and(tracelevel, instruction_trace) .ne. 0) then
c               print *,'Task ',me,' Perform op = ',(optable(k,iop),
c     *           k=1,loptable_entry),' segments: ',
c     *            (index_table(c_current_seg,k),k=1,nindex_table),
c     *            ' start_op, end_op = ',start_op,end_op,' iop = ',
c     *            iop,' line ',optable(c_lineno,iop)
c            endif

            if (opcode .eq. call_op) then
               call handle_call(optable, noptable, proctab, debug,
     *                start_op, end_op, iop)
            else if (opcode .eq. return_op) then
               call handle_return(optable, noptable, debug,
     *                start_op, end_op, iop)
            else if (opcode .eq. go_to_op .or.
     *               opcode .eq. jz_op) then
               call handle_goto(optable, noptable, debug,
     *                start_op, end_op, iop)
            else if (opcode .eq. do_op .or.
     *               opcode .eq. enddo_op) then
               call doloop(optable, noptable, iop, index_table, 
     *                  nindex_table, array_table, narray_table,
     *                  block_map_table, segment_table,
     *                  nsegment_table, 
     *                  debug, .false.,
     *                  start_op, end_op)
            else if (opcode .eq. pardo_op .or.
     *               opcode .eq. endpardo_op) then
               if (load_balance) then
                  call pardo_loadb(optable, noptable, iop, index_table,
     *                  nindex_table, array_table, narray_table,
     *                  block_map_table, segment_table,
     *                  nsegment_table,
     *                  comm, debug, .false.,
     *                  start_op, end_op)
               else
                  call pardo_loop(optable, noptable, iop, index_table,
     *                  nindex_table, array_table, narray_table,
     *                  block_map_table, 
     *                  segment_table, nsegment_table,
     *                  comm, debug, .false.,
     *                  start_op, end_op)
               endif
            else if (opcode .eq. exit_op) then
               call handle_exit(optable, noptable, debug,
     *                array_table, narray_table,
     *                index_table, nindex_table, 
     *                block_map_table,    
     *                start_op, end_op, iop)
            else if (opcode .eq. cycle_op) then
               call handle_cycle(optable, noptable, debug,
     *                start_op, end_op, iop)
            endif
c            print *,'AFTER DO_LOOP: start_op, end_op = ',
c     *             start_op, end_op

            if (iop .gt. noptable) go to 2000
            if (iopsave .ne. iop) then
               go to 1000
            endif

            if (iop .lt. start_op .or.
     *          iop .gt. end_op) go to 900 

            current_op = iop  ! save operation pointer for debugging  purposes.
            current_line = optable(c_lineno,iop)

c---------------------------------------------------------------------------
c   Perform the operation.
c---------------------------------------------------------------------------

            if (icount .ge. timer_count) then
               icount = 0
               t2 = mpi_wtime()

               if (t2-t1_heartbeat .gt. heartbeat_time) then
                  if (me .eq. 0) then
                     write (hbmsg, 2100) current_line
                     dbg = .true.
                     call prt_time(hbmsg)
                     dbg = .false.
                  endif
                  t1_heartbeat = t2 
               endif

               if (t2-t1 .gt. timeslice) then
c                  print *,'Task ',me,' Timeslice interrupt'
c                  call c_flush_stdout()
                  call exec_thread_server(0)
c                  print *,'Task ',me,' Back from Timeslice interrupt'
c                  call c_flush_stdout()

                  t1 = mpi_wtime()   ! reset timer
               endif
            endif 

            icount = icount + 1

            if (trace) then
               print *,'Task ',me,' line ',current_line,' iop ',
     *                 iop,' opcode = ',optable(c_opcode, iop)
               call c_flush_stdout()
            endif

            if (do_timer) then
               if (optable(c_opcode,iop) .ne. pardo_op) then
                  instruction_timer = optable(c_instr_timer,iop)
                  call timer_start(instruction_timer)
               endif
            endif

	    call compute_block(optable(1,iop), array_table, 
     *                narray_table, index_table, nindex_table, 
     *                block_map_table, nblock_map_table,
     *                segment_table, nsegment_table,
     *                scalar_table, nscalar_table, address_table,
     *                debug, validate, 
     *                flopcount, comm, comm_timer, 
     *                instruction_timer)
             if (do_timer) then
                if (optable(c_opcode,iop) .ne. pardo_op)
     *             call update_timer(instruction_timer)
             endif 
  900    continue
         iop = iop + 1 
c            print *,'AT 900: iop, start_op, end_op = ',
c     *         iop, start_op, end_op
         if (start_op .ne. 0 .or. iop .le. noptable) then
            if (iop .gt. end_op)  iop = start_op
c            print *,'   BRANCH TO iop = ',iop
            if (iop .le. 0 .or. iop .gt. noptable) then
               print *,'ERROR: iop out of range: iop = ',iop
               call abort_job()
            endif
            go to 1000
         endif
         
 2000    continue
      if (dbg) print *,'Processing of table is complete...'
      call server_takedown(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, optable )
      return
 2100 format('Heartbeat: line ',i6)
      end
