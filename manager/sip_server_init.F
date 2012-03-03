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
      subroutine sip_server_init(local_scratch, x, nx, msg_buffers,
     *                           blocksize, niocompany, do_timer)
c---------------------------------------------------------------------------
c   Initialization routine for the SIP server process.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'server_ckpt_data.h'
      include 'server_stat.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'machine_types.h'
      include 'dbugcom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer*8 c_loc64, ixx

      character*(*) local_scratch
      integer nx, msg_buffers, blocksize
      integer niocompany
      logical do_timer

      integer sip_server_datatag
      parameter (sip_server_datatag = 3335)
      integer status(mpi_status_size)

      integer i, nleft
      double precision x(nx)
      character*256 fn
      integer n
      integer str_trimlen
      integer ios
      integer server_table(1)
      integer*8 indtable, get_index_from_base
      integer nwtable, ndouble, nbtable, ierr
#ifdef ALTIX
      pointer (iptr, server_table)
#endif

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)

#ifdef ALTIX
      iptr = ishptr
#endif

c----------------------------------------------------------------------------
c   Receive the server_table from the master.
c----------------------------------------------------------------------------

      call mpi_recv(x, nx, mpi_integer, MPI_ANY_SOURCE, 
     *       sip_server_datatag, mpi_comm_world, status, ierr)

      call mpi_get_count(status, mpi_integer, nwtable, ierr)

c---------------------------------------------------------------------------
c   The second message from the master gives us the number of stacks used
c   and their respective blocksizes.
c----------------------------------------------------------------------------

      call mpi_recv(server_blocksizes, mx_server_blocksizes,
     *       mpi_integer, MPI_ANY_SOURCE,
     *       sip_server_datatag, mpi_comm_world, status, ierr)

      call mpi_get_count(status, mpi_integer, nserver_blocksizes, ierr)

c----------------------------------------------------------------------------
c   The server table is at the beginning of our allocated memory.
c   The memory blocks go at the end of the server table.
c----------------------------------------------------------------------------

      ixx = 1
      server_table_base_addr = c_loc64(x, ixx, bytes_per_double)
      nserver_table_entries = nwtable/lserver_table_entry

c----------------------------------------------------------------------------
c   Increase nwtable to be a multiple of the size of a double word.
c----------------------------------------------------------------------------

      nbtable = (nwtable * intsize + bytes_per_double - 1)/
     *           bytes_per_double * bytes_per_double
      ndouble = nbtable / bytes_per_double
      
c---------------------------------------------------------------------------
c   Calculate the base address of the memory blocks.
c---------------------------------------------------------------------------

      indtable = get_index_from_base(server_table_base_addr, 
     *                               server_table, 1)
      base_mem_addr = c_loc64(server_table, indtable+nbtable/intsize, 
     *                        intsize)
      nleft = nx - ndouble   ! num. of remaining doubles.

      nmessage_buffers      = msg_buffers
      server_mem_blocksize  = blocksize
      server_seqno          = 0
      server_work_list_head = 0
      server_work_list_tail = 0
      server_node_ptr       = 0
      dirty_list_head       = 0
      dirty_list_tail       = 0
      nserver_memblocks     = nleft/blocksize
      nserver_memblocks     = min(nserver_memblocks,mx_server_memblocks)
      dirty_threshold       = (nserver_memblocks-msg_buffers)/10
      max_backup            = min(msg_buffers/2,dirty_threshold)
      barrier_in_progress   = .false.
      barrier_msg_count     = 0
      nbarrier_msgs         = nprocs - niocompany

      do_stats = do_timer
      if (do_stats) call clear_server_stats()   ! initialize server statistics

      if (nmessage_buffers .ge. nserver_memblocks) then
         print *,'Error: Server requires at least ',nmessage_buffers,
     *     ' message buffers, but only ',nserver_memblocks,
     *     ' are available.'
         call abort_job()
      endif

      if (dbg) then
         print *,'Server ',me,' was initialized with ',
     *     nserver_memblocks,
     *    ' blocks of ',blocksize,' words each.'
         print *,'Backup threshold ',dirty_threshold,
     *         ' max_backup ',max_backup
      endif

c----------------------------------------------------------------------------
c   Open the server's scratch files.
c----------------------------------------------------------------------------

      do i = 1, nserver_blocksizes
         fn = ' '
         call generate_scratch_filename(fn, 'server'//char(0))
         n = str_trimlen(local_scratch)
         if (local_scratch(n:n) .eq. '/') then
            server_filename(i) = local_scratch(1:n)//fn
         else
            server_filename(i) = local_scratch(1:n) // '/' // fn
         endif

         n = str_trimlen(server_filename(i))
         call f_creatfile(server_filename(i)(1:n) // char(0), 
     *                       server_unit(i))
         if (dbg) print *,'Server ',me,' Opened file ',
     *              server_filename(i)(1:n),' file handle ',
     *              server_unit(i) 
         next_server_diskloc(i) = 0
      enddo

c----------------------------------------------------------------------------
c   Set up checkpoint files for fine-grained restart.  These files are 
c   only opened if a checkpoint message is processed by the server.
c----------------------------------------------------------------------------

      write (ckpt_ndx_filename, 100) me,'.ndx'
  100 format('./CKPT',i5,A4)
      write (ckpt_dat_filename, 100) me, '.dat'

c-----------------------------------------------------------------------------
c   Change blanks in filenames to 0's.
c-----------------------------------------------------------------------------

      call convert_chars(ckpt_ndx_filename,' ', '0')
      call convert_chars(ckpt_dat_filename,' ', '0')

c----------------------------------------------------------------------------
c   Initialize the message buffers.
c----------------------------------------------------------------------------

      do i = 1, nmessage_buffers
         server_msg(c_msg_msgbuffer,i) = i   ! fill in msg buffer field
         call push_work_node(i)
         server_table_ptr(i) = -1            ! mark block as msg_buffer
      enddo

      clean_block_ptr = 0  
      nclean_blocks = 0 

c---------------------------------------------------------------------------
c   Remaining blocks are all set up as "clean", i. e. ready for use.
c---------------------------------------------------------------------------

      do i = nmessage_buffers+1, nserver_memblocks
         server_table_ptr(i) = 0
         call push_clean_block(i, server_table, nserver_table_entries) 
      enddo

      return
      end

     
