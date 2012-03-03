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
      subroutine restart(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

c----------------------------------------------------------------------------
c   Fine-grain restart from checkpointed data.
c   This routine is executed before any instructions.
c   1. We check to determine if this is a restart job. 
c   2. If so, we read data from the master checkpoint data file.
c   3. A message is sent to each SIAL server process to indicate that it
c      should set up its data from its own checkpoint data files.
c   4. We establish the SIAL program context from this data.
c   5. Then we jump to the instruction indicated in the checkpoint file.
c   6. Normal program execution continues from that point.   
c
c   If this is not a restart job, the normal instruction sequence is  
c   performed.
c----------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'dbugcom.h'
      include 'checkpoint_data.h'
      include 'sip_tables.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer icom(1)
#ifdef ALTIX
      pointer (iptr, icom)
#else
      common icom
#endif

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer status(MPI_STATUS_SIZE)
      integer msg(5)
      integer i, j
      integer ierr
      integer company_comm, pst_get_company_comm
      integer array, array_type
      integer nscalar_table_file, nindex_table_file

      integer pst_get_company
      integer niocompany
      integer request(10000)
      integer statuses(MPI_STATUS_SIZE,10000)
      integer context(3)
      integer*8 addr, iopaddr, get_index_from_base
      integer instruction_timer, comm_timer, flopcount
      integer current_seg(nindex_table)
      logical file_exist

#ifdef ALTIX
      iptr = ishptr
#endif

      restart_status = .false.

      if (restart_job) then
         company_comm = pst_get_company_comm(me)
         restart_status = .true.

         if (.not. mpi_io_support) return   ! cannot use ckpt/restart.

      inquire (file='master_ckpt_file', exist=file_exist)
      if (.not. file_exist) then
         restart_job    = .false.
         restart_status = .false.
         return
      endif

      if (me .eq. 0) then

c--------------------------------------------------------------------------
c   Open the master checkpoint data file.
c--------------------------------------------------------------------------

         open(file='master_ckpt_file', unit=master_ckpt_unit,
     *     ACTION='READ',access='SEQUENTIAL', form='UNFORMATTED',
     *     STATUS='OLD',
     *     err = 1000, iostat = ierr)

c--------------------------------------------------------------------------
c   Read the data.
c--------------------------------------------------------------------------

         read (master_ckpt_unit) current_line, current_op
         read (master_ckpt_unit) nactive_allocate_table
         read (master_ckpt_unit) nactive_create_table
         read (master_ckpt_unit) nckpt_arrays
         read (master_ckpt_unit) nscalar_table_file,nindex_table_file
         read (master_ckpt_unit) (context(i),i=1,3)

         if (nscalar_table .ne. nscalar_table_file .or.
     *       nindex_table_file .ne. nindex_table) then
            print *,'Error: Restart failed sanity check.'
            print *,'nscalar_table ',nscalar_table 
            print *,'nscalar_table_file ',nscalar_table_file 
            print *,'nindex_table ',nindex_table 
            print *,'nindex_table_file ',nindex_table_file 
            call abort_job()
         endif

         do i = 1, nactive_allocate_table
            read (master_ckpt_unit) active_allocate_table(i),
     *                active_allocate_op(i)
         enddo

         do i = 1, nactive_create_table
            read (master_ckpt_unit) active_create_table(i),
     *               active_create_op(i)
         enddo

         do i = 1, nckpt_arrays
            read (master_ckpt_unit) ckpt_arrays(i), ckpt_diskaddr(i)
         enddo
   
         do i = 1, nscalar_table
            read (master_ckpt_unit) scalar_table(i)
         enddo

         do i = 1, nindex_table
            read (master_ckpt_unit) current_seg(i)
         enddo

         call restore_instruction_stack(master_ckpt_unit)

c--------------------------------------------------------------------------
c   Close the file.
c--------------------------------------------------------------------------

         close(master_ckpt_unit)

c--------------------------------------------------------------------------
c   Send restart messages to each server.
c--------------------------------------------------------------------------

         niocompany = 0
         msg(1) = sip_server_restart_signal
         do i = 1, nprocs
            if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Send the checkpoint message.
c---------------------------------------------------------------------------

               niocompany = niocompany + 1
               call mpi_isend(msg, 4,
     *               MPI_INTEGER, i-1,
     *               sip_server_message, mpi_comm_world,
     *               request(niocompany), ierr)
            endif
         enddo
      endif   ! me .eq. 0   

c---------------------------------------------------------------------------
c   Set up the new SIAL program context on each worker process.
c---------------------------------------------------------------------------

      call broadcast_instruction_stack(company_comm)
      call mpi_bcast(context, 3, MPI_INTEGER, 0, company_comm, ierr)
      call set_program_context(context(1), context(2), context(3))

c--------------------------------------------------------------------------
c   Restore the "current_seg" field for each index in the index_table.
c--------------------------------------------------------------------------

      call mpi_bcast(current_seg, nindex_table, MPI_INTEGER, 0,
     *               company_comm, ierr)

      do i = 1, nindex_table
         index_table(c_current_seg,i) = current_seg(i)
      enddo

c---------------------------------------------------------------------------
c   Initialize each doloop in the instruction stack.
c---------------------------------------------------------------------------

      iopaddr = get_index_from_base(optable_base, icom, 1)
      call restore_do_loops(icom(iopaddr), noptable_sip, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      current_seg, context(2), context(3))

c--------------------------------------------------------------------------
c   Execute all the necessary "create" instructions.
c--------------------------------------------------------------------------

      call mpi_bcast(nactive_create_table, 1, MPI_INTEGER, 0, 
     *               company_comm, ierr)
      call mpi_bcast(active_create_table, nactive_create_table,
     *               MPI_INTEGER, 0, company_comm, ierr)
      call mpi_bcast(active_create_op, nactive_create_table,
     *               MPI_INTEGER, 0, company_comm, ierr)

      do i = 1, nactive_create_table

c--------------------------------------------------------------------------
c   Form the address of the instruction, then compute its index relative
c   to blank common.
c--------------------------------------------------------------------------

         addr = optable_base + loptable_entry *
     *         (active_create_op(i)-1)*intsize
         iopaddr = get_index_from_base(addr, icom, 1)

         if (icom(iopaddr) .ne. create_op) then
            print *,'Error: Invalid create instruction'
            print 10101,(icom(iopaddr+j-1),j=1,loptable_entry)
10101 format(10(1x,i6))
            print *,'noptable_sip ',noptable_sip 
            call abort_job()
         endif

c--------------------------------------------------------------------------
c   Execute the "create" instruction.
c--------------------------------------------------------------------------

         instruction_timer = 0
         comm_timer        = 0
         flopcount         = 0
         call compute_block(icom(iopaddr), 
     *                      array_table, narray_table,
     *                      index_table, nindex_table,
     *                      block_map_table, nblock_map_table,
     *                      segment_table, nsegment_table,
     *                      scalar_table, nscalar_table,
     *                      address_table,
     *                      .false., .false.,
     *                      flopcount, company_comm, comm_timer,
     *                      instruction_timer)
      enddo

c--------------------------------------------------------------------------
c   Execute all the necessary "allocate" instructions.
c--------------------------------------------------------------------------

      call mpi_bcast(nactive_allocate_table, 1, MPI_INTEGER, 0,
     *               company_comm, ierr)
      call mpi_bcast(active_allocate_table, nactive_allocate_table,
     *               MPI_INTEGER, 0, company_comm, ierr)
      call mpi_bcast(active_allocate_op, nactive_allocate_table,
     *               MPI_INTEGER, 0, company_comm, ierr)

      do i = 1, nactive_allocate_table

c--------------------------------------------------------------------------
c   Form the address of the instruction, then compute its index relative 
c   to blank common.
c--------------------------------------------------------------------------

         addr = optable_base + loptable_entry *
     *         (active_allocate_op(i)-1)*intsize
         iopaddr = get_index_from_base(addr, icom, 1)

c---------------------------------------------------------------------------
c   Execute the instruction.
c---------------------------------------------------------------------------

         call allocate_instruction(array_table, narray_table,
     *                   index_table,
     *                   nindex_table, segment_table, nsegment_table,
     *                   block_map_table, nblock_map_table, 
     *                   icom(iopaddr))
      enddo

c---------------------------------------------------------------------------
c   Check to see if any distributed data was checkpointed.
c---------------------------------------------------------------------------

      call mpi_bcast(nckpt_arrays, 1, MPI_INTEGER, 0,
     *               company_comm, ierr)
      call mpi_bcast(ckpt_arrays, mx_ckpt_arrays, MPI_INTEGER, 0,
     *               company_comm, ierr)
      call mpi_bcast(ckpt_diskaddr, mx_ckpt_arrays, 
     *               MPI_INTEGER8, 0, company_comm, ierr)

      call restore_distributed_data(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op, company_comm)
 
c---------------------------------------------------------------------------
c   Wait for completion of server restart signals.
c---------------------------------------------------------------------------

      if (me .eq. 0) then
         if (niocompany .gt. 0) then
            call mpi_waitall(niocompany, request,
     *                    statuses, ierr)
            call prt_time('All server restart signals were received.')

            niocompany = 0
            do i = 1, nprocs
               if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Post a recv for the "restart complete" message.
c   Each server will send a 1-word indicator message using the
c   sip_server_checkpoint_signal as the message tag.
c---------------------------------------------------------------------------

                  niocompany = niocompany + 1
                  call mpi_irecv(msg, 1,
     *               MPI_INTEGER, i-1,
     *               sip_server_restart_signal, mpi_comm_world,
     *               request(niocompany), ierr)
               endif
            enddo

c-----------------------------------------------------------------------
c   Wait for all the acknowledgement messages.
c------------------------------------------------------------------------

            call mpi_waitall(niocompany, request, statuses,
     *                    ierr)
         endif   ! niocompany .gt. 0

         print *,'RESTARTING JOB AT LINE ',current_line
      endif   ! me .eq. 0

 1000    continue
         call mpi_barrier(company_comm, ierr)
      endif   ! restart_job

      restart_job = .false.
      return
      end

      subroutine restore_distributed_data(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op, company_comm)
c---------------------------------------------------------------------------
c   Restores each distributed array that was checkpointed data in a previously
c   run job.  The distributed arrays must have been created already.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'dbugcom.h'
      include 'checkpoint_data.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)
      integer company_comm

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      integer array, iarray, array_type
      integer mode, fh, ierr
      integer i, k, stack, blkndx, val1, val2, size
      integer nblk, nindex,iblock_map, iblock, nreq
      integer ind(mx_array_index), seg(mx_array_index)
      integer mx_mpi_req
      parameter (mx_mpi_req = 10000)
      integer request(mx_mpi_req)
      integer statuses(MPI_STATUS_SIZE, mx_mpi_req)
      integer*8 offset, indblk
      integer*8 get_block_index

      logical file_is_open

#ifdef ALTIX
      dptr = dshptr
#endif

      file_is_open = .false.
      do iarray = 1, nckpt_arrays
         array = ckpt_arrays(iarray)
         array_type = array_table(c_array_type, array)
         nindex     = array_table(c_nindex, array)
         if (array_type .eq. distributed_array) then

c----------------------------------------------------------------------------
c   Open the CKPT_DIST file.
c----------------------------------------------------------------------------

            if (.not. file_is_open) then
               file_is_open = .true.
               mode = MPI_MODE_RDONLY
               call mpi_file_open(company_comm, 'CKPT_DIST', mode,
     *                   MPI_INFO_NULL, fh, ierr)
               offset = 0
               call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION,
     *                       MPI_DOUBLE_PRECISION, 'native',
     *                       MPI_INFO_NULL, ierr)
            endif   ! .not. file_is_open   

            offset = ckpt_diskaddr(iarray)

            nblk = array_table(c_numblks, array)
            iblock_map = array_table(c_block_map, array)
            nreq = 0

            do iblock = 1, nblk

c--------------------------------------------------------------------------
c   Determine block size in words.
c--------------------------------------------------------------------------

               size = 1
               do k = 1, nindex
                  ind(k) = array_table(c_index_array1+k-1,array)
                  seg(k) = block_map_table(c_block_map_seg+k-1,
     *                       iblock_map+iblock-1)
                  call get_index_segment(ind(k), seg(k), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
                  size = size * (val2-val1+1)
               enddo

               blkndx = 
     *           block_map_table(c_bmap_blkndx,iblock_map+iblock-1)
               if (block_map_table(c_processor,iblock_map+iblock-1) 
     *             .eq. my_company_rank) then

c----------------------------------------------------------------------------
c   The block resides on this processor.   Read it in from disk.
c----------------------------------------------------------------------------

                  stack  = array_table(c_array_stack,array)
                  indblk = get_block_index(array, iblock, stack, blkndx,
     *                               x, .true.)

                  if (nreq .eq. mx_mpi_req) then
                     call mpi_waitall(nreq, request, statuses, ierr)
                     nreq = 0
                  endif

                  nreq = nreq + 1
                  call mpi_file_iread_at(fh, offset, x(indblk), size,
     *                              MPI_DOUBLE_PRECISION, request(nreq),
     *                              ierr)
               endif  ! block_map_table(c_processor...

               offset = offset + size
            enddo   ! iblock
         endif   !  array_type .eq. distributed_array
      enddo   ! iarray

c----------------------------------------------------------------------------
c   Wait on MPI requests to complete.
c---------------------------------------------------------------------------

      if (nreq .gt. 0) call mpi_waitall(nreq, request, statuses, ierr)

c----------------------------------------------------------------------------
c   Close file.
c----------------------------------------------------------------------------

      if (file_is_open) call mpi_file_close(fh, ierr)

      return
      end
