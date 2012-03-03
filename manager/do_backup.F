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
      subroutine do_backup(iblock, server_table, nserver_table)
c------------------------------------------------------------------------------
c   Finds the oldest dirty block that is not busy and backs it up to disk.
c-----------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'server_stat.h'
      include 'mpif.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer iblock
      integer i
      integer ptr, diskloc, ifile, size
      integer*8 indblk, get_index_from_base
      integer f_backupram
      integer find_free_diskloc
      integer istat
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif
      double precision tbackup

#ifdef ALTIX
      dptr = dshptr   ! set x pointer to shared mem.
#endif

      call find_oldest_dirty_block(iblock, server_table,
     *            nserver_table)
      if (iblock .gt. 0) then
         if (do_stats .and. stat_key .gt. 0) tbackup = mpi_wtime() 
         ptr = server_table_ptr(iblock)
         if (ptr .le. 0) then
            print *,'Server ',me,' Error: Dirty block ',iblock,
     *          ' found with ptr =',ptr
            call server_abort_job(server_table, nserver_table)
         endif

         diskloc = server_table(c_server_diskloc,ptr)
         ifile   = server_table(c_server_file,ptr)
 
         if (diskloc .eq. 0) then
            diskloc = find_free_diskloc(ifile, server_table, 
     *                  nserver_table)  ! check for a previously freed block
            if (diskloc .lt. 0) then
               next_server_diskloc(ifile) = 
     *                       next_server_diskloc(ifile) + 1
               server_table(c_server_diskloc,ptr) = 
     *                             next_server_diskloc(ifile)
               diskloc = next_server_diskloc(ifile)
            else
               server_table(c_server_diskloc,ptr) = diskloc
            endif
         endif
                 
         indblk = get_index_from_base(base_mem_addr, x, 2)
         indblk = (iblock-1)*server_mem_blocksize + indblk
         size = server_table(c_server_size,ptr)

         istat =  f_backupram(server_unit(ifile), diskloc, 
     *                 server_blocksizes(ifile), 
     *                 x(indblk), size)
         if (istat .ne. 0) 
     *     call server_abort_job(server_table, nserver_table)
         call mark_block_clean(iblock, server_table,
     *                      nserver_table)
         if (do_stats .and. stat_key .gt. 0) then
            tbackup = mpi_wtime() - tbackup 
            sstat_tbackup(stat_key) = sstat_tbackup(stat_key) + 
     *                         tbackup   
            sstat_tbackup2(stat_key) = sstat_tbackup2(stat_key) +
     *                         tbackup * tbackup
            sstat_nbackup(stat_key) = sstat_nbackup(stat_key) + 1
         endif
      endif

      return
      end
