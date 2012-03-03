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
      subroutine pardo_loadb_update_batch(iop, optable, index_table, 
     *                segment_table, nsegment_table,
     *                next_batch, last_batch) 
c---------------------------------------------------------------------------
c   Gets the next batch of the current pardo instruction from the
c   pardo master of the current rank.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'trace.h'
      include 'where_table.h'

      integer iop, noptable, next_batch, last_batch
      integer optable(loptable_entry, 1)
      integer index_table(lindex_table_entry, *)
      integer nsegment_table
      integer segment_table(lsegment_table_entry, *)

      integer indx(mx_array_index)
      integer ind(mx_array_index)
      integer bseg(mx_array_index)
      integer eseg(mx_array_index)
      integer seg(mx_array_index)
      integer i, j, max_batch, nind
      integer chunk_size
      integer last_chunk
      integer iwhere_save, nwhere_save 
      integer cluster_size, get_pardo_cluster_size
      integer cluster_start, cluster_end
      integer pardo_master, get_pardo_master

c---------------------------------------------------------------------------
c   This pardo loop may not have been entered on the current processor
c   at this time.  If that is the case, we must execute the code that sets 
c   the maximum number of batches for the loop.
c---------------------------------------------------------------------------

c      print *,'Task ',me,' PARDO_LOADB_UPDATE_BATCH: iop ',iop
c      print *,'   pardo_info: ',(optable(c_pardo_info+i-1,iop),i=1,6)
c      call c_flush_stdout()
      if (optable(c_pardo_max_batch,iop) .eq. 0) then
         do i = 1, mx_array_index
            indx(i) = optable(c_ind1+i-1,iop)
         enddo
      endif

c----------------------------------------------------------------------------
c   Get the workload distribution for the current pardo cluster.
c----------------------------------------------------------------------------

      call set_pardo_workload(optable(c_pardo_max_batch,iop), 
     *                        cluster_start, cluster_end)
      cluster_size = get_pardo_cluster_size()

      if (optable(c_pardo_next_batch_start,iop) .eq. 0) 
     *    optable(c_pardo_next_batch_start,iop) = cluster_start

      if (optable(c_pardo_next_batch_start,iop) .eq. -1) then
      
c--------------------------------------------------------------------------
c   If c_pardo_signal = cluster_size, then each processor has previously
c   received a "-1" signal from this pardo loop.  This means we are starting
c   a fresh loop execution, and the pardo counters must be reset.
c--------------------------------------------------------------------------

         cluster_size = get_pardo_cluster_size()
         if (optable(c_pardo_signal,iop) .eq. cluster_size) then
            optable(c_pardo_batch,iop)  = 0
            optable(c_pardo_next_batch_start,iop) = cluster_start
            optable(c_pardo_batch_end,iop)       = 0  
c            print *,'LOOP RESTART CODE: ',
c     *         (optable(c_pardo_info+i-1,iop),i=1,6)
         else
            
c--------------------------------------------------------------------------
c   Bump the signal counter, return a -1.
c--------------------------------------------------------------------------

            optable(c_pardo_signal,iop)=optable(c_pardo_signal,iop)+1
            next_batch = -1
            last_batch = 0
c            print *,'ALT RETURN: iop ',iop,' pardo_info ',
c     *       (optable(c_pardo_info+i-1,iop),i=1,6)
            return
         endif
      endif

c--------------------------------------------------------------------------
c   A new batch is needed.  Update the counters.
c--------------------------------------------------------------------------

         optable(c_pardo_batch,iop) = 
     *                   optable(c_pardo_next_batch_start,iop)
         last_chunk = optable(c_pardo_chunk_size,iop)

         pardo_master = get_pardo_master() 
         optable(c_pardo_chunk_size,iop) =
     *          chunk_size(optable(c_pardo_next_batch_start,iop),
     *                     cluster_end, cluster_size,
     *                     last_chunk)
         if (my_company_rank .eq. pardo_master) then
            optable(c_pardo_chunk_size,iop) = 
     *        max(optable(c_pardo_chunk_size,iop)/4, 1)
         endif

         optable(c_pardo_batch_end,iop) = optable(c_pardo_batch,iop)+
     *                   optable(c_pardo_chunk_size,iop) - 1
         if (optable(c_pardo_batch_end,iop) .gt. cluster_end)
     *      optable(c_pardo_batch_end,iop) = cluster_end
         optable(c_pardo_next_batch_start,iop) = 
     *     optable(c_pardo_next_batch_start,iop) +
     *     optable(c_pardo_chunk_size,iop)

c--------------------------------------------------------------------------
c   Check to determine if the loop maximum has been exceeded.  If so, a -1
c   batch number is returned.
c--------------------------------------------------------------------------

      if (optable(c_pardo_batch,iop) .gt. cluster_end) then
        optable(c_pardo_batch,iop) = -1
        optable(c_pardo_signal,iop) = 1
        optable(c_pardo_next_batch_start,iop) = -1
      endif

c---------------------------------------------------------------------------
c   Set the return value.
c---------------------------------------------------------------------------

      next_batch = optable(c_pardo_batch,iop)    
      last_batch = optable(c_pardo_batch_end, iop)
c      print *,'Task ',me,' PARDO_LOADB_UPDATE_BATCH: return ',
c     *   'next_batch, last_batch ',next_batch,last_batch
c      print *,'   pardo_info at end of routine ',
c     *   (optable(c_pardo_info+i-1,iop),i=1,6)
c      call c_flush_stdout()
      return
      end


      subroutine pardo_loadb_get_my_batch(iop, next_batch, last_batch)
c--------------------------------------------------------------------------
c   Determines the address of the optable, uses this to call 
c   pardo_loadb_update_batch.
c--------------------------------------------------------------------------
      implicit none
      include 'sip_tables.h'
      include 'parallel_info.h'
      include 'interpreter.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer iop, next_batch, last_batch
      integer optable(1)
      integer index_table(1)
      integer segment_table(1)
#ifdef ALTIX
      pointer (ptr_optable, optable)
      pointer (ptr_index_table, index_table)
      pointer (ptr_segment_table, segment_table)
#else
      common optable
#endif

      integer next_batch_save, last_batch_save
      integer inext, ilast

      integer*8 iiopt, get_index_from_base
      integer*8 iindext
      integer*8 isegt
      integer i, ipardoinfo

#ifdef ALTIX
      ptr_optable = ishptr
      ptr_index_table = ishptr
      ptr_segment_table = ishptr
#endif

      iiopt = get_index_from_base(optable_base, optable, 1)
      iindext = get_index_from_base(index_table_base, index_table, 1)
      isegt = get_index_from_base(segment_table_base, segment_table, 1)

c------------------------------------------------------------------------
c   Save table values needed by the master (rank 0) process.
c------------------------------------------------------------------------

      inext = (iop-1)*loptable_entry + c_pardo_batch
      ilast = (iop-1)*loptable_entry + c_pardo_batch_end
      next_batch_save = optable(iiopt+inext-1)
      last_batch_save = optable(iiopt+ilast-1)

c--------------------------------------------------------------------------
c   Update table entries.
c--------------------------------------------------------------------------

c      print *,'SERVER CALL to pardo_loadb_update_batch: saved ',
c     *   'next_batch, last_batch = ',next_batch_save,last_batch_save
c      ipardoinfo = (iop-1)*loptable_entry + c_pardo_info 
c      print *,'Pardo_info from calling routine ',
c     *   (optable(iiopt+ipardoinfo+i-2),i=1,6)
      call pardo_loadb_update_batch(iop, optable(iiopt), 
     *         index_table(iindext), segment_table(isegt),
     *         nsegment_table_sip,
     *         next_batch, last_batch)

c---------------------------------------------------------------------------
c   Restore the master's private values.
c---------------------------------------------------------------------------

      optable(iiopt+inext-1) = next_batch_save
      optable(iiopt+ilast-1) = last_batch_save
c      print *,'SERVER CALL pardo_info after call = ',
c     *     (optable(iiopt+ipardoinfo+i-2),i=1,6)
      return
      end
