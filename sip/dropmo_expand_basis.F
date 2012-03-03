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
      subroutine dropmo_expand_basis(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Dumps all blocks of an array to a list.
c   Uses new format (i. e. written with MPI_IO) BLOCK_INDEX/BLOCKDATA files.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blockdata_list.h'
      include 'dbugcom.h'
      include 'blkmgr.h'
      include 'dropmo.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'sip_tables.h'
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

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
      pointer (iptr, esegtable)
      pointer (iptr2, eindex_table)
      pointer (dptr2, scr_dropmo)
#else
      common x, scr_dropmo, esegtable, eindex_table
#endif
    
      integer proc, company_comm
      integer pst_get_company_comm, pst_get_company_rank
       integer*8 get_index_from_base
      integer*8 nmem_bytes
      integer   nmem_words
      integer*8 iblkmgr, idropmo
      integer*8 isegbase_exp, iexp_segtable, iexp_index_table
      integer   ensegment_table, ierr, ndropmo_blksize
      integer esegtable(1), eindex_table(1)
      double precision scr_dropmo(1)

      integer nestacks
      parameter (nestacks = max_stacks)
      integer estack_blocksizes(nestacks)

#ifdef ALTIX
      dptr = dshptr
      dptr2 = dshptr
      iptr = ishptr
      iptr2 = ishptr
#endif

      if (ndropmo_a .eq. 0 .and. ndropmo_b .eq. 0) return

      company_comm    = pst_get_company_comm(me)
      proc  = pst_get_company_rank(me)
      if (proc .eq. 0) then
         call mem_alloc_free(first_blk_addr, ierr)  ! release blkmgr memory.
         if (ierr .ne. 0) then
            print *,'Task ',me,' Error in mem_alloc_free'
            call abort_job()
         endif
 
         ndropmo_blksize = stack_blocksize(nblkmgr_stacks) ! save largest blksize
         call mem_alloc(scr_dropmo, ndropmo_blksize, bytes_per_double,
     *                  idropmo, .true., ierr)
         if (ierr .ne. 0) then
            print *,'Task ',me,' mem_alloc error in dropmo_expand_basis'
            print *,'Attempted to malloc ',ndropmo_blksize,' doubles.'
            call abort_job()
         endif  

         iexp_index_table = get_index_from_base(eindex_table_base, 
     *                       eindex_table, 1)

c---------------------------------------------------------------------------
c   Determine the proper segmentation that will be set up for the data
c   in the expanded basis set.
c---------------------------------------------------------------------------

         call get_expanded_basis_segsizes(enalpha_occupied,
     *           enalpha_virtual, enbeta_occupied, enbeta_virtual)

c---------------------------------------------------------------------------
c   Build a new segment table based on the expanded basis.
c---------------------------------------------------------------------------

         call build_segment_table(eindex_table(iexp_index_table), 
     *             nindex_table,
     *             end_nfps, nshells, sip_mx_ao_segsize,
     *             sip_mx_occ_segsize, sip_mx_virt_segsize,
     *             sip_sub_segsize, sip_sub_occ_segsize,
     *             sip_sub_virt_segsize, sip_sub_ao_segsize,
     *             enbeta_occupied, enalpha_occupied, enbeta_occupied,
     *             enalpha_virtual, enalpha_virtual, enbeta_virtual,
     *             isegbase_exp,ensegment_table, estack_blocksizes,
     *             nestacks, array_table, narray_table, .false.,
     *             atom)

          iexp_segtable = get_index_from_base(isegbase_exp, 
     *                                esegtable, 1)

          if (mpi_io_support) then
             call dropmo_expand_process_data_mpi_io(array_table, 
     *            narray_table,
     *            index_table, nindex_table,segment_table,
     *            nsegment_table, block_map_table, nblock_map_table,
     *            eindex_table(iexp_index_table), 
     *            esegtable(iexp_segtable), ensegment_table,
     *            scr_dropmo(idropmo), ndropmo_blksize)
          else
             call dropmo_expand_process_data_no_mpi_io(array_table,
     *            narray_table,
     *            index_table, nindex_table,segment_table,
     *            nsegment_table, block_map_table, nblock_map_table,
     *            eindex_table(iexp_index_table),
     *            esegtable(iexp_segtable), ensegment_table,
     *            scr_dropmo(idropmo), ndropmo_blksize)
          endif
      endif

      return
      end

      subroutine dropmo_expand_process_data_mpi_io(array_table, 
     *            narray_table,
     *            index_table, nindex_table,segment_table,
     *            nsegment_table, block_map_table, nblock_map_table,
     *            eindex_table, esegtable, ensegment_table,
     *            scr_dropmo, ndropmo_blksize)
c--------------------------------------------------------------------------
c   Performs the actual processing of the dropmo expansion for each 
c   array in the blockdata_list array.
c   Uses new format (i. e. MPI_IO-written) BLOCK_INDEX/BLOCKDATA files.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blockdata_list.h'
      include 'block_index.h'
      include 'dbugcom.h'
      include 'blkmgr.h'
      include 'dropmo.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table, ensegment_table
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer eindex_table(lindex_table_entry, nindex_table)
      integer esegtable(lsegment_table_entry, ensegment_table)

      integer ndropmo_blksize
      double precision scr_dropmo(ndropmo_blksize)

      integer list_unit, list_unit2, data_unit, data_unit2
      integer*8 data_ptr, data_ptr2, next_diskaddr
      integer ios
      logical is_file_open
 
      integer i, j, n, n_eblocks, max_eblocks, max_eblksize
      integer nindex, array, darray, iarray, ierr, map
      integer memloc
      integer dnindex, dblock
      integer index(mx_array_index), seg(mx_array_index)
      integer index_type(mx_array_index)
      integer dval1(mx_array_index), dval2(mx_array_index)  
      integer bseg(mx_array_index), eseg(mx_array_index), 
     *        nseg(mx_array_index)
      integer e_bseg(mx_array_index), e_eseg(mx_array_index), 
     *        e_nseg(mx_array_index)
      integer val1(mx_array_index), val2(mx_array_index), iblk,
     *        ieblk, nxeblocks
      integer*8 ixeblocks, i8size, malloced_previously
      integer*8 ieb, idiskloc, ixx
      integer*8 index_file_loc, index_file_loc2
      integer block_index_entry(lblk_index_entry)

      integer eblock_table(1)
      integer*8 diskloc_table(1)
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
      pointer (iptr, eblock_table)
      pointer (i8ptr, diskloc_table)
#else
      common x, diskloc_table, eblock_table
#endif
      logical restart

#ifdef ALTIX
      dptr = dshptr
      iptr = ishptr
      i8ptr = i8shptr
#endif

c---------------------------------------------------------------------------
c   Open files.
c---------------------------------------------------------------------------

      call f_openfile('BLOCK_INDEX'//char(0), list_unit)
      call f_openfile('BLOCKDATA' // char(0), data_unit)
      data_ptr = 0
      call f_creatfile('BLOCK_INDEX2'//char(0), list_unit2)
      call f_creatfile('BLOCKDATA2' // char(0), data_unit2)
      data_ptr2 = 0

c-------------------------------------------------------------------------
c   Determine max. number of expanded blocks, max. expanded blocksize.
c-------------------------------------------------------------------------

      max_eblocks = 0
      max_eblksize = 0
      do iarray = 1, narray_list
         array = array_list(iarray)
         nb = array_table(c_numblks, array)

         nindex = array_table(c_nindex, array)
         
         n_eblocks = 1
         do i = 1, nindex
            index(i) = array_table(c_index_array1+i-1,array)
            bseg(i)  = index_table(c_bseg, index(i))
            eseg(i)  = index_table(c_eseg, index(i)) 
            nseg(i)  = index_table(c_nsegments, index(i))

            e_bseg(i) = eindex_table(c_bseg, index(i)) 
            e_eseg(i) = eindex_table(c_eseg, index(i)) 
            e_nseg(i) = eindex_table(c_nsegments, index(i)) 

            n_eblocks = n_eblocks * e_nseg(i)
         enddo

c         print *,'DROPMO_EXPAND_BASIS: array ',array
c         print *,'number of blocks before expansion ',nb
c         print *,'number of blocks after expansion ',n_eblocks
         max_eblocks = max(max_eblocks, n_eblocks) 

c-------------------------------------------------------------------------
c   Loop over all blocks to determine max blocksize.
c-------------------------------------------------------------------------

         do i = 1, nindex
            seg(i) = e_bseg(i)
         enddo
         
         do iblk = 1, n_eblocks
            n = 1
            do i = 1, nindex
               call get_index_segment(index(i), seg(i), esegtable,
     *                             ensegment_table, eindex_table,
     *                             nindex_table, val1, val2)

               n = n * (val2(1) - val1(1) + 1) 
            enddo
      
            max_eblksize = max(max_eblksize, n)

c--------------------------------------------------------------------------
c   Increment segment indices for next block.
c--------------------------------------------------------------------------

            if (i .lt. n_eblocks) then 
               i = 1
  100          continue
               seg(i) = seg(i) + 1 
               if (seg(i) .gt. e_nseg(i)) then
                  seg(i)= e_bseg(i)
                  i = i + 1
                  go to 100
               endif
            endif
         enddo
 
      enddo    ! iarray

c      print *,'max blocks for all arrays ',max_eblocks
c      print *,'max blksize for all eblocks ',max_eblksize

c-------------------------------------------------------------------------
c   Allocate memory.
c-------------------------------------------------------------------------

      call mem_alloc(eblock_table, max_eblocks * leblock_entry, 
     *               intsize, ieb, .true., ierr)

      call mem_alloc(diskloc_table, max_eblocks, 8, idiskloc, 
     *               .true., ierr)
      
      call mem_alloc_query(malloced_previously)
      i8size = maxmem * 1024
      i8size = i8size * 1024
      n = (i8size - malloced_previously) / bytes_per_double
      call mem_alloc(x, n, bytes_per_double, ixeblocks, .true., ierr)
      nxeblocks = n / max_eblksize
 
c-------------------------------------------------------------------------
c   Main processing loop.
c-------------------------------------------------------------------------

      next_diskaddr = 0
      index_file_loc = 0
      index_file_loc2 = 0
      do iarray = 1, narray_list
         array  = array_list(iarray)
         nb     = array_table(c_numblks, array)
         nindex = array_table(c_nindex, array)
         map    = array_table(c_block_map, array)
c         print *,'Processing array ',iarray,' nb ',nb,' nindex ',
c     *         nindex,' indices = ',
c     *   (array_table(c_index_array1+i-1,array),i=1,nindex)

         n_eblocks = 1
         do i = 1, nindex
            index(i) = array_table(c_index_array1+i-1,array)
            index_type(i) = index_table(c_index_type, index(i))
            bseg(i)  = index_table(c_bseg, index(i))
            eseg(i)  = index_table(c_eseg, index(i))
            nseg(i)  = index_table(c_nsegments, index(i))

            e_bseg(i) = eindex_table(c_bseg, index(i))
            e_eseg(i) = eindex_table(c_eseg, index(i))
            e_nseg(i) = eindex_table(c_nsegments, index(i))

            n_eblocks = n_eblocks * e_nseg(i)
         enddo

c         print *,'e_bseg ',(e_bseg(i),i=1,nindex)
c         print *,'e_eseg ',(e_eseg(i),i=1,nindex)
c         print *,'e_nseg ',(e_nseg(i),i=1,nindex)

c-------------------------------------------------------------------------
c   Build the eblock tables.
c-------------------------------------------------------------------------

c         print *,'next_diskaddr for array ',array ,' is ',
c     *       next_diskaddr
         call build_eblock_table(eblock_table(ieb), max_eblocks,
     *         n_eblocks, e_bseg, e_eseg, e_nseg, nindex, 
     *         diskloc_table(idiskloc), index, esegtable, 
     *         ensegment_table, eindex_table, nindex_table,
     *         next_diskaddr)
       
c         print *,'Processing array ',array,' n_eblocks = ',n_eblocks

         restart = .true.
         do iblk = 1, nb

c-------------------------------------------------------------------------
c   Read a descriptor record from the BLOCK_INDEX file.
c-------------------------------------------------------------------------

            call f_read_disk_ints(list_unit, index_file_loc,
     *              block_index_entry, lblk_index_entry)
            index_file_loc = index_file_loc + lblk_index_entry

            darray =  block_index_entry(c_blk_index_array)
            dblock = block_index_entry(c_blk_index_blkno)
            dnindex = block_index_entry(c_blk_index_nind)
            n = block_index_entry(c_blk_index_size)
            do i = 1, nindex
               dval1(i) = block_index_entry(c_blk_index_bsegs+i-1)
               dval2(i) = block_index_entry(c_blk_index_esegs+i-1)
            enddo

c--------------------------------------------------------------------------
c   Read the data block from disk.
c--------------------------------------------------------------------------

            call f_read_disk(data_unit, data_ptr, scr_dropmo, n)
            data_ptr = data_ptr + n

c---------------------------------------------------------------------------
c   Form the contribution of this (dropmo) block into all eblocks to which
c   it contributes.
c---------------------------------------------------------------------------

            call contribute_eblocks(scr_dropmo, dval1, dval2, nindex,
     *                eblock_table(ieb), n_eblocks, 
     *                diskloc_table(idiskloc), data_unit2, 
     *                x(ixeblocks), nxeblocks, max_eblksize, 
     *                index_type, restart)
            restart = .false.
         enddo   ! iblk

c-------------------------------------------------------------------------
c   Write BLOCK_INDEX2 data for the array and flush out any unwritten
c   eblocks to BLOCKDATA2.
c-------------------------------------------------------------------------

         do ieblk = 1, n_eblocks

c--------------------------------------------------------------------------
c   Write out the BLOCK_INDEX2 data for the eblock.
c--------------------------------------------------------------------------

            ixx = ieb + (ieblk-1) * leblock_entry

            nindex = array_table(c_nindex, array)
            n      = eblock_table(ixx+c_eblock_blksize-1)

            do i = 1, nindex
               dval1(i) = eblock_table(ixx+c_eblock_range1+i-2)
               dval2(i) = eblock_table(ixx+c_eblock_range2+i-2)
            enddo
 
            block_index_entry(c_blk_index_array) = array
            block_index_entry(c_blk_index_blkno) = ieblk
            block_index_entry(c_blk_index_size)  = n
            block_index_entry(c_blk_index_nind)  = nindex
            do i = 1, nindex
               block_index_entry(c_blk_index_bsegs+i-1) = dval1(i)
               block_index_entry(c_blk_index_esegs+i-1) = dval2(i)
            enddo

            call f_write_disk_ints(list_unit2, index_file_loc2,
     *                  block_index_entry, lblk_index_entry)
            index_file_loc2 = index_file_loc2 + lblk_index_entry 
         enddo

         call flush_eblocks(eblock_table(ieb), n_eblocks,
     *              diskloc_table(idiskloc),
     *              data_unit2, x(ixeblocks), nxeblocks)

c         print *,'Processing of array ',array,' is complete.'
      enddo      ! iarray

c-------------------------------------------------------------------------
c   Close files.
c-------------------------------------------------------------------------

      call f_close_file(list_unit)
      call f_close_file(list_unit2)
      call f_close_file(data_unit)
      call f_close_file(data_unit2)

c--------------------------------------------------------------------------
c   Rename the files.
c--------------------------------------------------------------------------

      call f_renamefile('BLOCK_INDEX2'//char(0),'BLOCK_INDEX'//char(0))
      call f_renamefile('BLOCKDATA2'//char(0),'BLOCKDATA'//char(0))
      return
      end

      subroutine dropmo_expand_process_data_no_mpi_io(array_table, 
     *            narray_table,
     *            index_table, nindex_table,segment_table,
     *            nsegment_table, block_map_table, nblock_map_table,
     *            eindex_table, esegtable, ensegment_table,
     *            scr_dropmo, ndropmo_blksize)
c--------------------------------------------------------------------------
c   Performs the actual processing of the dropmo expansion for each 
c   array in the blockdata_list array.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blockdata_list.h'
      include 'dbugcom.h'
      include 'blkmgr.h'
      include 'dropmo.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table, ensegment_table
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer eindex_table(lindex_table_entry, nindex_table)
      integer esegtable(lsegment_table_entry, ensegment_table)

      integer ndropmo_blksize
      double precision scr_dropmo(ndropmo_blksize)

      integer list_unit, list_unit2, data_unit, data_unit2
      integer*8 data_ptr, data_ptr2, next_diskaddr
      integer ios
      logical is_file_open
 
      integer i, j, n, n_eblocks, max_eblocks, max_eblksize
      integer nindex, array, darray, iarray, ierr, map
      integer memloc
      integer dnindex, dblock
      integer index(mx_array_index), seg(mx_array_index)
      integer dseg(mx_array_index)
      integer index_type(mx_array_index), dindex_type(mx_array_index)
      integer dval1(mx_array_index), dval2(mx_array_index)  
      integer bseg(mx_array_index), eseg(mx_array_index), 
     *        nseg(mx_array_index)
      integer e_bseg(mx_array_index), e_eseg(mx_array_index), 
     *        e_nseg(mx_array_index)
      integer val1(mx_array_index), val2(mx_array_index), iblk,
     *        ieblk, nxeblocks
      integer*8 ixeblocks, i8size, malloced_previously
      integer*8 ieb, idiskloc, ixx

      integer eblock_table(1)
      integer*8 diskloc_table(1)
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
      pointer (iptr, eblock_table)
      pointer (i8ptr, diskloc_table)
#else
      common x, diskloc_table, eblock_table
#endif
      logical restart

#ifdef ALTIX
      dptr = dshptr
      iptr = ishptr
      i8ptr = i8shptr
#endif

c---------------------------------------------------------------------------
c   Open files.
c---------------------------------------------------------------------------

      list_unit = 29
      inquire (unit = list_unit, opened=is_file_open)
      if (is_file_open) then
         rewind list_unit
      else
         open (unit = list_unit, file='BLOCK_INDEX',err = 1000,
     *      iostat=ios, access='SEQUENTIAL', form='UNFORMATTED',
     *      status='OLD')
      endif

      call f_openfile('BLOCKDATA' // char(0), data_unit)
      data_ptr = 0

      list_unit2 = 30
      inquire (unit = list_unit2, opened=is_file_open)
      if (is_file_open) then
         rewind list_unit2
      else
         open (unit = list_unit2, file='BLOCK_INDEX2',err = 1000,
     *      iostat=ios, access='SEQUENTIAL', form='UNFORMATTED',
     *      status='NEW')
      endif

      call f_creatfile('BLOCKDATA2' // char(0), data_unit2)
      data_ptr2 = 0

c-------------------------------------------------------------------------
c   Determine max. number of expanded blocks, max. expanded blocksize.
c-------------------------------------------------------------------------

      max_eblocks = 0
      max_eblksize = 0
      do iarray = 1, narray_list
         array = array_list(iarray)
         nb = array_table(c_numblks, array)

         nindex = array_table(c_nindex, array)
         
         n_eblocks = 1
         do i = 1, nindex
            index(i) = array_table(c_index_array1+i-1,array)
            bseg(i)  = index_table(c_bseg, index(i))
            eseg(i)  = index_table(c_eseg, index(i)) 
            nseg(i)  = index_table(c_nsegments, index(i))

            e_bseg(i) = eindex_table(c_bseg, index(i)) 
            e_eseg(i) = eindex_table(c_eseg, index(i)) 
            e_nseg(i) = eindex_table(c_nsegments, index(i)) 

            n_eblocks = n_eblocks * e_nseg(i)
         enddo

c         print *,'DROPMO_EXPAND_BASIS: array ',array
c         print *,'number of blocks before expansion ',nb
c         print *,'number of blocks after expansion ',n_eblocks
         max_eblocks = max(max_eblocks, n_eblocks) 

c-------------------------------------------------------------------------
c   Loop over all blocks to determine max blocksize.
c-------------------------------------------------------------------------

         do i = 1, nindex
            seg(i) = e_bseg(i)
         enddo
         
         do iblk = 1, n_eblocks
            n = 1
            do i = 1, nindex
               call get_index_segment(index(i), seg(i), esegtable,
     *                             ensegment_table, eindex_table,
     *                             nindex_table, val1, val2)

               n = n * (val2(1) - val1(1) + 1) 
            enddo
      
            max_eblksize = max(max_eblksize, n)

c--------------------------------------------------------------------------
c   Increment segment indices for next block.
c--------------------------------------------------------------------------

            if (i .lt. n_eblocks) then 
               i = 1
  100          continue
               seg(i) = seg(i) + 1 
               if (seg(i) .gt. e_nseg(i)) then
                  seg(i)= e_bseg(i)
                  i = i + 1
                  go to 100
               endif
            endif
         enddo
 
      enddo    ! iarray

c      print *,'max blocks for all arrays ',max_eblocks
c      print *,'max blksize for all eblocks ',max_eblksize

c-------------------------------------------------------------------------
c   Allocate memory.
c-------------------------------------------------------------------------

      call mem_alloc(eblock_table, max_eblocks * leblock_entry, 
     *               intsize, ieb, .true., ierr)

      call mem_alloc(diskloc_table, max_eblocks, 8, idiskloc, 
     *               .true., ierr)
      
      call mem_alloc_query(malloced_previously)
      i8size = maxmem * 1024
      i8size = i8size * 1024
      n = (i8size - malloced_previously) / bytes_per_double
      call mem_alloc(x, n, bytes_per_double, ixeblocks, .true., ierr)
      nxeblocks = n / max_eblksize
 
c-------------------------------------------------------------------------
c   Main processing loop.
c-------------------------------------------------------------------------

      next_diskaddr = 0
      do iarray = 1, narray_list
         array  = array_list(iarray)
         nb     = array_table(c_numblks, array)
         nindex = array_table(c_nindex, array)
         map    = array_table(c_block_map, array)
c         print *,'Processing array ',iarray,' nb ',nb,' nindex ',
c     *         nindex,' indices = ',
c     *   (array_table(c_index_array1+i-1,array),i=1,nindex)

         n_eblocks = 1
         do i = 1, nindex
            index(i) = array_table(c_index_array1+i-1,array)
            index_type(i) = index_table(c_index_type, index(i))
            bseg(i)  = index_table(c_bseg, index(i))
            eseg(i)  = index_table(c_eseg, index(i))
            nseg(i)  = index_table(c_nsegments, index(i))

            e_bseg(i) = eindex_table(c_bseg, index(i))
            e_eseg(i) = eindex_table(c_eseg, index(i))
            e_nseg(i) = eindex_table(c_nsegments, index(i))

            n_eblocks = n_eblocks * e_nseg(i)
         enddo

c         print *,'e_bseg ',(e_bseg(i),i=1,nindex)
c         print *,'e_eseg ',(e_eseg(i),i=1,nindex)
c         print *,'e_nseg ',(e_nseg(i),i=1,nindex)

c-------------------------------------------------------------------------
c   Build the eblock tables.
c-------------------------------------------------------------------------

c         print *,'next_diskaddr for array ',array ,' is ',
c     *       next_diskaddr
         call build_eblock_table(eblock_table(ieb), max_eblocks,
     *         n_eblocks, e_bseg, e_eseg, e_nseg, nindex, 
     *         diskloc_table(idiskloc), index, esegtable, 
     *         ensegment_table, eindex_table, nindex_table,
     *         next_diskaddr)
       
c         print *,'Processing array ',array,' n_eblocks = ',n_eblocks

         restart = .true.
         do iblk = 1, nb

c-------------------------------------------------------------------------
c   Read a descriptor record from the BLOCK_INDEX file.
c-------------------------------------------------------------------------

c            print *,'   Read block ',iblk,' from BLOCKDATA'
            read (list_unit) darray, dblock, dnindex, n,
     *                       (dindex_type(i),i=1,nindex),
     *                       (dseg(i),i=1,nindex),
     *                       (dval1(i),i=1,nindex),
     *                       (dval2(i),i=1,nindex)

c--------------------------------------------------------------------------
c   Make sure the indices match with respect to number and type.
c--------------------------------------------------------------------------

            ierr = 0
            if (dnindex .ne. nindex) ierr = 1
            do i = 1, nindex
               if (dindex_type(i) .ne. index_type(i)) ierr = 1
            enddo

c---------------------------------------------------------------------------
c   Make sure the (dropmo) segment values match the segment values on disk.
c---------------------------------------------------------------------------

            n = 1
            do i = 1, nindex
               seg(i) = block_map_table(c_block_map_seg+i-1,
     *                                  map+dblock-1)
               if (seg(i) .ne. dseg(i)) ierr = 1

               call get_index_segment(index(i), seg(i), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
               if (dval1(i) .ne. val1(i)) ierr = 1
               if (dval2(i) .ne. val2(i)) ierr = 1
               n = n * (val2(i) - val1(i) + 1)
            enddo

            if (ierr .eq. 1) then
               print *,'Error in LIST_TO_BLOCKS: ',
     *               'Index data does not match data ',
     *              'for current job.'
               print *,'nindex ',nindex,' dnindex ',dnindex
               print *,'index_type: ',(index_type(i),i=1,nindex),
     *               ' dindex_type: ',(dindex_type(i),i=1,nindex)
               print *,'seg: ',(seg(i),i=1,nindex),' dseg ',
     *              (dseg(i),i=1,nindex)
               print *,'val1: ',(val1(i),i=1,nindex),' dval1: ',
     *                (dval1(i),i=1,nindex)
               print *,'val2: ',(val2(i),i=1,nindex),' dval2: ',
     *                (dval2(i),i=1,nindex)
               call abort_job()
            endif

c--------------------------------------------------------------------------
c   Read the data block from disk.
c--------------------------------------------------------------------------

            call f_read_disk(data_unit, data_ptr, scr_dropmo, n)
            data_ptr = data_ptr + n

c---------------------------------------------------------------------------
c   Form the contribution of this (dropmo) block into all eblocks to which
c   it contributes.
c---------------------------------------------------------------------------

            call contribute_eblocks(scr_dropmo, val1, val2, nindex,
     *                eblock_table(ieb), n_eblocks, 
     *                diskloc_table(idiskloc), data_unit2, 
     *                x(ixeblocks), nxeblocks, max_eblksize, 
     *                index_type, restart)
            restart = .false.
         enddo   ! iblk

c-------------------------------------------------------------------------
c   Write BLOCK_INDEX2 data for the array and flush out any unwritten
c   eblocks to BLOCKDATA2.
c-------------------------------------------------------------------------

         do ieblk = 1, n_eblocks

c--------------------------------------------------------------------------
c   Write out the BLOCK_INDEX2 data for the eblock.
c--------------------------------------------------------------------------

            ixx = ieb + (ieblk-1) * leblock_entry

            nindex = array_table(c_nindex, array)
            n      = eblock_table(ixx+c_eblock_blksize-1)

            do i = 1, nindex
               dseg(i) = eblock_table(ixx+c_eblock_segs+i-2)
               dval1(i) = eblock_table(ixx+c_eblock_range1+i-2)
               dval2(i) = eblock_table(ixx+c_eblock_range2+i-2)

               index(i) = array_table(c_index_array1+i-1,array) 
               dindex_type(i) = index_table(c_index_type, index(i))
            enddo
 
            write (list_unit2) array, ieblk, nindex, n,
     *                       (dindex_type(i),i=1,nindex),
     *                       (dseg(i),i=1,nindex),
     *                       (dval1(i),i=1,nindex),
     *                       (dval2(i),i=1,nindex)
         enddo

         call flush_eblocks(eblock_table(ieb), n_eblocks,
     *              diskloc_table(idiskloc),
     *              data_unit2, x(ixeblocks), nxeblocks)

c         print *,'Processing of array ',array,' is complete.'
      enddo      ! iarray

c-------------------------------------------------------------------------
c   Close files.
c-------------------------------------------------------------------------

      close (list_unit)
      close (list_unit2)
      call f_close_file(data_unit)
      call f_close_file(data_unit2)

c--------------------------------------------------------------------------
c   Rename the files.
c--------------------------------------------------------------------------

      call f_renamefile('BLOCK_INDEX2'//char(0),'BLOCK_INDEX'//char(0))
      call f_renamefile('BLOCKDATA2'//char(0),'BLOCKDATA'//char(0))
      return

 1000 continue
      print *,'Cannot open BLOCK_INDEX or BLOCKDATA files.'
      print *,'I/O status = ',ios
      call abort_job()
      end

      subroutine build_eblock_table(eblock_table, max_eblocks,
     *         n_eblocks, e_bseg, e_eseg, e_nseg, nindex, 
     *         diskloc_table, index, esegtable, ensegment_table,
     *         eindex_table, nindex_table, diskaddr)
c-------------------------------------------------------------------------
c   Initializes the table entries for the eblock control table.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'

      integer max_eblocks
      integer eblock_table(leblock_entry, max_eblocks)
      integer n_eblocks, nindex
      integer e_bseg(nindex), e_eseg(nindex), e_nseg(nindex)
      integer*8 diskloc_table(max_eblocks)
      integer index(nindex)
      integer ensegment_table
      integer esegtable(lsegment_table_entry,ensegment_table)
      integer nindex_table
      integer eindex_table(lindex_table_entry,nindex_table)

      integer i, j, k, n, seg(nindex)
      integer val1, val2
      integer*8 diskaddr

      do i = 1, nindex
         seg(i) = e_bseg(i)
      enddo

      do i = 1, n_eblocks
         eblock_table(c_eblock_memloc,i) = -1  ! uninitialized
         n = 1
         do j = 1, nindex
            eblock_table(c_eblock_segs+j-1,i) = seg(j)

c-------------------------------------------------------------------------
c   Look up the segment ranges and store them.
c-------------------------------------------------------------------------

            call get_index_segment(index(j), seg(j), esegtable,
     *                       ensegment_table, eindex_table,
     *                       nindex_table, val1, val2)
            eblock_table(c_eblock_range1+j-1,i) = val1
            eblock_table(c_eblock_range2+j-1,i) = val2  
            n = n * (val2 - val1 + 1)
         enddo

         eblock_table(c_eblock_blksize,i) = n
         diskloc_table(i) = diskaddr
c         print *,'   BUILD_EBLOCK: eblock ',i,' size ',n,' diskloc ',
c     *      diskloc_table(i)
         diskaddr = diskaddr + n

c--------------------------------------------------------------------------
c   Increment segment indices for next block.
c--------------------------------------------------------------------------

         if (i .lt. n_eblocks) then   ! no need for last eblock.
            j = 1
  100       continue
            seg(j) = seg(j) + 1
            if (seg(j) .gt. e_eseg(j)) then
               seg(j)= e_bseg(j)
               j = j + 1
               if (j .gt. nindex) then
                  print *,'ERROR: seg range err in build_eblock_table'
                  print *,'j, nindex ',j,nindex
                  print *,'seg ',(seg(k),k=1,nindex)
                  print *,'e_bseg ',(e_bseg(k),k=1,nindex)
                  print *,'e_eseg ',(e_eseg(k),k=1,nindex)
                  call abort_job()
               endif
 
               go to 100
            endif
         endif
      enddo
 
      return
      end

      subroutine contribute_eblocks(data, val1, val2, nindex,
     *                eblock_table, n_eblocks, diskloc_table,
     *                data_unit2, scr, nscr_blks, max_blksize,
     *                index_type, restart)
c--------------------------------------------------------------------------
c   Stores the contribution of one dropmo block into all eblocks to 
c   which it contributes.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'machine_types.h'
      include 'dropmo.h'

      double precision data(*)
      integer nindex
      integer val1(nindex), val2(nindex)
      integer n_eblocks
      integer eblock_table(leblock_entry, n_eblocks)
      integer*8 diskloc_table(n_eblocks)
      integer data_unit2
      integer nscr_blks, max_blksize
      double precision scr(*)
      integer index_type(nindex)
      
      integer i, j, k, l, noverlap, memloc
      integer eval1(nindex), eval2(nindex)
      integer oval1(nindex), oval2(nindex)
      integer ival1(nindex), ival2(nindex)

      integer iblkuse, next_block
      integer ncontrib_blocks
      logical restart
      save next_block

      if (restart) next_block = 1

c-------------------------------------------------------------------------
c   Remap the mo-typeindices into the expanded basis.  Non-mo-indices 
c   are passed through unchanged.
c-------------------------------------------------------------------------

      do i = 1, nindex
         if (index_type(i) .eq. moindex .or. 
     *       index_type(i) .eq. moaindex) then
            ival1(i) = moamap(val1(i))
            ival2(i) = moamap(val2(i))
         else if (index_type(i) .eq. mobindex) then
            ival1(i) = mobmap(val1(i))
            ival2(i) = mobmap(val2(i))
         else
            ival1(i) = val1(i)
            ival2(i) = val1(i)
         endif
      enddo
 
c      print *,'      ival = ',(ival1(j),ival2(j),j=1,nindex)

      ncontrib_blocks = 0
      do i = 1, n_eblocks
         
c-------------------------------------------------------------------------
c   Determine overlap between the current block and the ith eblock.
c-------------------------------------------------------------------------

         noverlap = 0
         do j = 1, nindex
            eval1(j) = eblock_table(c_eblock_range1+j-1,i)
            eval2(j) = eblock_table(c_eblock_range2+j-1,i)
            if (ival1(j) .le. eval2(j) .and. 
     *          ival2(j) .ge. eval1(j)) then
               noverlap = noverlap + 1
               oval1(j) = max(ival1(j), eval1(j))
               oval2(j) = min(ival2(j), eval2(j))
            else
               oval1(j) = 0
               oval2(j) = 0 
            endif
         enddo

         if (noverlap .eq. nindex) then
            ncontrib_blocks = ncontrib_blocks + 1

c-------------------------------------------------------------------------
c   The data block contributes to this eblock. Find the data.
c-------------------------------------------------------------------------

            memloc = eblock_table(c_eblock_memloc,i)
            iblkuse = 0
c            print *,'Block contributes to e_block ',i,' memloc ',
c     *          memloc
            if (memloc .le. 0) then

c--------------------------------------------------------------------------
c   Is there a block that has never been used?
c--------------------------------------------------------------------------
           
c               print *,'next_block, nscr_blks ',next_block,nscr_blks
               if (next_block .le. nscr_blks) then
                  iblkuse = next_block
                  next_block = next_block + 1  
                  go to 100
               endif

c--------------------------------------------------------------------------
c   The eblocks have been exhausted. Flush them to disk. 
c--------------------------------------------------------------------------

               call flush_eblocks(eblock_table, n_eblocks, 
     *              diskloc_table,
     *              data_unit2, scr, nscr_blks)
               iblkuse = 1

  100          continue
               if (iblkuse .eq. 0) then
                  print *,'Error: blocks are full in ',
     *               ' dropmo_expand_basis' 
                  print *,'e_block_table:'
                  do k = 1, n_eblocks
                     print *,'Entry ',k,' : ',(eblock_table(l,k),
     *                 l = 1,leblock_entry)
                  enddo
                  call abort_job()  
               endif
                
c-----------------------------------------------------------------------
c   Set the table values.
c-----------------------------------------------------------------------
 
               memloc = (iblkuse-1)*max_blksize + 1
c               print *,'Move data to memloc ',memloc,' iblkuse ',
c     *            iblkuse 
               if (eblock_table(c_eblock_memloc,i) .lt. 0) then
c                  print *,'Initialize iblkuse ',iblkuse,' to 0.'
c                  print *,'blksize = ',
c     *                  eblock_table(c_eblock_blksize,i)

c---------------------------------------------------------------------------
c   This is a totally new block.  Initialize it to 0's.
c---------------------------------------------------------------------------
 
                  do j = 1, eblock_table(c_eblock_blksize,i)
                     scr(memloc+j-1) = 0
                  enddo 
               else

c--------------------------------------------------------------------------
c   Read the previous contributions in from disk.
c--------------------------------------------------------------------------

c                  print *,'Read iblkuse ',iblkuse,' from disk loc ',
c     *                  diskloc_table(i),' blksize = ',
c     *                  eblock_table(c_eblock_blksize,i)
                  call f_read_disk(data_unit2, diskloc_table(i), 
     *              scr(memloc), eblock_table(c_eblock_blksize,i)) 
               endif
            endif

            eblock_table(c_eblock_memloc,i) = memloc

c--------------------------------------------------------------------------
c   Store into the eblock.
c--------------------------------------------------------------------------

            do j = 1, nindex
               eval1(j) = eblock_table(c_eblock_range1+j-1,i) 
               eval2(j) = eblock_table(c_eblock_range2+j-1,i) 
            enddo

            if (nindex .eq. 2) then
               call eblock2_contribution(
     *            data, val1(1), val2(1), val1(2), val2(2),
     *            scr(memloc),eval1(1),eval2(1),eval1(2),eval2(2),
     *            index_type)
            else if (nindex .eq. 4) then
               call eblock4_contribution(
     *              data, val1(1), val2(1), val1(2), val2(2),
     *                    val1(3), val2(3), val1(4), val2(4),
     *              scr(memloc), eval1(1),eval2(1), eval1(2), eval2(2),
     *                     eval1(3),eval2(3), eval1(4), eval2(4),
     *                     index_type)
            else if (nindex .eq. 6) then
               call eblock6_contribution(
     *              data, val1(1), val2(1), val1(2), val2(2),
     *                    val1(3), val2(3), val1(4), val2(4),
     *                    val1(5), val2(5), val1(6), val2(6),
     *              scr(memloc), eval1(1),eval2(1), eval1(2), eval2(2),
     *                     eval1(3),eval2(3), eval1(4), eval2(4),
     *                     eval1(5),eval2(5), eval1(6), eval2(6),
     *                     index_type)
            else
               print *,'Error: Invalid number of indices in ',
     *            'dropmo_expand_basis.  nindex = ',nindex,
     *            ' should be 2, 4, or 6.'
               call abort_job()
            endif
         endif
      enddo

c      print *,'Data block contributed to ',ncontrib_blocks,' e_blocks.'
      if (ncontrib_blocks .eq. 0) then
         print *,'Non-contributing data block.'
         print *,'data block segment ranges: ',(ival1(j),ival2(j),
     *             j = 1, nindex)
         print *,'eblock_table:'
         do i = 1, n_eblocks
            print *,'Entry ',i,':',(eblock_table(j,i),j=1,
     *            leblock_entry)
         enddo
         call abort_job()
      endif

      return
      end

      subroutine eblock2_contribution(data, a1, a2, b1, b2,
     *                 x, ax1, ax2, bx1, bx2, index_type)
      implicit none
      include 'interpreter.h'
      include 'dropmo.h'

      integer a1, a2, b1, b2, ax1, ax2, bx1, bx2
      double precision data(a1:a2,b1:b2)
      double precision x(ax1:ax2,bx1:bx2)
      integer index_type(*)
 
      integer a, b
      integer aa, bb
 
      do b = b1, b2
         if (index_type(2) .eq. moindex .or.
     *       index_type(2) .eq. moaindex) then
            bb = moamap(b)
         else if (index_type(2) .eq. mobindex) then
            bb = mobmap(b)
         else
            bb = b
         endif
 
      do a = a1, a2
         if (index_type(1) .eq. moindex .or.
     *       index_type(1) .eq. moaindex) then
            aa = moamap(a)
         else if (index_type(1) .eq. mobindex) then
            aa = moamap(a)
         else
            aa = a
         endif

         if (aa .ge. ax1 .and. aa .le. ax2 .and.
     *       bb .ge. bx1 .and. bb .le. bx2)
     *       x(aa,bb) = data(a,b)
      enddo
      enddo

      return
      end

      subroutine eblock4_contribution(data, a1, a2, b1, b2,
     *                    c1, c2, d1, d2,
     *                 x, ax1, ax2, bx1, bx2, 
     *                    cx1, cx2, dx1, dx2,
     *                    index_type)
      implicit none
      include 'interpreter.h'
      include 'dropmo.h'

      integer a1, a2, b1, b2, ax1, ax2, bx1, bx2
      integer c1, c2, d1, d2, cx1, cx2, dx1, dx2
      double precision data(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision x(ax1:ax2,bx1:bx2,cx1:cx2,dx1:dx2)
      integer index_type(*)
 
      integer a, b, c, d
      integer aa, bb, cc, dd
 
      do d = d1, d2
         if (index_type(4) .eq. moindex .or.
     *       index_type(4) .eq. moaindex) then
            dd = moamap(d)
         else if (index_type(4) .eq. mobindex) then
            dd = mobmap(d)
         else
            dd = d
         endif

      do c = c1, c2
         if (index_type(3) .eq. moindex .or.
     *       index_type(3) .eq. moaindex) then
            cc = moamap(c)
         else if (index_type(3) .eq. mobindex) then
            cc = mobmap(c)
         else
            cc = c
         endif
 
      do b = b1, b2
         if (index_type(2) .eq. moindex .or.
     *       index_type(2) .eq. moaindex) then
            bb = moamap(b)
         else if (index_type(2) .eq. mobindex) then
            bb = mobmap(b)
         else
            bb = b
         endif
 
      do a = a1, a2
         if (index_type(1) .eq. moindex .or.
     *       index_type(1) .eq. moaindex) then
            aa = moamap(a)
         else if (index_type(1) .eq. mobindex) then
            aa = moamap(a)
         else
            aa = a
         endif

         if (aa .ge. ax1 .and. aa .le. ax2 .and.
     *       bb .ge. bx1 .and. bb .le. bx2 .and.
     *       cc .ge. cx1 .and. cc .le. cx2 .and.
     *       dd .ge. dx1 .and. dd .le. dx2) 
     *       x(aa,bb,cc,dd) = data(a,b,c,d)
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine eblock6_contribution(data, a1, a2, b1, b2,
     *                    c1, c2, d1, d2, e1, e2, f1, f2,
     *                 x, ax1, ax2, bx1, bx2, 
     *                    cx1, cx2, dx1, dx2, ex1, ex2, fx1, fx2,
     *                    index_type)
      implicit none
      include 'interpreter.h'
      include 'dropmo.h'

      integer a1, a2, b1, b2, ax1, ax2, bx1, bx2
      integer c1, c2, d1, d2, cx1, cx2, dx1, dx2
      integer e1, e2, f1, f2, ex1, ex2, fx1, fx2
      double precision data(a1:a2,b1:b2,c1:c2,d1:d2,e1:e2,f1:f2)
      double precision x(ax1:ax2,bx1:bx2,cx1:cx2,dx1:dx2,
     *                    ex1:ex2,fx1:fx2)
      integer index_type(*)
 
      integer a, b, c, d, e, f
      integer aa, bb, cc, dd, ee, ff
 
      do f = f1, f2
         if (index_type(6) .eq. moindex .or.
     *       index_type(6) .eq. moaindex) then
            ff = moamap(f)
         else if (index_type(6) .eq. mobindex) then
            ff = mobmap(f)
         else
            ff = f
         endif

      do e = e1, e2
         if (index_type(5) .eq. moindex .or.
     *       index_type(5) .eq. moaindex) then
            ee = moamap(e)
         else if (index_type(5) .eq. mobindex) then
            ee = mobmap(e)
         else
            ee = e
         endif
 
      do d = d1, d2
         if (index_type(4) .eq. moindex .or.
     *       index_type(4) .eq. moaindex) then
            dd = moamap(d)
         else if (index_type(4) .eq. mobindex) then
            dd = mobmap(d)
         else
            dd = d
         endif

      do c = c1, c2
         if (index_type(3) .eq. moindex .or.
     *       index_type(3) .eq. moaindex) then
            cc = moamap(c)
         else if (index_type(3) .eq. mobindex) then
            cc = mobmap(c)
         else
            cc = c
         endif
 
      do b = b1, b2
         if (index_type(2) .eq. moindex .or.
     *       index_type(2) .eq. moaindex) then
            bb = moamap(b)
         else if (index_type(2) .eq. mobindex) then
            bb = mobmap(b)
         else
            bb = b
         endif
 
      do a = a1, a2
         if (index_type(1) .eq. moindex .or.
     *       index_type(1) .eq. moaindex) then
            aa = moamap(a)
         else if (index_type(1) .eq. mobindex) then
            aa = moamap(a)
         else
            aa = a
         endif

         if (aa .ge. ax1 .and. aa .le. ax2 .and.
     *       bb .ge. bx1 .and. bb .le. bx2 .and.
     *       cc .ge. cx1 .and. cc .le. cx2 .and.
     *       dd .ge. dx1 .and. dd .le. dx2 .and.
     *       ee .ge. ex1 .and. ee .le. ex2 .and.
     *       ff .ge. fx1 .and. ff .le. fx2)
     *       x(aa,bb,cc,dd,ee,ff) = data(a,b,c,d,e,f)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine flush_eblocks(eblock_table, n_eblocks, diskloc_table,
     *                         data_unit2, scr, nscr_blocks)
c---------------------------------------------------------------------------
c   Writes all eblocks to disk that currently have data in memory.
c--------------------------------------------------------------------------- 
      implicit none
      include 'interpreter.h'
      integer n_eblocks, data_unit2, nscr_blocks
      integer eblock_table(leblock_entry, n_eblocks)
      integer*8 diskloc_table(n_eblocks)
      double precision scr(*)
      double precision sum

      integer ieblk, memloc, i, n
      integer*8 data_ptr2

      do ieblk = 1, n_eblocks
         memloc = eblock_table(c_eblock_memloc,ieblk)
         if (memloc .gt. 0) then
            n = eblock_table(c_eblock_blksize,ieblk)
            data_ptr2 = diskloc_table(ieblk)

c            sum = 0.0
c            do i = 1, n
c               sum = sum + scr(memloc+i-1)*scr(memloc+i-1)
c            enddo

c            print *,'Writing eblock ',ieblk,' to disk: disk loc ',
c     *              data_ptr2,' n ',n, ' sum ',sum 
            call f_write_disk(data_unit2, data_ptr2,
     *                        scr(memloc),n)
 
            eblock_table(c_eblock_memloc,ieblk) = 0
         endif
      enddo
      return
      end

      subroutine get_expanded_basis_segsizes(enalpha_occupied, 
     *     enalpha_virtual, enbeta_occupied, enbeta_virtual)
c----------------------------------------------------------------------------
c   This routine determines the segment sizes in the expanded basis set.
c   It queries the ZMAT parameters.  If the parameters were coded by the
c   user, we do nothing.  If the segment size parameters were not coded
c   in the ZMAT, we determine the correct algorithm to use, and call the 
c   det_segsizes routine to calculate the new values.
c----------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      integer enalpha_occupied, enalpha_virtual
      integer enbeta_occupied, enbeta_virtual
      
      integer ao_segsize, occ_segsize, virt_segsize
      integer nao, nav, nbo, nbv
      integer auto_seg_switch
      character*80 auto_seg_algorithm

      ao_segsize   = 0
      occ_segsize  = 0
      virt_segsize = 0

      call igetparam('SIP_MX_SEGSIZE', 1, ao_segsize)
      call igetparam('SIP_MX_OCC_SEGSIZE', 1, occ_segsize)
      call igetparam('SIP_MX_VIRT_SEGSIZE', 1, virt_segsize)
      if (occ_segsize.lt.1.or.
     &    virt_segsize.lt.1.or.
     &    ao_segsize.lt.1) then
         auto_seg_algorithm = 'SEGMENT_OPTIMIZED'
         call cgetparam('AUTO_SEG_ALGORITHM', 1, auto_seg_algorithm)

         if (auto_seg_algorithm(1:16) .eq. 'MEMORY_OPTIMIZED') then
           auto_seg_switch = 1
         else if (auto_seg_algorithm(1:17) .eq.
     *                                     'SEGMENT_OPTIMIZED') then
           auto_seg_switch = 2
         else
           print *,'====>> Invalid value for AUTO_SEG_ALGORITHM'
           print *,'Value is ',auto_seg_algorithm
           call abort_job()
         endif

c---------------------------------------------------------------------------
c   Save current values of nalpha_occupied, etc., set common block to 
c   desired values.
c---------------------------------------------------------------------------

         nao = nalpha_occupied
         nbo = nbeta_occupied
         nav = nalpha_virtual
         nbv = nbeta_virtual

         nalpha_occupied = enalpha_occupied
         nbeta_occupied  = enbeta_occupied
         nalpha_virtual  = enalpha_virtual
         nbeta_virtual   = enbeta_virtual

c---------------------------------------------------------------------------
c   The new (expanded basis) segment sizes will be stored in int_gen_parms.h
c   by det_segsizes.
c---------------------------------------------------------------------------

         call det_segsizes(auto_seg_switch, my_company_size)

c---------------------------------------------------------------------------
c   Restore the previous values to common.
c---------------------------------------------------------------------------

         nalpha_occupied = nao
         nbeta_occupied  = nbo
         nalpha_virtual  = nav
         nbeta_virtual   = nbv 
      end if

      return
      end
