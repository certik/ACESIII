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
      subroutine init_eps(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table, 
     *                      op)
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'
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

      integer eps, i, n, ind
      integer ierr
      integer master, pst_get_master
      integer status(mpi_status_size)
      logical msg_present

      include 'int_gen_parms.h'
      include 'epsilon.h'

      integer company_comm, pst_get_company_comm

      double precision xxx(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, xxx)
#else
      common xxx
      equivalence (ibuf, xxx)
#endif

c-------------------------------------------------------------------------
c   Locate array address in array_table.
c-------------------------------------------------------------------------

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

      eps    = op(c_result_array)

c--------------------------------------------------------------------------
c   Determine size of epsilon array.
c--------------------------------------------------------------------------

      ind = array_table(c_index_array1,eps)
      n   = index_table(c_index_size,ind)
      master = pst_get_master()
      if (me .eq. master) then
         do i = 1, nepsa
            epsilon(i) = xxx(iepsa+i-1)
         enddo

         do i = 1, nepsb
            epsilonb(i) = xxx(iepsb+i-1)
         enddo 
      endif

      company_comm = pst_get_company_comm(me)
      call mpi_bcast(epsilon, nepsa, mpi_double_precision, master,
     *               company_comm, ierr)
      if (nepsb .gt. 0) 
     *   call mpi_bcast(epsilonb, nepsb, mpi_double_precision, master,
     *               company_comm, ierr)
      return
      end
