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
      subroutine init_c(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table, 
     *                      address_table, op)
      implicit none
      include 'interpreter.h'
      include 'proto_events.h'
      include 'mpif.h'
      include 'parallel_info.h'
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
      integer*8 address_table(narray_table)

      integer*8 addr, ixc, get_index_from_base
      integer*8 c_loc64
      integer c, master, ierr
      integer i, n
      integer pst_get_master
      integer status(mpi_status_size)
      logical msg_present
      double precision xxx(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, xxx)
#else
      common xxx
      equivalence (ibuf, xxx)
#endif

      include 'int_gen_parms.h'

      integer company_comm, pst_get_company_comm

c-------------------------------------------------------------------------
c   Locate array address in array_table.
c-------------------------------------------------------------------------

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

      c    = op(c_result_array)
      addr = address_table(c)
      ixc  = get_index_from_base(addr, xxx, 2)
      n    = index_table(c_index_size, 
     *                   array_table(c_index_array1, c))
      master = pst_get_master()

      if (c .eq. 1 .or. c .eq. 2) then
         if (me .eq. master) then
            do i = 1, nscfa
               xxx(ixc+i-1) = xxx(iscfa+i-1) 
            enddo
         endif
      else if (c .eq. 3) then
         if (me .eq. master) then
            do i = 1, nscfb
               xxx(ixc+i-1) = xxx(iscfb+i-1)
            enddo
         endif
      else
         print *,'init_c called with invalid array: ',c
         call abort_job()  
      endif

      company_comm = pst_get_company_comm(me)
      call mpi_bcast(xxx(ixc), nscfa, mpi_double_precision, master,
     *               company_comm, ierr)

c--------------------------------------------------------------------------
c   Reorder the coefficients in each column to correspond to ERD order.
c--------------------------------------------------------------------------

      if (intpkg .eq. flocke_package) 
     *    call apply_erd_order(xxx(ixc), n, ibuf(ierdind), 
     *                        xxx(iscale_fac))
       
      return
      end

      subroutine apply_erd_order(c, n, erd_ind, scale)
c---------------------------------------------------------------------------
c   Reorders the C array column by column into an order compatible with
c   integrals produced by the ERD integral package.
c---------------------------------------------------------------------------
      implicit none
      integer n, erd_ind(n)
      double precision c(n,n)
      double precision scale(n)
      double precision temp(n)

      integer i, j
 
      do j = 1, n
         do i = 1, n
            temp(i) = c(i,j)
         enddo

         do i = 1, n
            c(i,j) = temp(erd_ind(i)) * scale(erd_ind(i))
         enddo

      enddo

      return
      end
