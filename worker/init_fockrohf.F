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
      subroutine init_fockrohf(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table, 
     *                      address_table, op)
      implicit none
      include 'int_gen_parms.h'
      include 'interpreter.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'
      include 'epsilon.h'
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

      integer*8 addr, ixeps, get_index_from_base
      integer fock, i, j, n, ind
      integer ierr
      integer master, pst_get_master
      integer status(MPI_STATUS_SIZE)
      double precision xxx(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, xxx)
#else
      common xxx
      equivalence (ibuf, xxx)
#endif

      integer company_comm, pst_get_company_comm

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

c-------------------------------------------------------------------------
c   Locate array address in array_table.
c-------------------------------------------------------------------------

      fock    = op(c_result_array)
      addr = address_table(fock)
      ixeps  = get_index_from_base(addr, xxx, 2)

c--------------------------------------------------------------------------
c   Determine size of epsilon array.
c--------------------------------------------------------------------------

      ind = array_table(c_index_array1,fock)
      n   = index_table(c_index_size,ind)
      master = pst_get_master()

      do i = 1, n*n
         xxx(ixeps+i-1) = 0.
      enddo

      if (iuhf .eq. 2) then

c-------------------------------------------------------------------------
c   ROHF case: Recv FOCK array from data sent by master.
c-------------------------------------------------------------------------

         company_comm = pst_get_company_comm(me)
         if (fock .eq. 17) then
            call mpi_bcast(xxx(ixeps), nscfa, mpi_double_precision, 
     *              master, company_comm, ierr)
         else if (fock .eq. 18) then
            call mpi_bcast(xxx(ixeps), nscfb, mpi_double_precision, 
     *              master, company_comm, ierr)
         else
            print *,'Error: Invalid index type for ROHF FOCK matrix'
            print *,'Index = ',ind,' index_type = ',
     *       index_table(c_index_type,ind)
            call abort_job()
         endif 

c--------------------------------------------------------------------------
c   Reorder the Fock array to correspond to ERD order.
c--------------------------------------------------------------------------

         if (intpkg .eq. flocke_package)
     *       call apply_erd_fockorder(xxx(ixeps), n, ibuf(ierdind),
     *                                xxx(iscale_fac)) 

      endif  

      return
      end

