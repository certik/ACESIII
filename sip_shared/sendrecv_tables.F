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
      subroutine send_tables(proc, request, nreq)
c---------------------------------------------------------------------------
c   Distributes SIP-related table data to processor "proc".
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'sip_tables.h'
      include 'interpreter.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer proc
      integer me, ierr
      integer msg(4)
      integer n, request(*), nreq
      integer ix(1)
      double precision x(1)
#ifdef ALTIX
      pointer (iptr, ix)
      pointer (dptr, x)
#endif
      integer*8 ind, get_index_from_base

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      if (proc .eq. me) return

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   First, send a header message telling how much data is coming.
c---------------------------------------------------------------------------

      msg(1) = nindex_table_sip
      msg(2) = narray_table_sip
      msg(3) = noptable_sip
      msg(4) = nscalar_table_sip
      nreq = nreq + 1
      call mpi_isend(msg, 4, mpi_integer, proc, sip_table_header_tag,
     *               mpi_comm_world, request(nreq), ierr)

c---------------------------------------------------------------------------
c   Now send the data.
c---------------------------------------------------------------------------

      n = nindex_table_sip * lindex_table_entry 
      ind = get_index_from_base(index_table_base, ix, 1)
      nreq = nreq + 1
      call mpi_isend(ix(ind), n, mpi_integer, proc, 
     *               sip_index_table_send_tag, mpi_comm_world,
     *               request(nreq), ierr) 
 
      n = narray_table_sip * larray_table_entry
      ind = get_index_from_base(array_table_base, ix, 1)
      nreq = nreq + 1
      call mpi_isend(ix(ind), n, mpi_integer, proc,
     *               sip_array_table_send_tag, mpi_comm_world,
     *               request(nreq), ierr)

      n = noptable_sip * loptable_entry
      ind = get_index_from_base(optable_base, ix, 1)
      nreq = nreq + 1
      call mpi_isend(ix(ind), n, mpi_integer, proc,
     *               sip_optable_send_tag, mpi_comm_world,
     *               request(nreq), ierr)

      n = nscalar_table_sip
      ind = get_index_from_base(scalar_table_base, x, 2)
      nreq = nreq + 1
      call mpi_isend(x(ind), n, mpi_double_precision, proc,
     *               sip_scalar_table_send_tag, mpi_comm_world,
     *               request(nreq), ierr)

      return
      end

      subroutine recv_tables(proc)
c---------------------------------------------------------------------------
c   Receives SIP-related table data from the master processor.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'sip_tables.h'
      include 'interpreter.h'
      include 'pst.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer proc
      integer me, ierr
      integer msg(4)
      integer i, n, request(4), status(mpi_status_size,4)
      integer ix(1)
      double precision x(1)
#ifdef ALTIX
      pointer (iptr, ix)
      pointer (dptr, x)
#endif
      integer*8 ind, get_index_from_base

      integer create_index_table
      integer create_array_table
      integer create_optable

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      if (proc .eq. me) return

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif
c---------------------------------------------------------------------------
c   First, recv a header message telling how much data is coming.
c---------------------------------------------------------------------------

      call mpi_recv(msg, 4, mpi_integer, proc, sip_table_header_tag,
     *               mpi_comm_world, status, ierr)

c---------------------------------------------------------------------------
c   Create the tables.  Memory for the scalar_table is allocated at the
c   same time as the array_table, so a separate create step is not necessary.
c---------------------------------------------------------------------------
      
      ierr = create_index_table(msg(1))
      if (ierr .ne. 0) then
         print *,'Task ',me,' Error: Cannot create index_table, ',
     *           'nindex_table = ',msg(1)
         call c_exit(1)
      endif
 
      ierr = create_array_table(msg(2))
      if (ierr .ne. 0) then
         print *,'Task ',me,' Error: Cannot create array_table, ',
     *           'narray_table = ',msg(2)
         call c_exit(1)
      endif
 
      ierr = create_optable(msg(3))
      if (ierr .ne. 0) then
         print *,'Task ',me,' Error: Cannot create optable, ',
     *           'narray_table = ',msg(3)
         call c_exit(1)
      endif

c---------------------------------------------------------------------------
c   Now recv the data.
c---------------------------------------------------------------------------

      n = mx_nindex_table * lindex_table_entry 
      ind = get_index_from_base(index_table_base, ix, 1)
      call mpi_irecv(ix(ind), n, mpi_integer, proc, 
     *               sip_index_table_send_tag, mpi_comm_world,
     *               request(1), ierr) 
      nindex_table_sip = msg(1)
 
      n = mx_narray_table * larray_table_entry
      ind = get_index_from_base(array_table_base, ix, 1)
      call mpi_irecv(ix(ind), n, mpi_integer, proc,
     *               sip_array_table_send_tag, mpi_comm_world,
     *               request(2), ierr)
      narray_table_sip = msg(2)

      n = mx_noptable * loptable_entry
      ind = get_index_from_base(optable_base, ix, 1)
      call mpi_irecv(ix(ind), n, mpi_integer, proc,
     *               sip_optable_send_tag, mpi_comm_world,
     *               request(3), ierr)
      noptable_sip = msg(3)

      n = mx_scalar_table
      ind = get_index_from_base(scalar_table_base, x, 2)
      call mpi_irecv(x(ind), n, mpi_double_precision, proc,
     *               sip_scalar_table_send_tag, mpi_comm_world,
     *               request(4), ierr)
      nscalar_table_sip = msg(4)
       
      call mpi_waitall(4, request, status, ierr)  

      return
      end
