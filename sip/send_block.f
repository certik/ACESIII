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
      subroutine send_block(array, array_block, blkndx, dest, tarray,
     *                      tblock, tblkndx, nwblock, comm, opcode, 
     *                      request)
c---------------------------------------------------------------------------
c   Send block "array_block" of "array" to a given destination processor's
c   tarray, tblock.
c   Destination processor rank is relative to communicator "comm".
c---------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer array, array_block, blkndx, dest
      integer tarray, tblock, tblkndx
      integer nwblock
      integer comm
      integer opcode
      integer request

      integer ierr

      request = MPI_REQUEST_NULL
      if (my_company_rank .ne. dest) then
            if (opcode .eq. put_op) then
               call fsumdata(array, array_block, blkndx, dest, 
     *                    tarray, tblock, nwblock, comm, request) 
            else if (opcode .eq. put_replace_op) then
               call fcopydata(array, array_block, blkndx, dest, 
     *                    tarray, tblock, nwblock, comm, request) 
            endif
      endif

      return
      end
