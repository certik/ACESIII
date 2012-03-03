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
      subroutine set_opblock(array, block, blkndx, op)
c---------------------------------------------------------------------------
c   Handles the setting of the block into the c_opblock field of the
c   instruction table.
c
c   The block_computed_flag is checked first.  If it has already been
c   set, set_opblock does nothing.  If the block_computed_flag has not
c   yet been set, the c_opblock field of the instruction "op" is set
c   to the block number.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'trace.h'

      integer array, block, blkndx, op(loptable_entry)
      integer i
      integer value

      call get_block_computed_flag(array, block, blkndx, value)
      if (value .eq. 0) then
         op(c_opblock) = block
         op(c_opblkndx) = blkndx
      endif

      return
      end
