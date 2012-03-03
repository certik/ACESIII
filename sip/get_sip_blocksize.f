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
      subroutine get_sip_blocksize(blocksize, numblocks)
c-------------------------------------------------------------------------
c   Dummy routine to return blocksize, depending on value of 
c   algorithm_flag.
c-------------------------------------------------------------------------

      implicit none
      include 'int_gen_parms.h'

      integer blocksize, numblocks
      integer get_blkmgr_blocksize
      integer get_blkmgr_numblocks

      blocksize = get_blkmgr_blocksize(0)
      numblocks = get_blkmgr_numblocks(0)
      return
      end
  
