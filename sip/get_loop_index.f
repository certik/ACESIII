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
      subroutine get_loop_index(instr_index, index, subrange_flag)
c---------------------------------------------------------------------------
c   Decodes a raw instruction index into the actual loop index and a 
c   logical subrange_flag.  The subrange_flag indicates whether or not 
c   this index is used in a "do jj in j" syntax.
c---------------------------------------------------------------------------

      implicit none
      integer instr_index, index
      logical subrange_flag

      integer in_index_mask, mask_val
      parameter (in_index_mask = 65536)

c----------------------------------------------------------------------------
c  Mask out the subrange_flag with an and.
c----------------------------------------------------------------------------

      mask_val = and(instr_index, in_index_mask)      
      
c----------------------------------------------------------------------------
c   Set the subrange_flag return arg.
c----------------------------------------------------------------------------

      if (mask_val .eq. 0) then
         subrange_flag = .false.
         index         = instr_index
      else
         subrange_flag = .true.
         index         = xor(instr_index, in_index_mask)
      endif
  
      return
      end 
