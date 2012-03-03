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
      subroutine handle_cycle(optable, noptable, debug,
     *                        start_op, end_op, iop)
c----------------------------------------------------------------------------
c   Runtime code for the "cycle" instruction: The instruction cycles to the 
c   of the innermost pardo or do loop.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      integer noptable, start_op, end_op, iop
      integer optable(loptable_entry,noptable)
      logical debug

      iop = end_op 
      return
      end
