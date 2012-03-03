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
      subroutine set_program_context(instruction, start_context,
     *                               end_context)
c----------------------------------------------------------------------------
c   Sets the current program instruction execution context.
c----------------------------------------------------------------------------
      implicit none
      include 'context.h'
      integer instruction, start_context, end_context

      iop      = instruction
      start_op = start_context
      end_op   = end_context 
      return
      end

      subroutine get_program_context(instruction, start_context,
     *                               end_context)
c----------------------------------------------------------------------------
c   Returns the current program instruction execution context.
c----------------------------------------------------------------------------
      implicit none
      include 'context.h'
      integer instruction, start_context, end_context

      instruction   = iop
      start_context = start_op
      end_context   = end_op
      return
      end

