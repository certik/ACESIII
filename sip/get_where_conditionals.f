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
      subroutine get_where_conditionals(optable, start_op, end_op)
c---------------------------------------------------------------------------
c   Add any "where" conditions within the loop structure to the where_table.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'where_table.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer start_op, end_op
      integer optable(loptable_entry, *)
      integer i, j, iwhere_save
 
c----------------------------------------------------------------------------
c   Search the instruction table for WHERE or ENDDO or ENDPARDO.
c   Also check for DO or PARDO to account for nested loops.
c----------------------------------------------------------------------------

      iwhere_save = nwhere 
      do i = start_op+1, end_op
         if (optable(c_opcode, i) .eq. do_op .or.
     *       optable(c_opcode, i) .eq. enddo_op .or.
     *       optable(c_opcode, i) .eq. pardo_op .or.
     *       optable(c_opcode, i) .eq. endpardo_op .or.
     *       optable(c_opcode, i) .eq. do_op) go to 100

         if (optable(c_opcode, i) .eq. where_op) then
            
c-----------------------------------------------------------------------------
c   Add the WHERE condition to the table.
c-----------------------------------------------------------------------------

            if (nwhere .eq. max_where_table) then
               print *,'Error: Too many where conditionals at line ',
     *                 current_line
               print *,'nwhere ',nwhere,' max_where_table ',
     *              max_where_table
               print *,' iwhere ',iwhere 
               print *,'where_table: '
               do j = 1, max_where_table
                  print *,where_table(j,c_where_ind1),
     *                     where_table(j,c_where_cond),
     *                     where_table(j,c_where_ind2)     
               enddo  
     
               call abort_job()
            endif

            nwhere = nwhere + 1
            where_table(nwhere,c_where_ind1) = optable(c_op1_array,i)
            where_table(nwhere,c_where_cond) = optable(c_op2_array,i)
            where_table(nwhere,c_where_ind2) = optable(c_result_array,i)
         endif
      enddo

c----------------------------------------------------------------------------
c   If we have added entries, pick up the beginning index so we can roll 
c   back at the end of the loop.
c----------------------------------------------------------------------------

  100 continue
      if (iwhere_save .ne. nwhere) iwhere = iwhere_save + 1
      return
      end
