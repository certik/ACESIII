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
      character*12 function get_operation_desc(opcode)
c------------------------------------------------------------------------------
c   Translate an opcode into a character string.
c------------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'saved_data.h'

      character*12 table(contraction_op:last_opcode)
      integer i
      integer opcode

      save table

      if (first_operation_desc) then
         first_operation_desc = .false.
         table( contraction_op) = contraction_opdesc
         table( sum_op)         = sum_opdesc
         table( reindex_op)     = reindex_opdesc
         table( do_op)          = do_opdesc
         table( enddo_op)       = enddo_opdesc
         table( get_op)         = get_opdesc
         table( user_sub_op)    = user_sub_opdesc
         table( put_op)         = put_opdesc
         table( go_to_op)       = go_to_opdesc
         table( create_op)      = create_opdesc
         table( delete_op)      = delete_opdesc
         table( call_op)        = call_opdesc
         table( return_op)      = return_opdesc
         table( jz_op)          = jz_opdesc
         table( stop_op)        = stop_opdesc
         table( sp_add_op)      = sp_add_opdesc
         table( sp_sub_op)      = sp_sub_opdesc
         table( sp_mult_op)     = sp_mult_opdesc
         table( sp_div_op)      = sp_div_opdesc
         table( sp_equal_op)    = sp_equal_opdesc
         table( sp_nequal_op)   = sp_nequal_opdesc
         table( sp_ge_op)       = sp_ge_opdesc
         table( sp_le_op)       = sp_le_opdesc
         table( sp_gt_op)       = sp_gt_opdesc
         table( sp_lt_op)       = sp_lt_opdesc
         table( sp_ldi_op)      = sp_ldi_opdesc
         table( sp_ldindex_op)  = sp_ldindex_opdesc
         table( pardo_op)       = pardo_opdesc
         table( endpardo_op)    = endpardo_opdesc
         table( exit_op)        = exit_opdesc
         table( assignment_op)  = assignment_opdesc
         table( self_multiply_op) = self_multiply_opdesc
         table( subtract_op)    = subtract_opdesc
         table( collective_sum_op) = collective_sum_opdesc
         table( divide_op)      = divide_opdesc
         table( prepare_op)     = prepare_opdesc
         table( request_op)     = request_opdesc
         table( prequest_op)    = prequest_opdesc
         table( compute_integrals_op) = compute_integrals_opdesc
         table( put_replace_op) = put_replace_opdesc
         table( tensor_op)      = tensor_opdesc
         table( prepare_increment_op) = prepare_increment_opdesc
         table(allocate_op)     = allocate_opdesc
         table(deallocate_op)   = deallocate_opdesc
         table(sp_ldi_sym_op)   = sp_ldi_sym_op_opdesc
      endif

      if (opcode .ge. contraction_op .and.
     *    opcode .le. last_opcode) then
         get_operation_desc = table(opcode)
      else
         get_operation_desc = ' '
      endif

      return
      end

      integer function get_operation_timer_type(opcode)
c------------------------------------------------------------------------------
c   Translate an opcode into a timer type.
c------------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'timerz.h'

      integer table(contraction_op:last_opcode)
      integer i
      integer opcode
      logical first

      save
      data first/.true./

      if (first) then
         first = .false.

         do i = contraction_op, last_opcode
            table(i) = elapsed_time_timer
         enddo

c         table( contraction_op) = cpu_timer
c         table( sum_op)         = cpu_timer
c         table( assignment_op)  = cpu_timer
c         table( subtract_op)    = cpu_timer
c         table( divide_op)      = cpu_timer
c         table( compute_integrals_op) = cpu_timer

      endif

      get_operation_timer_type = table(opcode)
      return
      end

