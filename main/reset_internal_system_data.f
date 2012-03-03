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
      subroutine reset_internal_system_data()
c---------------------------------------------------------------------------
c   Resets internal data used by subroutines back to their original state.
c---------------------------------------------------------------------------
      implicit none
      include 'saved_data.h'

      first_master_send_work        = .true.
      first_sip_get_next_work       = .true.
      terminated_worker_termination = .false.
      printwarn_worker_work         = .true.
      flag_internal_event           = .false.
      master_error_processing       = .false.
      master_called                 = .false.
      managers_notified             = .false.
      first_timer                   = .true.
      first_aces_to_erd             = .true.
      first_contract_blocks         = .true.
      first_dcont2                  = .true.
      first_der_int_setup           = .true.
      first_operation_desc          = .true.
      first_blk_event_log           = .true.
      cwork_debug                   = .false.
      first_diis_setup              = .true.
      first_blocks_to_list          = .true.
      first_list_to_blocks          = .true.
      n_saved_requests              = 0
      cwork_iprint                  = 0
      i_been_called_list_to_blocks  = 0
      i_been_called_shmem_list      = 0

c-------------------------------------------------------------------------
c   Remove entries from user_sub table.
c-------------------------------------------------------------------------

      call clear_user_sub()

      return
      end
