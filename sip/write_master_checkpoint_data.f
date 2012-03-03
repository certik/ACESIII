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
      subroutine write_master_checkpoint_data(scalar_table, 
     *               nscalar_table, index_table, nindex_table)
c----------------------------------------------------------------------------
c   Writes checkpoint data needed to restart a job to the master checkpoint
c   data file.  
c----------------------------------------------------------------------------

      implicit none
      include 'checkpoint_data.h'
      include 'trace.h'
      include 'interpreter.h'

      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)

      integer i, ios
      integer iop, start_op, end_op

c---------------------------------------------------------------------------
c   Open the master file.
c---------------------------------------------------------------------------

      open(file='master_ckpt_file2', unit=master_ckpt_unit,
     *     ACTION='WRITE',access='SEQUENTIAL', form='UNFORMATTED',
     *     err = 200, iostat = ios)
         
c--------------------------------------------------------------------------
c   Write the data.
c--------------------------------------------------------------------------

      write (master_ckpt_unit) current_line, current_op
      write (master_ckpt_unit) nactive_allocate_table
      write (master_ckpt_unit) nactive_create_table
      write (master_ckpt_unit) nckpt_arrays
      write (master_ckpt_unit) nscalar_table, nindex_table

c---------------------------------------------------------------------------
c   Save the current program execution context.
c---------------------------------------------------------------------------

      call get_program_context(iop, start_op, end_op)
      write (master_ckpt_unit) iop, start_op, end_op

c---------------------------------------------------------------------------
c   Save necessary tables.
c---------------------------------------------------------------------------

      do i = 1, nactive_allocate_table
         write (master_ckpt_unit) active_allocate_table(i),
     *                active_allocate_op(i)
      enddo

      do i = 1, nactive_create_table
         write (master_ckpt_unit) active_create_table(i),
     *               active_create_op(i)
      enddo

      do i = 1, nckpt_arrays
         write (master_ckpt_unit) ckpt_arrays(i), ckpt_diskaddr(i)
      enddo

c----------------------------------------------------------------------------
c   Save the values of the scalar table on the master processor.
c----------------------------------------------------------------------------

      do i = 1, nscalar_table
         write (master_ckpt_unit) scalar_table(i)
      enddo

c---------------------------------------------------------------------------
c   Write the "current_seg" fields from the index_table.
c---------------------------------------------------------------------------

      do i = 1, nindex_table
         write (master_ckpt_unit) index_table(c_current_seg,i)
      enddo

c--------------------------------------------------------------------------
c   Save the instruction stack.
c--------------------------------------------------------------------------

      call checkpoint_instruction_stack(master_ckpt_unit)

c--------------------------------------------------------------------------
c   Close the file.
c--------------------------------------------------------------------------

      close(master_ckpt_unit)

c---------------------------------------------------------------------------
c   Data has been safely written to disk. Rename file to 
c   'master_checkpoint_file'.  
c---------------------------------------------------------------------------

      call f_renamefile('master_ckpt_file2' // char(0),
     *                  'master_ckpt_file' // char(0)) 
      return

  200 continue
      print *,'Error: Cannot open master checkpoint file: iostat = ',
     *          ios
      call abort
      end
