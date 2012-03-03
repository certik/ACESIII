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
      subroutine init_stack(depth)
c---------------------------------------------------------------------------
c   Sets up a stack for the (start, end) operation pairs that define a 
c   do loop.
c
c   "Depth" is the maximum stack depth.
c---------------------------------------------------------------------------
      implicit none
      include 'machine_types.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      
      integer x(2) 
#ifdef ALTIX
      pointer (iptr, x)
#else
      common x
#endif
      integer depth
      integer size, ierr
      integer*8 ixstack
      integer sp, stack_depth
      common /dostack/ixstack, sp, stack_depth

      stack_depth = depth
      sp = 1

#ifdef ALTIX
      iptr = ishptr
#endif

      call mem_alloc(x, 2*depth, intsize, ixstack, .true., ierr)
      if (ierr .ne. 0) then
         print *,'Error: Unable to allocate memory for do stack.'
         call abort_job()
      endif
      return
      end

      integer function push_do(start, end)
c--------------------------------------------------------------------------
c   Push a pair of start, end values onto the do stack.
c   Stack overflow --> function returns -1, else 0.
c--------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer x(2) 
#ifdef ALTIX
      pointer (iptr, x)
#else 
      common x
#endif
      integer*8 ixstack
      integer sp, stack_depth
      common /dostack/ixstack, sp, stack_depth

      integer*8 ix
      integer start, end
 
#ifdef ALTIX
      iptr = ishptr
#endif

      if (sp .gt. stack_depth) then
         push_do = -1
      else
         ix = ixstack + 2*(sp-1)
         x(ix)   = start   
         x(ix+1) = end
         sp      = sp + 1
         push_do = 0
      endif

c      print *,'Task ',me,' PUSH ix ',ix,' start, end ',start,end
c      call c_flush_stdout()
      return
      end

      integer function pop_do(start, end)
c---------------------------------------------------------------------------
c   Pops a pair of start, end values off the stack, returns 0 if successful,
c   or -1 if stack underflow results.
c---------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer x(2)
#ifdef ALTIX
      pointer (iptr, x)
#else
      common x
#endif
      integer*8 ixstack
      integer sp, stack_depth
      common /dostack/ixstack, sp, stack_depth

      integer start, end
      integer*8 ix

#ifdef ALTIX
      iptr = ishptr
#endif

      if (sp .eq. 1) then
         pop_do = -1
      else
         pop_do = 0
         sp = sp - 1
         ix = ixstack + 2*(sp-1)
         start = x(ix)
         end   = x(ix+1)
      endif

c      print *,'Task ',me,' POP ix ',ix,' start end ',start,end
      return
      end

  
      subroutine checkpoint_instruction_stack(unit)
c---------------------------------------------------------------------------
c   Saves the instruction stack to a checkpoint file.
c---------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer x(2)
#ifdef ALTIX
      pointer (iptr, x)
#else
      common x
#endif
      integer*8 ixstack
      integer sp, stack_depth
      common /dostack/ixstack, sp, stack_depth

      integer i, n, unit

#ifdef ALTIX
      iptr = ishptr
#endif

      write (unit) sp
      if (sp .gt. 1) then
         n = (2*(sp-1))/4
         write (unit) (x(ixstack+i-1),i=1,4*n)
      endif
 
      return
      end

      subroutine restore_instruction_stack(unit)
c---------------------------------------------------------------------------
c   Reads the instruction stack from the checkpoint file.
c---------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer x(2)
#ifdef ALTIX
      pointer (iptr, x)
#else
      common x
#endif
      integer*8 ixstack
      integer sp, stack_depth
      common /dostack/ixstack, sp, stack_depth

      integer i, n, unit

#ifdef ALTIX
      iptr = ishptr
#endif

      read (unit) sp
      if (sp .gt. 1) then
         n = (2*(sp-1)) / 4    ! each entry is startop, endop, iwhere, nwhere
         read (unit) (x(ixstack+i-1),i=1,4*n)
      endif

      return
      end

      subroutine broadcast_instruction_stack(comm)
c---------------------------------------------------------------------------
c   Broadcasts the instruction stack to each process in the communicator.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer x(2)
#ifdef ALTIX
      pointer (iptr, x)
#else
      common x
#endif
      integer*8 ixstack
      integer sp, stack_depth
      common /dostack/ixstack, sp, stack_depth

      integer i, comm, ierr

#ifdef ALTIX
      iptr = ishptr
#endif

      call mpi_bcast(x(ixstack), 2*stack_depth, MPI_INTEGER, 0,
     *               comm, ierr)
      call mpi_bcast(sp, 1, MPI_INTEGER, 0, comm, ierr)
      return
      end

      subroutine restore_do_loops(optable, noptable, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      checkpoint_seg, ckpt_start_op, ckpt_end_op) 
c---------------------------------------------------------------------------
c   Reinitializes each doloop instruction in the instruction stack during 
c   a checkpoint/restart operation.
c
c   The normal do_loop_init routine is first called, then the current_seg 
c   entry in the index_table is set to the corresponding entry from the 
c   checkpoint_seg array.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'parallel_info.h'
      include 'interpreter.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer noptable, nindex_table, nsegment_table
      integer optable(loptable_entry,noptable)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer checkpoint_seg(nindex_table)
      integer ckpt_start_op, ckpt_end_op

      integer x(2)
#ifdef ALTIX
      pointer (iptr, x)
#else
      common x
#endif
      integer*8 ixstack
      integer sp, stack_depth
      common /dostack/ixstack, sp, stack_depth

      integer i, n, comm, ierr
      integer index, start_op, end_op, iop
      logical subrange_flag

#ifdef ALTIX
      iptr = ishptr
#endif

      n = (2*(sp-1))/4   ! number of stack entries.

      do i = 1, n
         iop = x(ixstack+4*(i-1))

         if (optable(c_opcode,iop) .eq. do_op) then
            optable(c_oploop, iop) = 1  ! loop init flag
            call get_loop_index(optable(c_ind1,iop), index, 
     *                  subrange_flag)
         endif
      enddo
      
c------------------------------------------------------------------------
c   Carry out the same procedure for the checkpointed start_op and 
c   end_op.  It is possible that this is also a doloop instruction, and
c   therefore that loop must also be initialized to execute correctly.
c------------------------------------------------------------------------

      iop = ckpt_start_op
      if (optable(c_opcode,iop) .eq. do_op) then
         optable(c_oploop, iop) = 1  ! loop init flag
         call get_loop_index(optable(c_ind1,iop), index,
     *                  subrange_flag)
      endif
 
      return
      end

