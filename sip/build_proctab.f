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
      subroutine build_proctab(optable, noptable, proctab, nproctab)
c---------------------------------------------------------------------------
c   Builds entries for proctab (table used in proc calls.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'

      integer i, j, target, opcode, last
      integer noptable, nproctab
      integer optable(loptable_entry,noptable)
      integer proctab(2,*)

c---------------------------------------------------------------------------
c   First instruction is a jump past the procs.  This is used as the last 
c   range.
c---------------------------------------------------------------------------

      opcode = optable(c_opcode, 1)
      if (opcode .ne. go_to_op) return       ! no procs

      last = optable(c_result_array,1) - 1   ! jump target address.
      if (last .eq. 1) return                ! no procs 

c--------------------------------------------------------------------------
c   Procs are present.  Create an entry for each unique call.
c--------------------------------------------------------------------------
 
      nproctab = 0
      do i = 1, noptable
         if (optable(c_opcode,i) .eq. call_op) then
            target = optable(c_result_array,i)
     
c----------------------------------------------------------------------------
c   Search for previous entry with this target.
c----------------------------------------------------------------------------

            do j = 1, nproctab
               if (proctab(1,j) .eq. target) go to 100
            enddo

c--------------------------------------------------------------------------
c   This is a call to a new proc.  Add it to the proctab.
c--------------------------------------------------------------------------

            nproctab = nproctab + 1
            proctab(1,nproctab) = target
            proctab(2,nproctab)   = last
  100       continue
         endif
      enddo

      return
      end

