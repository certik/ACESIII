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
      subroutine init_prefetch_context()
c----------------------------------------------------------------------------
c   Initializes the prefetch context data structures.
c----------------------------------------------------------------------------
      implicit none
      include 'prefetch.h'
      integer i

      do i = 1, mx_prefetch_context
         ptr_context(i) = 0
      enddo

      nptr_context = 0

      return
      end

      subroutine set_prefetch_context(ind, nind)
c----------------------------------------------------------------------------
c   Saves the array of indices as the "prefetch index context" for the 
c   current loop.
c----------------------------------------------------------------------------
      implicit none
      include 'prefetch.h'
      include 'parallel_info.h'
      include 'trace.h'
      integer nind
      integer ind(nind)
      integer i

      nptr_context = nptr_context + 1
      if (nptr_context .gt. mx_prefetch_context) then
         print *,'Task ',me,' Error: Prefetch context is exhausted.'
         call abort_job() 
      endif

      do i = 1, nind
         context(i,nptr_context) = ind(i)
      enddo

      do i = nind+1,mx_array_index
         context(i,nptr_context) = 0
      enddo

      prefetch_flag = .true.
      return
      end

      logical function is_prefetch_index(ind)
c----------------------------------------------------------------------------
c   Determines whether a given index is in the current prefetch context.
c----------------------------------------------------------------------------
      implicit none
      include 'prefetch.h'
      integer ind
      integer i
  
      is_prefetch_index = .false.
     
      if (nptr_context .eq. 0) return
 
      do i = 1, mx_array_index
         if (ind .eq. context(i,nptr_context)) then
            is_prefetch_index = .true.
            return
         endif 
      enddo
      return
      end

      subroutine unset_prefetch_context()
c----------------------------------------------------------------------------
c   Removes the current prefetch index context, replaces it with the last 
c   previous context.
c----------------------------------------------------------------------------
      implicit none
      include 'prefetch.h'
      include 'parallel_info.h' 
      include 'trace.h'

      nptr_context = nptr_context - 1
      if (nptr_context .lt. 0) then
         print *,'Task ',me,' Error: Prefetch context underflow'
         call abort_job()
      endif

      return
      end
