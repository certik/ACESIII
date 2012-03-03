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
      subroutine mem_alloc(base, nwords, element_size, ixmem, 
     *                     heap, ierr)
c--------------------------------------------------------------------------
c   Allocates "nwords*element_size" bytes of memory using MPI_ALLOC_MEM,
c   and returns an index "ixmem" relative to the base address.  Upon 
c   return, the memory may be referenced by Fortran code as base(ixmem+i-1),
c   i = 1, 2, 3, ...,nwords.
c
c   Arguments:
c	base		Base array to reference the returned index.
c	nwords		Number of words of memory to allocate.
c	element_size	Number of bytes per element.
c 	ixmem		Returned index used in conjunction with "base" to
c			reference the allocated memory. This argument 
c			should also be declared integer*8.
c       heap            Logical variable.  heap = .true. requests heap  
c                       allocation, .false. requests allocation on 
c                       shared heap. 
c	ierr		Return code.  0 = Successful allocation, otherwise
c			the subroutine was unable to obtain the required
c			amount of memory.
c	
c--------------------------------------------------------------------------

      implicit none
      include 'machine_types.h'
      include 'mem_alloc_data.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer base(*)
      integer nwords, element_size, ierr
      integer*8 ixmem
      integer*8 addr1, addr2
      integer*8 get_index_from_base
      integer*8 size
      integer type
      integer sheap_flag 

      integer*8 ixx, c_loc64

      logical heap

      ierr = 0

      if (nwords .lt. 0) then
         print *,'Error: mem_alloc was called with nwords = ',nwords
         call abort_job()
      endif

      if (element_size .eq. 8) then

c-------------------------------------------------------------------------
c   Force alignment on a 8-byte boundary.
c-------------------------------------------------------------------------

         if (memnxt/8*8 .ne. memnxt) then
            memnxt = (memnxt+7)/8*8 
         endif
      endif

      size = nwords
      size = size * element_size

c      print *,'Task ',me,' MEM_ALLOC: nwords, element_size, size ',
c     *   nwords, element_size, size,' memnxt ',memnxt,
c     *   ' mem_malloced ',mem_malloced
      if (memnxt + size .gt. mem_malloced) then
         print *,'Memory exhausted.'
         print *,'mem_malloced: ',mem_malloced,' memnxt ',memnxt,
     *        ' size requested ',element_size,'*',nwords
         ierr = 1
         return
      endif

c      print *,'Task ',me,' memix, intsize = ',memix,intsize
#ifdef ALTIX
      addr1 = %loc(ishared_heap)
#else
      addr1 = c_loc64(membase, memix, intsize)
#endif
      addr2 = addr1 + memnxt
      if (element_size .eq. 4) then
         type = 1
      else if (element_size .eq. 8) then
         type = 2
      else
         print *,'MEM_ALLOC: invalid element size ',element_size
      endif

c      print *,'Task ',me,' ADDR1, ADDR2 ',ADDR1,ADDR2
#ifdef ALTIX
      ixmem = get_index_from_base(addr2, ishared_heap, type)
#else
      ixmem = get_index_from_base(addr2, base, type) 
#endif
c      print *,'Task ',me,' IXMEM RETURNED FROM MEM_ALLOC: ',
c     *       ixmem,' @base(ixmem) ',
#ifdef ALTIX
c     *   %loc(ishared_heap(ixmem)), %loc(dshared_heap(ixmem))
#else
c     *   c_loc64(base, ixmem, element_size)
#endif
      memnxt = memnxt + size
c      print *,'Task ',me,' New value of memnxt = ',memnxt

      return
      end

      subroutine mem_alloc_init(mem_needed, sheap_flag, ierr)
c------------------------------------------------------------------
c   Performs an actual malloc of the amount of memory needed.  This
c   will only be done one time in the life of an execution, even
c   if we are called multiple times.
c------------------------------------------------------------------
      implicit none
      include 'mem_alloc_data.h'
      include 'machine_types.h'
      include 'parallel_info.h'
c      include 'nothreads.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      
      integer mem_needed
      integer nwords, ierr
      integer*8 c_loc64
      integer sheap_flag
      integer*8 temp

      logical first
      save first

      data first/.true./

      ierr = 0

      if (first) then
         first = .false.

         call init_machine_types()   ! must be initialized.
         temp = mem_needed * 1024
         mem_malloced = temp * 1024

         nwords = mem_malloced / intsize
#ifdef ALTIX
         call fallocate_it(sheap_flag, nwords, ierr )
#else

c--------------------------------------------------------------------------
c   "Normal" memory allocation via malloc.
c-------------------------------------------------------------------------

         call malloc_wrapper(nwords, intsize, sheap_flag,
     *                    membase, memix, ierr)
#endif

         memnxt = 0
c         print *,'MEM_ALLOC_INIT: mem_needed ',mem_needed,
c     *          ' mem_malloced ',mem_malloced
c         print *,'Task ',me,' MEM BASE ADDR :',
#ifdef ALTIX
c     *      ishptr
#else
c     *      c_loc64(membase,memix,intsize)
#endif
c         print *,'Task ',me,
c     *     ' memix at end of mem_alloc_init= ',memix
      endif

      return
      end
      
      subroutine mem_alloc_reset()
c-----------------------------------------------------------------------
c   Resets the "next" pointer for the mem_alloc block.
c----------------------------------------------------------------------
      implicit none
      include 'mem_alloc_data.h'
      
      memnxt = 0
      return
      end

      subroutine mem_alloc_query(nbused)
c--------------------------------------------------------------------------
c   Returns the total bytes used so far.
c--------------------------------------------------------------------------
      implicit none
      include 'mem_alloc_data.h'

      integer*8 nbused
      nbused = memnxt
      return
      end 

      integer*8 function get_mem_base_addr()
      implicit none
c--------------------------------------------------------------------------
c   Returns the address of the processor's heap memory.
c--------------------------------------------------------------------------

      include 'mem_alloc_data.h'
      include 'machine_types.h'
      integer*8 c_loc64

#ifdef ALTIX
      include 'sheap.h'

      get_mem_base_addr = ishptr
#else
      get_mem_base_addr = c_loc64(membase, memix, intsize)
#endif
      return
      end 

      subroutine mem_alloc_free(addr, ierr)
c--------------------------------------------------------------------------
c   Reclaims all memory allocated beyond "addr" in mem_alloc_init.
c   This memory can then be claimed again via mem_alloc.
c--------------------------------------------------------------------------
      include 'mem_alloc_data.h' 
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer*8 addr
      integer ierr

      integer*8 c_loc64, base
       
      ierr = 0 
#ifdef ALTIX
      base = %loc(ishared_heap)
#else
      base = c_loc64(membase, memix, intsize)
#endif
      if (addr .lt. base) then
         print *,'Error in mem_alloc_free: called with addr ',
     *       addr,' but the base address is ',base
         ierr = 1
      endif
      memnxt = addr - base
      return
      end

      integer*8 function get_max_heap_usage()
c--------------------------------------------------------------------------
c   Returns the number of bytes malloced for heap memory.
c--------------------------------------------------------------------------
      implicit none
      include 'mem_alloc_data.h'

      get_max_heap_usage = mem_malloced
      return
      end
