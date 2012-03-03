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
      subroutine memreq_init(size, heap_flag, ierr)
c-------------------------------------------------------------------------
c   Initializes the memory tables, allocates "size" bytes.
c-------------------------------------------------------------------------
      implicit none
      include 'mem.h'
      include 'mpif.h'
      include 'machine_types.h'

      integer ierr
      integer*8 size
      integer nints

      logical heap_flag

      if (size .lt. 0) then
         max_memsize = -1
         return
      endif

      nints = (size+intsize-1)/intsize
      max_memsize = nints * intsize

      call mem_alloc(anchor, nints, intsize, mem_index,
     *                     heap_flag, ierr)
      if (ierr .ne. 0) then
         print *,'Error attempting to allocate ',size,' bytes.'
         ierr = 1
      endif

      call mem_reset_subsystem()
      return
      end

      subroutine mem_reset_subsystem()
c------------------------------------------------------------------------
c  Resets all memory tables.  Any handles and memory objects that have 
c  previously been defined are no longer valid after calling this routine.
c------------------------------------------------------------------------
      implicit none
      include 'mem.h'
      integer i, j

      do i = 1, max_mem_handle
         memtable(i,memflag) = -1
         memtable(i,memtype)  = -1
         perm_offset(i)       = -1
         temp_offset(i)       = -1
      enddo

      do i = 1, max_config_handle
         do j = 1, max_mem_handle
            config_table(i,j) = 0
         enddo
      enddo

      temp_ptr = 0

      mem_context = 0
      return
      end
      
      integer function memreq(nwords, type, flag)
c------------------------------------------------------------------------
c   Makes a request to the memory subsystem for an object of nwords.
c   Flag = 1 indicates the object is to be allocated permanently, 
c   flag = 0 indicates the object is for temporary storage.
c
c   Type = 1 for integer, 2 for real (float), 3 for double precision.
c
c   The function returns a handle to which the object may be referred
c   in subsequent calls to the memory subsystem.  If the request is too
c   big for available memory, a null handle (0) is returned.
c-------------------------------------------------------------------------
      implicit none
      include 'mem.h'

      integer*8 nwords
      integer  type, flag
      integer i
      integer*8 to_bytes, len
      integer*8 type8

      memreq = 0

      if (nwords .le. 0) then
         return
      endif

      if (flag .ne. perm_memory .and. 
     *    flag .ne. temp_memory) then
         print *,'memreq called with invalid flag: ',flag
         return
      endif

      if (type .ne. 1 .and.
     *    type .ne. 2 .and.
     *    type .ne. 3) then
         print *,'Invalid type for memreq call.'
         return
      endif

      type8 = type
      len = to_bytes(nwords, type8)
      if (len .le. max_memsize) then
         
c--------------------------------------------------------------------------
c   Search for an available handle.
c--------------------------------------------------------------------------

         do i = 1, max_mem_handle
            if (memtable(i,memflag) .eq. -1) then
               memtable(i,memflag) = flag
               memtable(i,memsize) = nwords
               memtable(i,memtype) = type

               memreq = i
               return
            endif
         enddo

c--------------------------------------------------------------------------
c   Too many objects for the table.
c--------------------------------------------------------------------------

         print *,'Memory package contains too many objects for table.'
      endif
      
      return
      end

      integer function mem_config(handles, nhandles)
c--------------------------------------------------------------------------    
c   Registers a memory configuration with the memory subsystem.  The 
c   configuration consists of a number of handles returned from earlier
c   calls to "memreq".
c
c   If a configuration is too large to fit available memory, a null handle
c   is returned.  Otherwise a configuration handle is returned, which must
c   be used to set the memory context before using any of the memory 
c   objects defined by the memory handles.
c--------------------------------------------------------------------------
      implicit none
      include 'mem.h'

      integer nhandles
      integer handles(nhandles)
      integer i, j
      integer*8 mem_total
      integer*8 to_bytes 

      mem_config = 0
      if (nhandles .le. 0) return     ! error, return a null handle
      if (max_memsize .le. 0) then
         print *,'Attempt to configure memory before memreq_init.'
         return  ! attempted config before initialization
      endif

c---------------------------------------------------------------------------
c   Find the next free configuration handle.
c---------------------------------------------------------------------------

      do i = 1, max_config_handle
         if (config_table(i,1) .eq. 0) then

c---------------------------------------------------------------------------
c   Check the total memory for the configuration before saving it.
c---------------------------------------------------------------------------

            mem_total = 0
            do j = 1, nhandles
               if (handles(j) .gt. max_mem_handle) then
                  print *,'Invalid memreq handle in configuration'
                  return
               endif

               mem_total = mem_total + 
     *                to_bytes(memtable(handles(j), memsize), 
     *                         memtable(handles(j), memtype))
            enddo

            if (mem_total .gt. max_memsize) return  ! null handle

c--------------------------------------------------------------------------
c   Save the handles in the configuration.
c--------------------------------------------------------------------------

            do j = 1, nhandles
               config_table(i,j) = handles(j)
            enddo
            
            mem_config = i
            return
         endif
      enddo

      print *,'Max. number of memory configurations has been exceeded.'
      return
      end

      subroutine set_memory_context(config, ierr)
c---------------------------------------------------------------------------
c   Sets up a memory context using the memory objects referenced through
c   a configuration handle.
c
c   Once the context has been set, the memory indices for each memory 
c   object may be retrieved using the function "get_mem".
c
c   Only one memory context may be used at any given time.
c   When the memory configuration needs to be changed, call set_memory_context
c   again with a different configuration, and retrieve the indices for 
c   the new context. 
c 
c   Permanent memory is preserved across context changes.
c
c   ierr = 0 --> All is well, ierr = non-zero --> error.
c---------------------------------------------------------------------------
      implicit none
      include 'mem.h'
      
      integer config, ierr
      integer*8 len, ptr
      integer i, handle
      integer*8 to_bytes

      ierr = 0
      if (config .le. 0 .or.
     *    config .gt. max_config_handle) then
         print *,'Attempt to set an invalid configuration handle: ',
     *            config
         ierr = 1
         return
      endif

      if (config_table(config, 1) .eq. 0) then
         print *,'Attempt to set an undefined configuration handle: ',
     *            config
         ierr = 1
         return
      endif

      mem_context = config

c-------------------------------------------------------------------------
c   Set the byte offsets for the current context.
c   Set up all permanent memory first.
c-------------------------------------------------------------------------

      do i = 1, max_mem_handle
         handle = config_table(config, i)
         if (handle .eq. 0) go to 100     ! no more handles in config.
  
         if (memtable(handle, memflag) .eq. perm_memory .and.
     *       perm_offset(handle) .eq. -1) then
            len = to_bytes(memtable(handle, memsize),
     *                     memtable(handle, memtype))
            perm_offset(handle) = temp_ptr
            temp_ptr = temp_ptr + len
         endif 
      enddo

  100 continue

c-------------------------------------------------------------------------
c   All temp objects are reinitialized by the context change.
c-------------------------------------------------------------------------

      do i = 1, max_mem_handle
         temp_offset(i) = -1
      enddo

c-------------------------------------------------------------------------
c   Now set up all temp memory objects.
c-------------------------------------------------------------------------

      ptr = temp_ptr
      do i = 1, max_mem_handle
         handle = config_table(config, i)
         if (handle .eq. 0) go to 200     ! no more handles in config.

         if (memtable(handle, memflag) .eq. temp_memory) then
            temp_offset(handle) = ptr
            len = to_bytes(memtable(handle, memsize),
     *                     memtable(handle, memtype))
            ptr = ptr + len
         endif
      enddo

  200 continue

      return
      end
 
      integer function get_memory_context()
c------------------------------------------------------------------------
c   Returns the current memory context handle.
c------------------------------------------------------------------------
      implicit none
      include 'mem.h'

      get_memory_context = mem_context
      return
      end

      integer*8 function get_mem(handle, base)
c-------------------------------------------------------------------------
c   Returns the index (relative to base) of the memory for the "handle" 
c   object.
c
c   The "base" must be declared the same type as was used in the memreq
c   call for the object.
c
c   After calling get_mem, the memory is referenced as base(imem+i-1),
c   where imem is the index returned by get_mem, and i goes from 1 to the
c   number of words declared in the memreq call.
c-------------------------------------------------------------------------
      implicit none
      include 'mem.h'
      include 'mpif.h'
      integer me

      integer ierr, handle
      integer*8 base
      integer*8 index, offset 
      integer type

      if (handle .le. 0) then
         call mpi_comm_rank(mpi_comm_world, me, ierr)
         print *,'Task ',me,' GET_MEM was called with a handle of ',
     *            handle
         call mpi_abort(mpi_comm_world, 1, ierr)
      endif

c-------------------------------------------------------------------------
c   Find the memory offset.
c-------------------------------------------------------------------------

      if (memtable(handle, memflag) .eq. perm_memory) then
         offset = perm_offset(handle)
         if (offset .eq. -1) then 
            call mpi_comm_rank(mpi_comm_world, me, ierr)
            print *,'Task ',me,
     *              ' Perm. memory offset is unset for handle ',handle
            index = -10101
            call dump_mem_tables()
            call mpi_abort(mpi_comm_world, 1, ierr)
         endif
      else  
         offset = temp_offset(handle)
         if (offset .eq. -1) then
            call mpi_comm_rank(mpi_comm_world, me, ierr)
            print *,'Task ',me,
     *               ' Temp. memory offset is unset for handle ',
     *               handle
            index = -20202
            call dump_mem_tables()
            call mpi_abort(mpi_comm_world, 1, ierr)
         endif
      endif

      type = memtable(handle, memtype)
      call convert_to_base_index(anchor(mem_index), offset, base,
     *                      type, index)
      get_mem = index 
      return
      end

      integer*8 function get_mem_size(handle, bytes)
c-------------------------------------------------------------------------
c   Returns the size (in bytes) of the memory object referenced by "handle". 
c   Bytes = .true. ---> return size in bytes, else size is returned in 
c   words.
c   A negative size indicates an error.
c-------------------------------------------------------------------------
      implicit none
      include 'mem.h'
      
      integer ierr, handle
      logical bytes
      integer*8 to_bytes

      get_mem_size = -1
      if (handle .le. 0) return

      if (bytes) then
         get_mem_size = to_bytes(memtable(handle,memsize),
     *                           memtable(handle,memtype))
      else
         get_mem_size = memtable(handle,memsize)
      endif
      return
      end

      integer*8 function mem_query()
c---------------------------------------------------------------------------
c   Returns the maximum amount of temp. memory left (in bytes) 
c   in the current memory context.
c
c   A return value of -1 means the memory subsystem has not been initialized.
c---------------------------------------------------------------------------
      implicit none
      include 'mem.h'

      integer i, handle
      integer*8 len, mem_left

      integer*8 to_bytes

      if (max_memsize .lt. 0) then
         mem_query = -1
         return
      endif

      mem_left = max_memsize
      if (mem_context .eq. 0) then
         mem_query = mem_left
         return
      endif

      do i = 1, max_mem_handle
         handle = config_table(mem_context,i)
         if (handle .ne. 0) then

c--------------------------------------------------------------------------
c   We have a handle from the current configuration.
c   Get the amount of memory used from the memtable.
c--------------------------------------------------------------------------

            if (memtable(handle, memflag) .eq. temp_memory) then
               len = to_bytes(memtable(handle, memsize),
     *                        memtable(handle, memtype))
               mem_left = mem_left - len
            endif
         endif
      enddo

c---------------------------------------------------------------------------
c   Subtract off any permanent memory.
c---------------------------------------------------------------------------

      do i = 1, max_mem_handle
         if (perm_offset(i) .ne. -1) then
            len = to_bytes(memtable(i, memsize), 
     *                     memtable(i, memtype))
            mem_left = mem_left - len
         endif
      enddo

      mem_query = mem_left
      return
      end

      subroutine dump_mem_tables()
c--------------------------------------------------------------------------
c   Prints the current state of the mem tables.
c--------------------------------------------------------------------------

      implicit none
      include 'mem.h'
      integer i,j

      print *,'Current handles:'
      do i = 1, max_mem_handle
         if (memtable(i,1) .gt. 0) then
            print *,'   Handle ',i,' size: ',memtable(i,memsize),
     *     ' type: ',memtable(i,memtype),' flag: ',memtable(i,memflag),
     *     ' perm_offset: ',perm_offset(i),' temp_offset: ',
     *       temp_offset(i)
         endif
      enddo

      print *,'Active memory context is ',mem_context
      print *,'Configuration table: '
      do i = 1, max_config_handle
         if (config_table(i,1) .gt. 0) then
            print *,'Configuration ',i,' handles: ',(config_table(i,j),
     *               j = 1, max_mem_handle)
         endif
      enddo
      return
      end
