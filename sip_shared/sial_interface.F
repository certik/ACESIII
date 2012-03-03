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
      integer function create_optable(nentries)
c---------------------------------------------------------------------------
c   Creates an operation table with "nentries" possible entries.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'sip_tables.h'
      include 'machine_types.h'
      include 'dbugcom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
 
      integer nentries
      integer iii(2), size, ierr
      integer*8 ioptable, c_loc64

      call mem_alloc(iii, loptable_entry*nentries, intsize, ioptable, 
     *               .true., ierr)
      create_optable = ierr
      if (ierr .ne. 0) return

#ifdef ALTIX
      optable_base = %loc(ishared_heap(ioptable))
#else
      optable_base = c_loc64(iii, ioptable, intsize)
#endif
      noptable_sip = 0
      mx_noptable = nentries
      if (dbg) print *,'Operation table has a maximum of ',
     *     nentries,' entries.'
      return
      end 

      integer function create_index_table(nentries)
c---------------------------------------------------------------------------
c   Creates an index table with "nentries" possible entries.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'sip_tables.h'
      include 'machine_types.h'
      include 'dbugcom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nentries
      integer iii(2), size, ierr
      integer*8 itable, c_loc64

      call mem_alloc(iii, lindex_table_entry*nentries, intsize, itable,
     *               .true., ierr)
      create_index_table = ierr
      if (ierr .ne. 0) return

#ifdef ALTIX
      index_table_base = %loc(ishared_heap(itable))
#else
      index_table_base = c_loc64(iii, itable, intsize)
#endif
      nindex_table_sip = 0
      mx_nindex_table = nentries

      if (dbg) print *,'Index table has a maximum of ',
     *     nentries,' entries.'
      return
      end

      integer function create_array_table(nentries)
c---------------------------------------------------------------------------
c   Creates an array table with "nentries" possible entries.
c   Also creates a table for scalar values.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'sip_tables.h'
      include 'machine_types.h'
      include 'dbugcom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nentries
      integer iii(2), size, ierr
      double precision xxx(2)
      integer*8 itable
      integer*8 c_loc64

      call mem_alloc(iii, larray_table_entry*nentries, intsize, itable,
     *               .true., ierr)
      create_array_table = ierr
      if (ierr .ne. 0) return

#ifdef ALTIX
      array_table_base = %loc(ishared_heap(itable))
#else
      array_table_base = c_loc64(iii, itable, intsize)
#endif
      narray_table_sip = 0
      mx_narray_table = nentries
      if (dbg) print *,'Array table has a maximum of ',
     *     nentries,' entries.'

      call mem_alloc(xxx, nentries, bytes_per_double, itable,
     *               .true., ierr)
      create_array_table = ierr
      if (ierr .ne. 0) return

#ifdef ALTIX
      scalar_table_base = %loc(dshared_heap(itable))
#else
      scalar_table_base = c_loc64(xxx,itable,bytes_per_double)
#endif
      nscalar_table_sip = 0
      mx_scalar_table = nentries

      return
      end

      integer function add_optable(opcode, op1_array, op2_array,
     *                             result_array, index_array,
     *                             user_sub_index, lineno)
c-------------------------------------------------------------------------
c   Adds an entry to an existing optable.
c-------------------------------------------------------------------------

      implicit none
      include 'sip_tables.h'
      include 'interpreter.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer opcode, op1_array, op2_array, result_array, 
     *        user_sub_index, lineno
      integer index_array(mx_array_index)

      integer optable(1)
#ifdef ALTIX
      pointer (ptr, optable)
#else
      common optable
#endif
      integer*8 get_index_from_base
      integer*8 ioptable, ktable
      integer i

      if (noptable_sip .ge. mx_noptable) then
         add_optable = -1    ! too many entries in optable
         return
      endif

#ifdef ALTIX
      ptr = ishptr    ! point to shared heap 
#endif
      ioptable = get_index_from_base(optable_base, optable, 1)
      noptable_sip = noptable_sip + 1
      ktable   = ioptable + (noptable_sip-1)*loptable_entry

      do i = 1, loptable_entry
         optable(ktable+i-1) = 0    ! set all fields to 0.
      enddo

c--------------------------------------------------------------------------
c   Fill in relevelant data.
c--------------------------------------------------------------------------

      optable(ktable+c_opcode-1)       = opcode
      optable(ktable+c_op1_array-1)    = op1_array
      optable(ktable+c_op2_array-1)    = op2_array
      optable(ktable+c_result_array-1) = result_array
      do i = 1, mx_array_index
         optable(ktable+c_ind1+i-2) = index_array(i)
      enddo
      optable(ktable+c_user_sub-1)     = user_sub_index
      optable(ktable+c_lineno-1)       = lineno

      add_optable = noptable_sip
      return
      end

c      integer function add_index_table(label, index_size, nsegments,
      integer function add_index_table(index_type, bseg, eseg)
c-------------------------------------------------------------------------
c   Adds an entry to an existing optable.
c-------------------------------------------------------------------------

      implicit none
      include 'sip_tables.h'
      include 'interpreter.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer bseg, eseg, nsegments, index_type

c      character*(*) label
      integer index_table(1)
#ifdef ALTIX
      pointer (ptr, index_table)
#else
      common index_table
#endif
      integer*8 get_index_from_base
      integer*8 itable, ktable
      integer i

      if (nindex_table_sip .ge. mx_nindex_table) then
         add_index_table = -1    ! too many entries in index_table
         return
      endif

#ifdef ALTIX
      ptr = ishptr
#endif
      itable = get_index_from_base(index_table_base,
     *                             index_table(1), 1)
      nindex_table_sip = nindex_table_sip + 1
      ktable = itable + (nindex_table_sip-1)*lindex_table_entry

      do i = 1, lindex_table_entry
         index_table(ktable+i-1) = 0    ! set all fields to 0.
      enddo

c--------------------------------------------------------------------------
c   Fill in relevelant data.
c--------------------------------------------------------------------------

      index_table(ktable+c_index_size-1)    = 0    ! filled in later
      index_table(ktable+c_nsegments-1)     = 0    ! filled in later
      index_table(ktable+c_bseg-1)          = bseg 
      index_table(ktable+c_eseg-1)          = eseg 
      index_table(ktable+c_index_type-1)    = index_type 
      index_table(ktable+c_next_seg-1)      = 0
      index_table(ktable+c_subindex_ptr-1)  = 0

      add_index_table = nindex_table_sip
      return
      end 

c      integer function add_array_table(label, nindex, array_type, 
      integer function add_array_table(nindex, array_type, 
     *                                 numblks, index_array, 
     *                                 svalue)
c-------------------------------------------------------------------------
c   Adds an entry to an existing array table
c-------------------------------------------------------------------------

      implicit none
      include 'sip_tables.h'
      include 'interpreter.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nindex, array_type, numblks
      integer index_array(mx_array_index)
      double precision svalue

c      character*(*) label
      integer array_table(1)
      double precision scalar_table(1)
#ifdef ALTIX
      pointer (iptr, array_table)
      pointer (dptr, scalar_table)
#else
      common scalar_table, array_table
#endif
      integer*8 get_index_from_base
      integer*8 itable, ktable
      integer i
      integer*8 c_loc64

      if (narray_table_sip .ge. mx_narray_table) then
         add_array_table = -1    ! too many entries in index_table
         return
      endif

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif
      itable = get_index_from_base(array_table_base, 
     *                             array_table, 1)
      narray_table_sip = narray_table_sip + 1
      ktable = itable + (narray_table_sip-1)*larray_table_entry
      
      do i = 1, larray_table_entry
         array_table(ktable+i-1) = 0    ! set all fields to 0.
      enddo

c--------------------------------------------------------------------------
c   Fill in relevelant data.
c--------------------------------------------------------------------------

      array_table(ktable+c_nindex-1)     = nindex
      array_table(ktable+c_array_type-1) = array_type
      array_table(ktable+c_numblks-1)    = numblks
      do i = 1, mx_array_index
         array_table(ktable+c_index_array1+i-2) = index_array(i)
      enddo

c---------------------------------------------------------------------------
c   If the type is scalar, add it to the scalar table.
c---------------------------------------------------------------------------

      if (array_type .eq. scalar_value) then
         nscalar_table_sip = nscalar_table_sip + 1
         if (nscalar_table_sip .ge. mx_scalar_table) then
            print *,'Error: Max. scalars = ',mx_scalar_table
            add_array_table = -1    ! too many entries in scalar_table
            return
         endif

         array_table(ktable+c_scalar_index-1) = nscalar_table_sip

         itable = get_index_from_base(scalar_table_base,
     *                             scalar_table, 2)
         scalar_table(itable+nscalar_table_sip-1) = svalue
      endif

      add_array_table = narray_table_sip
      return
      end 

      subroutine read_tables(obj_filename)
c----------------------------------------------------------------------------
c   Reads pre-compiled tables needed for SIP from the ".sio" object file.
c   Memory for the index_table, array_table, and operation_table is 
c   allocated, and the data is read into the tables.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'sip_tables.h'
      include 'dbugcom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer index_table(lindex_table_entry,1)
      integer array_table(larray_table_entry,1)
      integer optable(loptable_entry,1)
      double precision scalar_table(1)
#ifdef ALTIX
      pointer (ptr_index, index_table)
      pointer (ptr_array, array_table)
      pointer (ptr_optable, optable)
      pointer (ptr_scalar, scalar_table) 
#else
      common scalar_table, index_table, array_table, optable
#endif

      integer*8 itable
      integer*8 get_index_from_base
      integer i, j, n, fid, fversion, frelease, ierr
      integer nscalar_table_save
      integer create_array_table
      integer create_index_table
      integer create_optable

      character*(*) obj_filename
      character*256 siofile, aces_exe_path

      integer version
      integer release
      integer id
      integer iun, ios
      integer lenval

      logical file_exists
      save version, release, id

      data id/70707/
      data version/1/
      data release/0/

c--------------------------------------------------------------------------
c   Search for the null terminator at the end of the "progname" string.
c--------------------------------------------------------------------------

      n = 1
  100 continue    
      if (obj_filename(n:n) .ne. char(0) .and.
     *    obj_filename(n:n) .ne. ' ') then
         n = n + 1
         go to 100
      endif
     
      n = n - 1

c--------------------------------------------------------------------------
c   Check first for the existence of the file in the current directory.
c--------------------------------------------------------------------------

      inquire(file = obj_filename(1:n), exist = file_exists)
      if (file_exists) then
         siofile = obj_filename(1:n)
      else
 
c--------------------------------------------------------------------------
c   Check for the existence of the file in the directory defined by the
c   ACES_EXE_PATH environment variable.
c--------------------------------------------------------------------------

         call c_getenv('ACES_EXE_PATH'//char(0),
     *            aces_exe_path, lenval,ierr)
         if (ierr .ne. 0) then
            print *,'Error attempting to read ACES_EXE_PATH ',
     *               'env. variable, ierr = ',ierr
            call c_exit(1)
         endif
        
         siofile = aces_exe_path(1:lenval) // '/sio/' // 
     *             obj_filename(1:n)
      endif
c--------------------------------------------------------------------------
c   Open the output file name and read the file header.
c--------------------------------------------------------------------------

      if (dbg) print *,'Reading SIAL object file ',siofile
 
      iun = 20
      open (unit = iun, file=siofile, status = 'OLD', 
     *      form = 'UNFORMATTED', err = 200, iostat = ios)
      go to 300
  200 continue
      if (ios .ne. 0) then
         print *,'File open error: Cannot open object file ',
     *           siofile
         print *,'I/O status = ',ios
         call c_exit(1)
      endif

  300 continue
      read (iun) fid, fversion, frelease, mx_nindex_table,
     *             mx_narray_table, mx_noptable, nscalar_table_sip
      if (fid .ne. id .or.
     *    fversion .ne. version .or.
     *    frelease .ne. release) then
         print *,'Error: Invalid object header data:'
         print *,'       File id, version, release = ',
     *                    fid, fversion, release
         print *,'       Current program id, version, release = ',
     *             id, version, release
         call c_exit(1) 
      endif

      if (dbg) then
         print *,'Object file contains the following tables:'
         print *,'     Index table:     ',mx_nindex_table,' entries.'
         print *,'     Array table:     ',mx_narray_table,' entries.'
         print *,'     Operation table: ',mx_noptable,' entries.'
         print *,'     Scalar table:    ',nscalar_table_sip,
     *                            ' entries.'
      endif

      nscalar_table_save = nscalar_table_sip

c---------------------------------------------------------------------------
c   Index table
c---------------------------------------------------------------------------

      n    = mx_nindex_table
      ierr = create_index_table(mx_nindex_table) 
      nindex_table_sip = n
#ifdef ALTIX
      ptr_index = ishptr
      ptr_array = ishptr
      ptr_optable = ishptr
      ptr_scalar = dshptr
#endif
      itable = get_index_from_base(index_table_base, index_table, 1)
      do i = 1, nindex_table_sip
         read (iun) (index_table(itable+j-1,i),
     *               j = 1, lindex_table_entry)
      enddo

c---------------------------------------------------------------------------
c   Array table
c---------------------------------------------------------------------------

      n    = mx_narray_table
      ierr = create_array_table(n)
      narray_table_sip = n
      itable = get_index_from_base(array_table_base, array_table, 1)
      do i = 1, narray_table_sip
         read (iun) (array_table(itable+j-1,i),
     *               j = 1, larray_table_entry)
      enddo

c---------------------------------------------------------------------------
c   Operation table
c---------------------------------------------------------------------------

      n    = mx_noptable
      ierr = create_optable(n)
      noptable_sip = n
      itable = get_index_from_base(optable_base, optable, 1)
      do i = 1, noptable_sip
         read (iun) (optable(itable+j-1,i),
     *               j = 1, loptable_entry)
      enddo
  
c---------------------------------------------------------------------------
c   Scalar table.  Memory has been allocated in create_array_table, so we
c   don't need to do it here.
c---------------------------------------------------------------------------

      nscalar_table_sip = nscalar_table_save
      itable = get_index_from_base(scalar_table_base, scalar_table, 2)
      read (iun) (scalar_table(itable+i-1),i=1,nscalar_table_sip)

      close (iun)
      return
      end

      subroutine free_sip_tables()
 
      implicit none
      include 'sip_tables.h'

      integer*8 ind, get_index_from_base
      integer ix(1)
      double precision x(1)

#ifdef ALTIX
      return
#endif

      ind = get_index_from_base(index_table_base, ix, 1)
      call c_mem_free(ix(ind))

      ind = get_index_from_base(array_table_base, ix, 1)
      call c_mem_free(ix(ind))

      ind = get_index_from_base(optable_base, ix, 1)
      call c_mem_free(ix(ind))

      ind = get_index_from_base(scalar_table_base, x, 2)
      call c_mem_free(x(ind))
      return
      end

      subroutine write_tables(obj_filename)
c----------------------------------------------------------------------------
c   Writes tables needed for SIP to disk in a binary format.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'sip_tables.h'
      include 'dbugcom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer index_table(lindex_table_entry,1)
      integer array_table(larray_table_entry,1)
      integer optable(loptable_entry,1)
      double precision scalar_table(1)
#ifdef ALTIX
      pointer (ptr_index, index_table)
      pointer (ptr_array, array_table)
      pointer (ptr_optable, optable)
      pointer (ptr_scalar, scalar_table)
#else
      common scalar_table, index_table, array_table, optable
#endif

      integer*8 itable
      integer*8 get_index_from_base
      integer i, j, n

      character*(*) obj_filename

      integer version
      integer release
      integer id
      integer iun, ios
      save version, release, id

      data id/70707/
      data version/1/
      data release/0/

c--------------------------------------------------------------------------
c   Search for the null terminator at the end of the "progname" string.
c--------------------------------------------------------------------------

      n = 1
  100 continue    
      if (obj_filename(n:n) .ne. char(0) .and.
     *    obj_filename(n:n) .ne. ' ') then
         n = n + 1
         go to 100
      endif
     
      n = n - 1

c--------------------------------------------------------------------------
c   Open the output file name and write the file header.
c--------------------------------------------------------------------------

      if (dbg) print *,'Creating SIAL object file ',obj_filename(1:n)
 
      iun = 20
      open (unit = iun, file=obj_filename(1:n), 
     *      form = 'UNFORMATTED', iostat = ios)
      if (ios .ne. 0) then
         print *,'File open error: Cannot open object file ',
     *            obj_filename(1:n)
         print *,'I/O status = ',ios
         call c_exit(1)
      endif

      write (iun) id, version, release, nindex_table_sip,
     *             narray_table_sip, noptable_sip, 
     *             nscalar_table_sip

c---------------------------------------------------------------------------
c   Index table
c---------------------------------------------------------------------------

#ifdef ALTIX
      ptr_index = ishptr
      ptr_array = ishptr
      ptr_optable = ishptr
      ptr_scalar = dshptr
#endif
      itable = get_index_from_base(index_table_base, index_table, 1)
      do i = 1, nindex_table_sip
         write (iun) (index_table(itable+j-1,i),
     *               j = 1, lindex_table_entry)
      enddo

c---------------------------------------------------------------------------
c   Array table
c---------------------------------------------------------------------------

      itable = get_index_from_base(array_table_base, array_table, 1)
      do i = 1, narray_table_sip
         write (iun) (array_table(itable+j-1,i),
     *               j = 1, larray_table_entry)
      enddo

c---------------------------------------------------------------------------
c   Operation table
c---------------------------------------------------------------------------

      itable = get_index_from_base(optable_base, optable, 1)
      do i = 1, noptable_sip
         write (iun) (optable(itable+j-1,i),
     *               j = 1, loptable_entry)
      enddo
  
c---------------------------------------------------------------------------
c   Scalar table
c---------------------------------------------------------------------------

      itable = get_index_from_base(scalar_table_base, scalar_table, 2)
      write (iun) (scalar_table(itable+i-1),i=1,nscalar_table_sip)

      close (iun)
      return
      end
