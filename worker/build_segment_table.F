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
      subroutine build_segment_table(index_table, nindex_table, 
     *                        end_nfps, nshells, mx_ao_seg_size,
     *                        mx_occ_seg_size, mx_virt_seg_size,
     *                        sip_sub_segsize, sip_sub_occ_segsize,
     *                        sip_sub_virt_segsize, sip_sub_ao_segsize,
     *                        nocc_orbitals, naocc_orbitals,
     *                        nbocc_orbitals, 
     *                        nvirt_orbitals, navirt_orbitals,
     *                        nbvirt_orbitals, 
     *                        isegbase, nsegment_table,
     *                        stack_blocksizes, nstacks,
     *                        array_table, narray_table, 
     *                        atom_based, atoms)
c--------------------------------------------------------------------------
c   Builds the segment table from the indices defined in the index_table 
c   and the basis function shell definitions.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'machine_types.h'
      include 'symbolic_constants.h'
c      include 'int_gen_parms.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      
      integer nindex_table, nshells, nsegment_table
      integer nocc_orbitals, naocc_orbitals, nbocc_orbitals
      integer nvirt_orbitals, navirt_orbitals, nbvirt_orbitals 
      integer index_table(lindex_table_entry,nindex_table)
      integer end_nfps(nshells)
      integer*8 isegbase
      integer*8 c_loc64
      integer nstacks
      integer stack_blocksizes(nstacks)
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer atoms(nshells)

      integer segtable(2)
#ifdef ALTIX
      pointer (iptr, segtable)
#endif
      integer naosegs
      integer aosegs(nshells)
      integer norbitals, nsegs, i, j, k, l, type, size, next
      integer bseg, eseg, bseg_table, eseg_table
      integer brange, erange, ierr
      integer minrange, maxrange
      integer*8 itable, ktable
      integer boccval, eoccval, bvirtval, evirtval
      integer baoccval, eaoccval, bavirtval, eavirtval
      integer bboccval, eboccval, bbvirtval, ebvirtval
      integer mx_ao_seg_size, mx_occ_seg_size, mx_virt_seg_size
      integer sip_sub_segsize, sip_sub_occ_segsize,
     *        sip_sub_virt_segsize, sip_sub_ao_segsize
      integer superindex, superindex_type
      integer factor, factore
      integer segsize, range, nsubseg, norb_subseg, nrem, brange_s,
     *        erange_s, nseg_total, subsegment_counter
      integer get_subindex_segsize

      integer translate_symbolic_constant
      logical atom_based

#ifdef ALTIX
      iptr = ishptr
#endif

c---------------------------------------------------------------------------
c   Map the ao shells defined by basis function data into segments.
c--------------------------------------------------------------------------- 

c     norbitals = end_nfps(nshells)
      norbitals = naocc_orbitals + navirt_orbitals   
      call shells_to_segments(end_nfps, nshells, mx_ao_seg_size,
     *                        aosegs, naosegs, atom_based, atoms)

c----------------------------------------------------------------------------
c   Translate the symbolic constants needed for index segment checking.
c----------------------------------------------------------------------------

      boccval = translate_symbolic_constant(bocc, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      eoccval = translate_symbolic_constant(eocc, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      baoccval = translate_symbolic_constant(baocc, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      eaoccval = translate_symbolic_constant(eaocc, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      bboccval = translate_symbolic_constant(bbocc, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      eboccval = translate_symbolic_constant(ebocc, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      bvirtval = translate_symbolic_constant(bvirt, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      evirtval = translate_symbolic_constant(evirt, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      bavirtval = translate_symbolic_constant(bavirt, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      eavirtval = translate_symbolic_constant(eavirt, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      bbvirtval = translate_symbolic_constant(bbvirt, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
      ebvirtval = translate_symbolic_constant(ebvirt, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)

c----------------------------------------------------------------------------
c   Store some values in the "symbolic_constants" table.
c-----------------------------------------------------------------------------

      symbolic_constant_table(1) = naosegs               ! norb
      symbolic_constant_table(2) = eoccval-boccval+1     ! nocc
      symbolic_constant_table(3) = evirtval-bvirtval+1   ! nvirt
      symbolic_constant_table(4) = boccval               ! bocc
      symbolic_constant_table(5) = eoccval               ! eocc
      symbolic_constant_table(6) = bvirtval              ! bvirt
      symbolic_constant_table(7) = evirtval              ! evirt
      symbolic_constant_table(8) = eaoccval-baoccval+1   ! naocc
      symbolic_constant_table(9) = eboccval-bboccval+1   ! nbocc
      symbolic_constant_table(10) = eavirtval-bavirtval+1   ! navirt
      symbolic_constant_table(11) = ebvirtval-bbvirtval+1   ! nbvirt
      symbolic_constant_table(12) = baoccval                ! baocc
      symbolic_constant_table(13) = bboccval                ! bbocc
      symbolic_constant_table(14) = eaoccval                ! eaocc
      symbolic_constant_table(15) = eboccval                ! ebocc
      symbolic_constant_table(16) = bavirtval               ! bavirt
      symbolic_constant_table(17) = bbvirtval               ! bbvirt
      symbolic_constant_table(18) = eavirtval               ! eavirt
      symbolic_constant_table(19) = ebvirtval               ! ebvirt
      symbolic_constant_table(20) = nocc_orbitals           ! noccorb
      symbolic_constant_table(21) = nvirt_orbitals          ! nvirtorb
      symbolic_constant_table(22) = 1                       ! boccorb
      symbolic_constant_table(23) = nocc_orbitals           ! eoccorb
      symbolic_constant_table(24) = nocc_orbitals+1         ! bvirtorb
      symbolic_constant_table(25) = norbitals               ! evirtorb
      symbolic_constant_table(26) = naocc_orbitals          ! naoccorb
      symbolic_constant_table(27) = nbocc_orbitals          ! nboccorb
      symbolic_constant_table(28) = navirt_orbitals         ! navirtorb
      symbolic_constant_table(29) = nbvirt_orbitals         ! nbvirtorb
      symbolic_constant_table(30) = 1                       ! baoccorb
      symbolic_constant_table(31) = 1                       ! bboccorb
      symbolic_constant_table(32) = naocc_orbitals          ! eaoccorb
      symbolic_constant_table(33) = nbocc_orbitals          ! eboccorb
      symbolic_constant_table(34) = naocc_orbitals+1        ! bavirtorb
      symbolic_constant_table(35) = nbocc_orbitals+1        ! bbvirtorb
      symbolic_constant_table(36) = norbitals               ! eavirtorb
      symbolic_constant_table(37) = norbitals               ! ebvirtorb

c----------------------------------------------------------------------------
c   Make a pass through the index table to count the number of segments needed.
c----------------------------------------------------------------------------

      nsegs = 0
      do 100 i = 1, nindex_table

c------------------------------------------------------------------------------
c   Set the superindex in the c_subindex_ptr field if this is a subindex type.
c------------------------------------------------------------------------------

         if (index_table(c_index_type,i) .eq. subindex) then

c------------------------------------------------------------------------------
c   The compiler has filled in the superindex in the bseg field of the table
c   as it is written to disk.  We must store it in the subindex_ptr field.
c------------------------------------------------------------------------------

            superindex = index_table(c_bseg, i) + 1 ! convert to a Fortran index
            superindex_type = index_table(c_index_type,superindex)
            index_table(c_subindex_ptr,i) = superindex

c------------------------------------------------------------------------------
c   Fill in the ranges from the superindex into the index_table entry for the 
c   subindex.
c------------------------------------------------------------------------------

            index_table(c_bseg,i) = index_table(c_bseg,superindex)
            index_table(c_eseg,i) = index_table(c_eseg,superindex)

c----------------------------------------------------------------------------
c   Compiler puts "norb" in table for some moindex indices.
c   Change the symbolic "norb" value to one that is unambiguous.
c----------------------------------------------------------------------------

            if (eseg_table .eq. norb) then
               if (superindex_type .eq. moindex) then
                  eseg_table = evirt
               else if (superindex_type .eq. moaindex) then
                  eseg_table = eavirt
               else if (superindex_type .eq. mobindex) then
                  eseg_table = ebvirt
               endif
            endif
         else    ! not a subindex
            index_table(c_subindex_ptr,i) = 0
         endif

         bseg_table = index_table(c_bseg,i)
         eseg_table = index_table(c_eseg,i)

c----------------------------------------------------------------------------
c   Compiler puts "norb" in table for some moindex indices.
c   Change the symbolic "norb" value to one that is unambiguous.
c----------------------------------------------------------------------------

         if (eseg_table .eq. norb) then
            if (index_table(c_index_type,i) .eq. moindex) then
               eseg_table = evirt
            else if (index_table(c_index_type,i) .eq. moaindex) then
               eseg_table = eavirt
            else if (index_table(c_index_type,i) .eq. mobindex) then
               eseg_table = ebvirt
            endif
         endif

c----------------------------------------------------------------------------
c   If the table segment is a symbolic segment value, translate it now.
c----------------------------------------------------------------------------
 
         if (bseg_table .le. 0) then
            bseg = translate_symbolic_constant(bseg_table, naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals,
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
         else
            bseg = bseg_table
            call check_index_range(bseg, index_table(c_index_type,i),
     *               i, naosegs, boccval, eoccval, bvirtval, evirtval,
     *               baoccval, eaoccval, bavirtval, eavirtval,
     *               bboccval, eboccval, bbvirtval, ebvirtval)
         endif

         if (eseg_table .le. 0) then
            eseg = translate_symbolic_constant(eseg_table,  naosegs,
     *                       norbitals, nocc_orbitals, naocc_orbitals, 
     *                       nbocc_orbitals, mx_ao_seg_size,
     *                       mx_occ_seg_size, mx_virt_seg_size)
         else
            eseg = eseg_table
            call check_index_range(eseg, index_table(c_index_type,i),
     *               i, naosegs, boccval, eoccval, bvirtval, evirtval,
     *               baoccval, eaoccval, bavirtval, eavirtval,
     *               bboccval, eboccval, bbvirtval, ebvirtval)
         endif

         if (index_table(c_index_type,i) .eq. subindex) then

c---------------------------------------------------------------------------
c   Loop over the segments of the superindex.
c---------------------------------------------------------------------------

            bseg = index_table(c_bseg,superindex)
            eseg = index_table(c_eseg,superindex)

            nseg_total = 0
            do j = bseg, eseg

c---------------------------------------------------------------------------
c   Determine the subindex segsize for this segment.
c---------------------------------------------------------------------------

               segsize = get_subindex_segsize(superindex_type, j, 
     *                      sip_sub_segsize,
     *                      sip_sub_ao_segsize, sip_sub_occ_segsize,
     *                      sip_sub_virt_segsize, boccval, eoccval,
     *                      baoccval, eaoccval, bboccval, eboccval)
            
c----------------------------------------------------------------------------
c   Find the range of orbitals in the superindex segment.
c----------------------------------------------------------------------------

               call segment_get_range(superindex_type, j, 
     *                       aosegs, naosegs,
     *                       norbitals, nocc_orbitals,
     *                       naocc_orbitals, nbocc_orbitals,
     *                       mx_ao_seg_size, mx_occ_seg_size,
     *                       mx_virt_seg_size,
     *                       brange, erange)

c----------------------------------------------------------------------------
c   Calculate the number of subindex segments for this superindex segment.
c----------------------------------------------------------------------------

               range   = erange - brange + 1 
               nsubseg = range / segsize 
               if (nsubseg * segsize .ne. range) nsubseg = nsubseg + 1 
               nseg_total = nseg_total + nsubseg
            enddo

            bseg = 1
            eseg = nseg_total
         endif

c---------------------------------------------------------------------------
c   Set remaining fields of index_table entry.
c---------------------------------------------------------------------------

         index_table(c_bseg,i) = bseg
         index_table(c_eseg,i) = eseg
         index_table(c_nsegments,i) = eseg - bseg + 1

c---------------------------------------------------------------------------
c   Count the number of segments.
c---------------------------------------------------------------------------

         nsegs = nsegs + eseg - bseg + 1
  100 continue

c--------------------------------------------------------------------------
c   Allocate the memory for the segment_table.
c--------------------------------------------------------------------------

      call mem_alloc(segtable, lsegment_table_entry*nsegs, intsize, 
     *               itable, .true., ierr)
      if (ierr .ne. 0) then
         print *,'Memory allocation error: ',
     *     ' Cannot allocate segment_table. nsegs = ', nsegs
         call abort_job()
      endif

      isegbase = c_loc64(segtable, itable, intsize)

c--------------------------------------------------------------------------
c   Make a second pass through the index_table, creating segment_table 
c   entries for each segment of each regular index.
c--------------------------------------------------------------------------

      nsegment_table = 0
      do 200 i = 1, nindex_table
         bseg = index_table(c_bseg,i)
         eseg = index_table(c_eseg,i)
         type = index_table(c_index_type,i)
         if (type .eq. subindex) go to 200 ! do normal index entries first 

         minrange = 1000000
         maxrange = 0
         do j = bseg, eseg
            call segment_get_range(type, j, aosegs, naosegs,
     *                       norbitals, nocc_orbitals, 
     *                       naocc_orbitals, nbocc_orbitals, 
     *                       mx_ao_seg_size, mx_occ_seg_size,
     *                       mx_virt_seg_size, 
     *                       brange, erange)
            ktable = itable + nsegment_table * lsegment_table_entry
            segtable(ktable + c_index - 1)   = i
            segtable(ktable + c_segment - 1) = j
            segtable(ktable + c_range1 - 1)  = brange
            segtable(ktable + c_range2 - 1)  = erange
            segtable(ktable + c_subseg1 - 1) = 0
            segtable(ktable + c_subseg2 - 1) = 0
            if (brange .gt. erange) then
               print *,'Error: Index ',i,' segment ',j,' has ',
     *            'brange =',brange,' erange = ',erange,
     *            ' index type ',type
               print *,'Index table:'
               do k = 1, nindex_table
                  print *,'Entry ',k,': ',
     *               (index_table(l,k),l=1,lindex_table_entry)
               enddo
               print *,'norbitals ',norbitals,' nocc_orbitals ',
     *               nocc_orbitals,' naocc_orbitals ',
     *               naocc_orbitals,' nbocc_orbitals ',
     *               nbocc_orbitals
               call abort_job()
            endif
            nsegment_table = nsegment_table + 1
            minrange = min(brange, minrange)
            maxrange = max(erange, maxrange)
         enddo

         index_table(c_index_size,i) = maxrange - minrange + 1
  200 continue

c------------------------------------------------------------------------------
c   Make another pass through the index table, creating the segtable entries
c   for each subindex.
c------------------------------------------------------------------------------

      do 300 i = 1, nindex_table
         type = index_table(c_index_type,i)
         if (type .ne. subindex) go to 300

c-----------------------------------------------------------------------------
c   Determine the superindex and its type, use the superindex size as the
c   subindex size.
c-----------------------------------------------------------------------------

         superindex = index_table(c_subindex_ptr,i)
         type       = index_table(c_index_type, superindex)
         index_table(c_index_size,i) = 
     *                index_table(c_index_size,superindex)

c----------------------------------------------------------------------------
c   Loop over the segments of the superindex.
c----------------------------------------------------------------------------

         bseg = index_table(c_bseg, superindex)
         eseg = index_table(c_eseg, superindex)
         subsegment_counter = 0
         do j = bseg, eseg

c----------------------------------------------------------------------------
c   Recalculate the brange, erange of the superindex segment.
c----------------------------------------------------------------------------

            call segment_get_range(type, j, aosegs, naosegs,
     *                       norbitals, nocc_orbitals,
     *                       naocc_orbitals, nbocc_orbitals,
     *                       mx_ao_seg_size, mx_occ_seg_size,
     *                       mx_virt_seg_size,
     *                       brange, erange)

c----------------------------------------------------------------------------
c   Determine the correct segment size to use for this subindex segment.
c----------------------------------------------------------------------------

            segsize = get_subindex_segsize(type, j, sip_sub_segsize,
     *                      sip_sub_ao_segsize, sip_sub_occ_segsize,
     *                      sip_sub_virt_segsize, boccval, eoccval,
     *                      baoccval, eaoccval, bboccval, eboccval)

c----------------------------------------------------------------------------
c   Determine the number of subindex segments needed to fill up the 
c   superindex segment.
c----------------------------------------------------------------------------

            range   = erange-brange+1
            nsubseg = range/segsize
            if (segsize * nsubseg .ne. range) 
     *           nsubseg = nsubseg + 1  
            norb_subseg = range / nsubseg    ! orbitals per subsegment
            nrem        = range - norb_subseg * nsubseg   ! num extra segs

c----------------------------------------------------------------------------
c   Fill in the subindex segment table entries.
c----------------------------------------------------------------------------

            brange_s = brange
            do k = 1, nsubseg
               erange_s = brange_s + norb_subseg - 1
               if (k .le. nrem) erange_s = erange_s + 1
               subsegment_counter = subsegment_counter + 1
 
               ktable = itable + nsegment_table * lsegment_table_entry
               segtable(ktable + c_index - 1)   = i
               segtable(ktable + c_segment - 1) = subsegment_counter
               segtable(ktable + c_range1 - 1)  = brange_s
               segtable(ktable + c_range2 - 1)  = erange_s
               segtable(ktable + c_subseg1 - 1) = 0
               segtable(ktable + c_subseg2 - 1) = 0

               nsegment_table = nsegment_table + 1
               brange_s = erange_s + 1 
            enddo    ! k
         enddo       ! j
  300 continue       ! i

c-----------------------------------------------------------------------------
c   Now make another pass, indexing the subsegment ranges for each superindex
c   segment.  This will be used later in "do ii in i" constructs.      
c-----------------------------------------------------------------------------

      call set_subsegment_ranges(index_table, nindex_table, 
     *                           segtable(itable), nsegment_table)

c------------------------------------------------------------------------------
c   Determine the stack size for each individual type of data used.
c------------------------------------------------------------------------------

      call determine_stack_blocksizes(index_table, nindex_table,
     *         segtable(itable), nsegment_table, boccval, eoccval,
     *         baoccval, eaoccval, bboccval, eboccval,
     *         stack_blocksizes, nstacks)
      return
      end

      integer function translate_symbolic_constant(iseg, naosegs,
     *                              norbitals, nocc_orbitals,
     *                              naocc_orbitals, nbocc_orbitals, 
     *                              mx_ao_seg_size, mx_occ_seg_size,
     *                              mx_virt_seg_size)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      integer iseg, naosegs, norbitals, nocc_orbitals
      integer mx_ao_seg_size, mx_occ_seg_size, mx_virt_seg_size
      integer naocc_orbitals, nbocc_orbitals
      integer seg, bo, eo, bv, ev
      integer bao, eao, bav, eav
      integer bbo, ebo, bbv, ebv
      
      bo = 1
      eo = (nocc_orbitals - 1) / mx_occ_seg_size + 1
      bv = eo + 1
      ev = (norbitals - nocc_orbitals - 1) / mx_virt_seg_size + bv

      bao = 1
      eao = (naocc_orbitals - 1) / mx_occ_seg_size + 1
      bav = eao + 1
      eav = (norbitals -naocc_orbitals - 1) / mx_virt_seg_size + bav

      bbo = 1
      ebo = (nbocc_orbitals - 1) / mx_occ_seg_size + 1
      bbv = ebo + 1
      ebv = (norbitals - nbocc_orbitals - 1) / mx_virt_seg_size + bbv

      if (iseg .eq. norb) then
         seg = naosegs
      else if (iseg .eq. bocc) then
         seg = bo
      else if (iseg .eq. eocc) then
         seg = eo
      else if (iseg .eq. bvirt) then
         seg = bv
      else if (iseg .eq. evirt) then
         seg = ev
      else if (iseg .eq. naocc) then
         seg = naocc_orbitals
      else if (iseg .eq. baocc) then
         seg = bao
      else if (iseg .eq. eaocc) then
         seg = eao
      else if (iseg .eq. bavirt) then
         seg = bav
      else if (iseg .eq. eavirt) then
         seg = eav
      else if (iseg .eq. nbocc) then
         seg = nbocc_orbitals
      else if (iseg .eq. bbocc) then
         seg = bbo
      else if (iseg .eq. ebocc) then
         seg = ebo
      else if (iseg .eq. bbvirt) then
         seg = bbv
      else if (iseg .eq. ebvirt) then
         seg = ebv
      else if (iseg .eq. noccorb) then
         seg = nocc_orbitals
      else if (iseg .eq. nvirtorb) then
         seg = norbitals - nocc_orbitals
      else if (iseg .eq. boccorb) then
         seg = 1
      else if (iseg .eq. eoccorb) then
         seg = nocc_orbitals
      else if (iseg .eq. bvirtorb) then
         seg = nocc_orbitals + 1
      else if (iseg .eq. evirtorb) then
         seg = norbitals
      else if (iseg .eq. naoccorb) then
         seg =  naocc_orbitals
      else if (iseg .eq. nboccorb) then
         seg = nbocc_orbitals
      else if (iseg .eq. navirtorb) then
          seg = norbitals - naocc_orbitals
      else if (iseg .eq. nbvirtorb) then
          seg = norbitals - nbocc_orbitals 
      else if (iseg .eq. baoccorb) then
          seg = 1
      else if (iseg .eq. bboccorb) then
          seg = 1
      else if (iseg .eq. eaoccorb) then
          seg = naocc_orbitals
      else if (iseg .eq. eboccorb) then
          seg = nbocc_orbitals
      else if (iseg .eq. bavirtorb) then
          seg = naocc_orbitals + 1
      else if (iseg .eq. bbvirtorb) then
          seg = nbocc_orbitals + 1
      else if (iseg .eq. eavirtorb) then
          seg = norbitals
      else if (iseg .eq. ebvirtorb) then
          seg = norbitals 
      else if (iseg .eq. cc_iter_cons) then
         seg = cc_iter
      else if (iseg .eq. cc_hist_cons) then
         seg = cc_hist
      else if (iseg .eq. cc_beg_cons) then
         seg = cc_beg
      else if (iseg .eq. scf_iter_cons) then
         seg = scf_iter
      else if (iseg .eq. scf_hist_cons) then
         seg = scf_hist
      else if (iseg .eq. scf_beg_cons) then
         seg = scf_beg
      else if (iseg .eq. natoms_cons) then
         seg = ncenters
      else if (iseg .eq. subb_parm) then
         seg = subb
      else if (iseg .eq. sube_parm) then
         seg = sube
      else if (iseg .eq. sip_sub_segsize_parm) then
         seg = sip_sub_segsize
      else if (iseg .eq. sip_sub_occ_segsize_parm) then
         seg = sip_sub_occ_segsize
      else if (iseg .eq. sip_sub_virt_segsize_parm) then 
         seg = sip_sub_virt_segsize
      else if (iseg .eq. sip_sub_ao_segsize_parm) then 
          seg = sip_sub_ao_segsize 
      else  
         print *,'Error in translate_symbolic_constant: ',
     *         'Unknown symbolic constant ',iseg
         call abort_job()
      endif

      translate_symbolic_constant = seg
      return
      end                       

      subroutine segment_get_range(type, seg, aosegs, naosegs,
     *                       norbitals, nocc_orbitals, 
     *                       naocc_orbitals, nbocc_orbitals, 
     *                       mx_ao_seg_size, mx_occ_seg_size,
     *                       mx_virt_seg_size,
     *                       brange, erange)
      implicit none
      include 'interpreter.h'
      integer type, seg, naosegs, norbitals, nocc_orbitals
      integer naocc_orbitals, nbocc_orbitals
      integer no_per_seg, no_left, nv_per_seg, nv_left, nvorbitals
      integer nosegs, nvsegs
      integer mx_ao_seg_size, mx_occ_seg_size, mx_virt_seg_size,
     *        brange, erange
      integer aosegs(naosegs)

      if (type .eq. aoindex) then
         if (seg .eq. 1) then
            brange = 1
         else
            brange = aosegs(seg-1) + 1
         endif

         erange = aosegs(seg)
      else if (type .eq. moindex) then
         nosegs = (nocc_orbitals - 1) / mx_occ_seg_size + 1
         if (seg .le. nosegs) then
            no_per_seg = (nocc_orbitals - 1)/nosegs
            no_left    = nocc_orbitals - nosegs * no_per_seg
            if (seg .le. no_left) then
               brange = (seg-1) * (no_per_seg+1) + 1
               erange = min0(nocc_orbitals, brange + no_per_seg)
            else
               brange = no_left + (seg-1)*no_per_seg + 1
               erange = min0(nocc_orbitals, brange + no_per_seg-1)
            endif
         else
            nvorbitals = norbitals - nocc_orbitals
            nvsegs = (norbitals - nocc_orbitals - 1) / 
     *                                  mx_virt_seg_size + 1
            nv_per_seg = (nvorbitals - 1)/nvsegs
            nv_left    = nvorbitals - nvsegs * nv_per_seg

            if (seg .le. nosegs + nv_left) then
               brange = (seg-nosegs-1) * (nv_per_seg + 1) + 
     *                   nocc_orbitals + 1
               erange = min0(norbitals, brange + nv_per_seg)
            else
               brange = nv_left + (seg-nosegs-1)*nv_per_seg +
     *                  nocc_orbitals + 1
               erange =  min0(norbitals, brange + nv_per_seg - 1)
            endif
         endif

      else if (type .eq. moaindex) then
         nosegs = (naocc_orbitals - 1) / mx_occ_seg_size + 1
         nvsegs = (norbitals - naocc_orbitals - 1) / 
     *                               mx_virt_seg_size + 1
         if (seg .le. nosegs) then
            no_per_seg = (naocc_orbitals - 1)/nosegs
            no_left    = naocc_orbitals - nosegs * no_per_seg
            if (seg .le. no_left) then
               brange = (seg-1) * (no_per_seg+1) + 1
               erange = min0(naocc_orbitals, brange + no_per_seg)
            else
               brange = no_left + (seg-1)*no_per_seg + 1
               erange = min0(naocc_orbitals, brange + no_per_seg-1)
            endif
         else
            nvorbitals = norbitals - naocc_orbitals
            nv_per_seg = (nvorbitals - 1)/nvsegs
            nv_left    = nvorbitals - nvsegs * nv_per_seg

            if (seg .le. nosegs + nv_left) then
               brange = (seg-nosegs-1) * (nv_per_seg + 1) + 
     *                   naocc_orbitals + 1
               erange = min0(norbitals, brange + nv_per_seg)
            else
               brange = nv_left + (seg-nosegs-1)*nv_per_seg +
     *                  naocc_orbitals + 1
               erange =  min0(norbitals, brange + nv_per_seg - 1)
            endif
         endif

      else if (type .eq. mobindex) then
         nosegs = (nbocc_orbitals - 1) / mx_occ_seg_size + 1
         nvsegs = (norbitals - nbocc_orbitals - 1) / 
     *                                  mx_virt_seg_size + 1
         if (seg .le. nosegs) then
            no_per_seg = (nbocc_orbitals - 1)/nosegs
            no_left    = nbocc_orbitals - nosegs * no_per_seg
            if (seg .le. no_left) then
               brange = (seg-1) * (no_per_seg+1) + 1
               erange = min0(nbocc_orbitals, brange + no_per_seg)
            else
               brange = no_left + (seg-1)*no_per_seg + 1
               erange = min0(nbocc_orbitals, brange + no_per_seg-1)
            endif
         else
            nvorbitals = norbitals - nbocc_orbitals
            nv_per_seg = (nvorbitals - 1)/nvsegs
            nv_left    = nvorbitals - nvsegs * nv_per_seg

            if (seg .le. nosegs + nv_left) then
               brange = (seg-nosegs-1) * (nv_per_seg + 1) + 
     *                   nbocc_orbitals + 1
               erange = min0(norbitals, brange + nv_per_seg)
            else
               brange = nv_left + (seg-nosegs-1)*nv_per_seg +
     *                  nbocc_orbitals + 1
               erange =  min0(norbitals, brange + nv_per_seg - 1)
            endif
         endif
      else if (type .eq. simple_index) then
         brange = seg
         erange = seg
      else if (type .eq. laindex) then
         brange = 1
         erange = 1
      else if (type .eq. subindex) then
         brange = seg
         erange = seg
      else 
         print *,'Error in segment_get_range: Unknown type ',type
         call abort_job()
      endif
     
      return
      end

      subroutine check_index_range(seg, type, index_id, naosegs, 
     *                       ibocc, ieocc, ibvirt, ievirt, 
     *                       ibaocc, ieaocc, ibavirt, ieavirt,
     *                       ibbocc, iebocc, ibbvirt, iebvirt)
c--------------------------------------------------------------------------
c   Verifies that a given segment is within the actual range for its
c   index type. If an error is detected, an error message is printed,
c   and the job aborts.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'

      integer seg, type, naosegs, ibocc, ieocc, ibvirt, ievirt
      integer ibavirt, ieavirt, ibaocc, ieaocc
      integer ibbvirt, iebvirt, ibbocc, iebocc
      integer seg1, seg2, ierr
      integer index_id

      ierr = 0
      if (type .eq. aoindex) then
         if (seg .lt. 1 .or. seg .gt. naosegs) then
            ierr = 1
            seg1 = 1
            seg2 = naosegs
         endif
      else if (type .eq. moindex) then
         if (seg .lt. ibocc .or. seg .gt. ievirt) then
            ierr = 1
            seg1 = ibocc
            seg2 = ievirt
         endif
      else if (type .eq. moaindex) then
         if (seg .lt. ibaocc .or. seg .gt. ieavirt) then
            ierr = 1
            seg1 = ibaocc
            seg2 = ieavirt
         endif
      else if (type .eq. mobindex) then
         if (seg .lt. ibbocc .or. seg .gt. iebvirt) then
            ierr = 1
            seg1 = ibbocc
            seg2 = iebvirt
         endif
      else if (type .eq. simple_index) then
         if (seg .le. 0) then
            print *,'Error: Invalid range for index number ',index_id
            print *,'Range must be > 0, code has index value of ',
     *           seg
            call abort_job() 
         endif
      else if (type .eq. laindex) then
      else if (type .eq. subindex) then
      else
         print *,'Error: Invalid index type in check_index_range: ',
     *     'Index = ',index_id, ' Type = ',type
         call abort_job()
      endif

      if (ierr .eq. 1) then
         print *,'Error: Index range error for index ',index_id,
     *               ' type = ',type
         print *,'Index should be between ',seg1,' and ',seg2
         call abort_job()
      endif
      return
      end

      integer function get_subindex_segsize(superindex_type, seg, 
     *                      sip_sub_segsize,
     *                      sip_sub_ao_segsize, sip_sub_occ_segsize,
     *                      sip_sub_virt_segsize, boccval, eoccval,
     *                      baoccval, eaoccval, bboccval, eboccval)
c---------------------------------------------------------------------------
c   Function to return the correct segment size for a subindex segment 
c   corresponding to segment "seg" of a superindex of type "superindex_type".
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'

      integer superindex_type, seg
      integer sip_sub_segsize, sip_sub_ao_segsize, sip_sub_occ_segsize,
     *        sip_sub_virt_segsize
      integer boccval, eoccval, baoccval, eaoccval, bboccval, eboccval

      if (superindex_type .eq. aoindex) then
         get_subindex_segsize = sip_sub_ao_segsize
      else if (superindex_type .eq. simple_index) then
         get_subindex_segsize = sip_sub_segsize
      else if (superindex_type .eq. moindex) then
         if (seg .ge. boccval .and. seg .le. eoccval) then
            get_subindex_segsize = sip_sub_occ_segsize
         else
            get_subindex_segsize = sip_sub_virt_segsize
         endif
      else if (superindex_type .eq. moaindex) then
         if (seg .ge. baoccval .and. seg .le. eaoccval) then
            get_subindex_segsize = sip_sub_occ_segsize
         else
            get_subindex_segsize = sip_sub_virt_segsize
         endif
      else if (superindex_type .eq. mobindex) then
         if (seg .ge. bboccval .and. seg .le. eboccval) then
            get_subindex_segsize = sip_sub_occ_segsize
         else
            get_subindex_segsize = sip_sub_virt_segsize
         endif
      endif

      return
      end 

      subroutine set_subsegment_ranges(index_table, nindex_table,
     *                 segment_table, nsegment_table)
c-------------------------------------------------------------------------
c   For each superindex segment, fills in the beginning and ending 
c   segment of its corresponding subindex.
c-------------------------------------------------------------------------
     
      implicit none
      include 'interpreter.h'
      integer nsegment_table
      integer nindex_table 
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)

      integer i, j, k
      integer subindex_seg, superindex, super_bseg, super_eseg
      integer bseg, eseg

      do k = 1, nindex_table
         if (index_table(c_index_type, k) .eq. subindex) then

c----------------------------------------------------------------------------
c   We have a subindex.  Process all the segments of its superindex.
c----------------------------------------------------------------------------

            superindex = index_table(c_subindex_ptr,k)
            super_bseg = index_table(c_bseg, superindex)
            super_eseg = index_table(c_eseg, superindex)
    
            do j = super_bseg, super_eseg

c---------------------------------------------------------------------------
c   Find the bseg and eseg of the subindex corresponding to superindex
c   segment j.
c---------------------------------------------------------------------------

               call get_subrange_segments(superindex, j, k,
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           bseg, eseg)

c----------------------------------------------------------------------------
c   Now look up the segment_table entry for this superindex segment, and
c   store the bseg and eseg.  This data may be used in a "do ii in i"
c   type loop.
c----------------------------------------------------------------------------

               do i = 1, nsegment_table
                  if (segment_table(c_index,i) .eq. superindex .and.
     *                segment_table(c_segment,i) .eq. j) then
                     segment_table(c_subseg1,i) = bseg
                     segment_table(c_subseg2,i) = eseg
                     go to 100
                  endif
               enddo   ! i

  100          continue
            enddo  !  j
         endif
      enddo    ! k  

      return
      end


