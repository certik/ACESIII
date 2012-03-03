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
      subroutine dcont2(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Performs a contraction of a block of the gamma matrix with the 
c   2-electron derivative integrals on an atom-by-atom basis.  The result of
c   each contraction is summed into its appropriate component of the 
c   gradient array.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'machine_types.h'
      include 'saved_data.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer i, j, n
      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer array, ind(mx_array_index), nindex, ierr
      integer block, blkndx, seg
      integer*8 iderblk(4,3), indblk,  
     *           get_block_data_index
      integer stack, handle, handle2, icenter, match
      integer find_current_block
 
      integer val1(mx_array_index), val2(mx_array_index)
      integer iatom, igrad, flopcount
      integer component
      
      double precision x(1), y
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      integer dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,
     *        derint_counter
      integer derint_array(12)
      common /derint/dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,
     *               derint_counter
      equivalence (derint_array(1), dx1)

      if (derint_counter .ne. 12) then
         print *,'Error in DCONT2: derivative_integrals mapping ',
     *     'is not set properly. Number of arrays is ',derint_counter
         call abort_job()
      endif

#ifdef ALTIX
      dptr = dshptr
#endif

c----------------------------------------------------------------------------
c   Locate the data for the intermediate block.
c---------------------------------------------------------------------------

      array = op(c_result_array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
         seg = index_table(c_current_seg,ind(i)) 
         call get_index_segment(ind(i), seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
      enddo

      a1 = val1(1)
      a2 = val2(1)
      b1 = val1(2)
      b2 = val2(2)
      c1 = val1(3)
      c2 = val2(3)
      d1 = val1(4)
      d2 = val2(4)

      block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

      stack = array_table(c_array_stack,array)
      indblk = get_block_data_index(array, block, stack,
     *                                        blkndx, x)
 
c--------------------------------------------------------------------------
c   Get the address of each of the derivative integral blocks.
c---------------------------------------------------------------------------

      n = 0
      do i = 1, 4
      do j = 1, 3
         n = n + 1 
         array = derint_array(n)

         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack = array_table(c_array_stack,array)
         iderblk(i,j) = get_block_data_index(array, block, stack,
     *                                        blkndx, x)
      enddo
      enddo

      do iatom = 1, ncenters
         do component = 1, 3
         do icenter   = 1, 4
  
c-----------------------------------------------------------------------------
c   Process each atom's contribution to the gradient.
c-----------------------------------------------------------------------------

            call sum_atom_data(iatom, x(indblk), 
     *                      x(iderblk(icenter,component)),
     *                      a1,a2,b1,b2,c1,c2,d1,d2, 
     *                      icenter, match, y)

c-----------------------------------------------------------------------------
c   Sum the contribution into its proper component for this atom.
c-----------------------------------------------------------------------------

            igrad = (iatom-1)*3 + component
            gradient_data(igrad) = gradient_data(igrad) + y
         enddo
         enddo 
      enddo

      return
      end

      subroutine sum_atom_data(iatom, data, derv, 
     *                  a1,a2,b1,b2,c1,c2,d1,d2, 
     *                  icenter, match, sum)
c---------------------------------------------------------------------------
c   Extracts the derivative integrals  and density data for a particular atom 
c   into blocks. 
c---------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      include 'parallel_info.h'

      integer iatom, icenter, match
      integer a1,a2,b1,b2,c1,c2,d1,d2
      double precision data(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision derv(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision sum

      integer a, b, c, d
      logical atom_match
      logical aflag(a1:a2)
      logical bflag(b1:b2)
      logical cflag(c1:c2) 
      logical dflag(d1:d2)
     
      sum = 0.
      match = 0
      if (icenter .eq. 1) then
         do a = a1, a2
            aflag(a) = atom_match(a,iatom)
            if (aflag(a)) match = match + 1
         enddo

         if (match .eq. 0) return

         do d = d1, d2
         do c = c1, c2
         do b = b1, b2
         do a = a1, a2
            if (aflag(a)) sum = sum + derv(a,b,c,d) * data(a,b,c,d)
         enddo
         enddo
         enddo
         enddo
      else if (icenter .eq. 2) then
         do b = b1, b2
            bflag(b) = atom_match(b,iatom)
            if (bflag(b)) match = match + 1
         enddo

         if (match .eq. 0) return

         do b = b1, b2
            if (bflag(b)) then
               do d = d1, d2
               do c = c1, c2
               do a = a1, a2
                  sum = sum + derv(a,b,c,d) * data(a,b,c,d)
               enddo
               enddo
               enddo
            endif
         enddo
      else if (icenter .eq. 3) then
         do c = c1, c2
            cflag(c) = atom_match(c,iatom)
            if (cflag(c)) match = match + 1
         enddo

         if (match .eq. 0) return

         do c = c1, c2
            if (cflag(c)) then
               do d = d1, d2
               do b = b1, b2
               do a = a1, a2
                  sum = sum + derv(a,b,c,d) * data(a,b,c,d)
               enddo
               enddo
               enddo
            endif
         enddo
      else if (icenter .eq. 4) then
         do d = d1, d2
            dflag(d) = atom_match(d,iatom)
            if (dflag(d)) match = match + 1
         enddo

         if (match .eq. 0) return

         do d = d1, d2
            if (dflag(d)) then
               do c = c1, c2
               do b = b1, b2
               do a = a1, a2
                  sum = sum + derv(a,b,c,d) * data(a,b,c,d)
               enddo
               enddo
               enddo
            endif
         enddo
      endif

      return
      end

      logical function atom_match(ao_orbital, iatom)
c--------------------------------------------------------------------------
c   Returns .true. if the ao_orbital comes from atom "iatom", .false.
c   otherwise.
c--------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      integer ao_orbital, iatom
      integer m

c---------------------------------------------------------------------------
c   Convert the AO orbital into a shell block.
c---------------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, ao_orbital, m)

c---------------------------------------------------------------------------
c   Check the atom corresponding to the shell block index for a match.
c---------------------------------------------------------------------------

      if (atom(m) .eq. iatom) then
         atom_match = .true.
      else
         atom_match = .false.
      endif
      return
      end

