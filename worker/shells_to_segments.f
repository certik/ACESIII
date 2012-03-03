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
      subroutine shells_to_segments(end_nfps, nshells, mx_seg_size,
     *                        segs, nsegs, atom_based, atoms)
c---------------------------------------------------------------------------
c   Fit the shells into index segments, based on the maximum segment size.
c
c   If the "atom_based" flag is .true., we will use a segmentation scheme 
c   that forces the segments to occur on atomic boundaries.  Otherwise,
c   the shells are arranged for most efficient packing within the segment
c   boundaries.
c---------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'

      integer nshells, mx_seg_size, nsegs
      integer end_nfps(nshells), segs(nshells)
      integer atoms(nshells), last_atom
      logical atom_based

      integer i, ishell, nsh, sum 

      ishell = 1
      nsegs  = 0
      nsh    = end_nfps(1)
      sum = 0
      last_atom = atoms(1)
  100 continue
      if (nsh .gt. mx_seg_size) then
         print *,'Error: Shell ',ishell,' > max. seg. size of ',
     *           mx_seg_size
         print *,'   Shell contains ',nsh,' functions.'
         call abort_job()
	endif

      if (atom_based) then
         if (atoms(ishell) .ne. last_atom) then

c-----------------------------------------------------------------------------
c   This is the boundary of a new atom.  Force a segment at this point.
c-----------------------------------------------------------------------------

            last_atom   = atoms(ishell)   ! thsi is the atom for the next seg.
            ishell      = ishell - 1
            nsegs       = nsegs + 1
            segs(nsegs) = end_nfps(ishell) 
            sum         = 0
            go to 150
         endif
      endif

      sum = sum + nsh
      if (sum .gt. mx_seg_size) then
         last_atom   = atoms(ishell)
         ishell      = ishell - 1
         nsegs       = nsegs + 1         ! increment the number of segments
         segs(nsegs) = end_nfps(ishell)  ! set the segment size
         sum         = 0                 ! reset counter
      else if (sum .eq. mx_seg_size) then
         nsegs       = nsegs + 1         ! increment the number of segments
         segs(nsegs) = end_nfps(ishell)  ! set the segment size
         sum         = 0                 ! reset counter
      endif

  150 continue
      ishell    = ishell + 1
      if (ishell .gt. nshells) then
         if (sum .gt. 0) then
            nsegs = nsegs + 1
            segs(nsegs) = end_nfps(nshells)
         endif
         go to 200
      endif
      nsh = end_nfps(ishell) - end_nfps(ishell-1)
      go to 100

  200 continue

      return
      end
