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
      subroutine identfy_fragments(nratms,ibndto,connected_atoms,
     &                             nfrags, ifraglen,ifragmat)
C     This routine separates all of the atoms into disconnected fragments.
C     Each column of ifragmat will have the atom numbers of all the atoms
C     in that fragment.  Fragment lengths are stored in ifraglen, and 
C     the number of fragments is nfrags.
C     This routine calls the logical function is_connected, which determines
C     whether two atoms have a path connecting them.
C     AGT 11/06
      implicit none
C     Input variables
      integer nratms
      integer ibndto(nratms,nratms),connected_atoms(nratms)
C     Output variables
      integer nfrags,ifraglen(nratms),ifragmat(nratms,nratms)
C     Function variables
      logical is_connected
C     External Procedures
      external is_connected
C     Local Variables
      integer iatoms,jatoms,unfrag,unfrag2,ii
      logical iam_connected
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     Initial number of fragments = 0
      nfrags=0
      call izero(ifragmat,nratms*nratms)
      Call izero(ifraglen,nratms)
C     Initialize all atoms in first column of ifragmat
      do iatoms=1,nratms
         ifragmat(1,iatoms)=iatoms
      end do
C     Unfrag is a counter for the atoms where we haven't determined what
C     fragment they are on.
      unfrag=nratms
C     Maximum number of fragments is the number of atoms
C
      do while ((unfrag.ne.0).and.(nfrags.lt.nratms))
C     If all atoms have been assigned to fragments, stop.
         nfrags=nfrags+1
         ifraglen(nfrags)=1
         iatoms=ifragmat(nfrags,1)
C     unfrag2: Temporary counter for the changing number of unassigned
C              atoms
         unfrag2=unfrag-1
         do ii=2,unfrag
            jatoms=ifragmat(nfrags,ii)
            iam_connected = is_connected(iatoms,jatoms,ibndto,nratms,
     &                       connected_atoms)
            if (iam_connected) then
C     If connected, put it in the correct place and in/de - crement counters
               ifraglen(nfrags)=ifraglen(nfrags)+1
               ifragmat(nfrags,ifraglen(nfrags))=jatoms
               unfrag2=unfrag2-1
            else
C     Otherwise, move atom to the next column, and we will repeat
               ifragmat(nfrags+1,ii-(unfrag-unfrag2))=jatoms               
            endif
         end do
         do ii=ifraglen(nfrags)+1,unfrag
C     To make sure there is no confusion, zero out extraneous elements
            ifragmat(nfrags,ii)=0
         end do
         unfrag=unfrag2
      end do
c
      RETURN
      END
