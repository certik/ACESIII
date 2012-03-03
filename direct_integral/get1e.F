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
      SUBROUTINE GET1E(LUINT, type, x, LDIM,ILNBUF)
C
      IMPLICIT none
C
      integer luint, type, ilnbuf
      integer int, i, j, ldim, nxt
      double precision x(ldim,ldim)
      double precision BUF(ILNBUF)
      double precision scr(ldim*(ldim+1)/2)
#ifdef ACES64
      integer*8 IBUF(ILNBUF)
      integer*8 nut
#else
      integer IBUF(ILNBUF)
      integer nut
#endif

      character*8 ctype

      do i = 1, ldim
         scr(i) = 0.
      enddo

      if (type .eq. 1) then
         return   ! summed in with kinetic
      else if (type .eq. 2) then
         ctype = 'ONEHAMIL'   ! sum of kinetic and attraction
      else
         ctype = 'OVERLAP '
      endif

      CALL LOCATE(LUINT, ctype)
      NUT = ILNBUF
      DO WHILE (NUT.EQ.ILNBUF)
         READ(LUINT) BUF, IBUF, NUT
         DO INT = 1, NUT
            scr(IBUF(INT)) = BUF(INT)
         END DO
      END DO

c------------------------------------------------------------------------
c   Expand the data to its full (symmetric) form.
c------------------------------------------------------------------------

      nxt = 0
      do j = 1, ldim
      do i = 1, j
         nxt = nxt + 1
         x(i,j) = scr(nxt)
         x(j,i) = scr(nxt)
      enddo
      enddo

      RETURN
      END


      subroutine locate(iunit,labels)

c Finds a label in a vmol integral file and positions pointer to the
c first record after the label.

      integer iunit
      character*(*) labels

      character*8 lab,label

      label=labels
      rewind (iunit)
   10 read (iunit,end=900,err=910) lab
      if (lab(1:2).ne.'**') goto 10
      read (iunit,end=900) lab
      if (lab.eq.label) return
      goto 10

  900 write (6,9000) label,iunit
 9000 format(t3,'@LOCATE-F, Label ',a8,' not found on unit ',i3,'.')
      call abort_job()
      return

  910 write (6,9010) iunit,label
 9010 format(t3,'@LOCATE-F, I/O error on unit ',i3,' looking for ',
     &    'label ',a8,'.')
      call abort_job()
      return
      end

