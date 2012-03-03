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
      subroutine dputrec(iFlag,szArchive,szRecName,iRecLen,dSrc)
c-----------------------------------------------------------------------------
c   Writes double precision data to JOBARC file.
c
c   The iRecLen argument is the number of double precision words contained 
c   in dSrc.
c-----------------------------------------------------------------------------

      implicit none

      integer iFlag, iRecLen
      character*(*) szArchive, szRecName

      double precision dSrc(*) 
      integer i, iRecLen_int
      integer iconvert

#include "machsp.com"

      iconvert    = ifltln / 8
      iRecLen_int = iRecLen * iconvert
      call putrec(iFlag,szArchive,szRecName,iRecLen_int,dSrc)
      return
      end
