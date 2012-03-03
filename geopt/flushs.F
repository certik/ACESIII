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
C
C FLUSHES SYMMETRY INFORMATION TO FILE SYMINF.
C
      SUBROUTINE FLUSHS(IORDER,IORBIT,SYMOPS,IPTR,MEMBER,ORBPOP,
     &                  CLSTYP,NATOMS,LABEL,PTGRP,NCLASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER CLSTYP,ORBPOP
#include "mxatms.par"
      CHARACTER*4 PTGRP
      CHARACTER*4 LABEL
      CHARACTER*8 FULLAB,szJunk,STSYM(MXATMS)
      character*(8*mxatms) szStSymTmp
      DIMENSION SYMOPS(9*IORDER),IPTR(NATOMS,IORDER),MEMBER(NATOMS)
      DIMENSION ORBPOP(IORBIT),CLSTYP(IORBIT)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /LOCAL/ STSYM
C      COMMON /LOCAL/ STSYM
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)

      szJunk(1:4)=PTGRP
      szJunk(5:8)=' '
C
C PUT STUFF OUT ON JOBARC.
C
C
C WRITE OUT POINT GROUP, ORDER AND NUMBER OF CLASSES.
C

csb
      if(label.eq.'FULL' .and. iflags(46).gt.0) then
        write(6,1000)
        write(6,1002) szJunk
        write(6,1004) iorder
        write(6,1008)
        write(6,1009) (orbpop(i),i=1,natoms)
        write(6,1010)
        write(6,1009) (member(i),i=1,natoms)
        write(6,1012)
      endif
 1000 format('*** ORBIT')
 1002 format(' PTGRP: ',a4)
 1004 format(' ORDR : ',i4)
 1008 format(' POPV : ')
 1009 format(15(i5))
 1010 format(' MEMB : ')
 1012 format('*** END')

      ISIZE=1
      FULLAB=LABEL//'PTGP'
CSS      CALL IGETREC(0, 'JOBARC', FULLAB, length, Ijunk)
CSS      Print*, "The length", length, FULLAB, szjunk
      CALL PUTCREC(20,'JOBARC',FULLAB,8,szJunk)
      FULLAB=LABEL//'ORDR'
CSS      CALL IGETREC(0, 'JOBARC', FULLAB, length, Ijunk) 
      CALL IPUTREC(20,'JOBARC',FULLAB,ISIZE,IORDER)
      FULLAB=LABEL//'NIRR'
      CALL IPUTREC(20,'JOBARC',FULLAB,ISIZE,NCLASS)
C
C WRITE OUT THE NUMBER OF ORBITS
C
      FULLAB=LABEL//'NORB'
      CALL IPUTREC(20,'JOBARC',FULLAB,ISIZE,IORBIT)
C
C WRITE OUT SYMMETRY OPERATIONS
C
      FULLAB=LABEL//'SYOP'
      ISIZE=9*IORDER
      CALL DPUTREC(20,'JOBARC',FULLAB,ISIZE,SYMOPS)
C
C WRITE OUT PERMUTATION VECTORS
C
      FULLAB=LABEL//'PERM'
      ISIZE=IORDER*NATOMS
      CALL IPUTREC(20,'JOBARC',FULLAB,ISIZE,IPTR)
C
C WRITE OUT ORBIT CENSUS VECTOR AND ORBIT POPULATION VECTOR.
C
      FULLAB=LABEL//'MEMB'
      ISIZE=NATOMS
      CALL IPUTREC(20,'JOBARC',FULLAB,ISIZE,MEMBER)
      FULLAB=LABEL//'POPV'
c      ISIZE=IORBIT
      CALL IPUTREC(20,'JOBARC',FULLAB,ISIZE,ORBPOP)
C
C WRITE OUT CLASS TYPE VECTOR
C
      FULLAB=LABEL//'CLSS'
      ISIZE=IORDER
      CALL IPUTREC(20,'JOBARC',FULLAB,ISIZE,CLSTYP)
C
C WRITE OUT SITE GROUP VECTOR
C
c   o 'expand' STSYM into a 1-dimensional char array
      FULLAB=LABEL//'STGP'
c      ISIZE=IORBIT
      ISIZE=NATOMS
      iNdx = 1
      do i = 1, iSize
         szStSymTmp(iNdx:iNdx+7) = STSYM(i)(1:8)
         iNdx = iNdx + 8
      end do
      iSize = iSize*8
      CALL PUTCREC(20,'JOBARC',FULLAB,iSize,szStSymTmp(1:iSize))

      CLOSE(UNIT=30,STATUS='KEEP')
      RETURN
      END
