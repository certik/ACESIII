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

      SUBROUTINE STRIPDUM(IPTR,ATMASS,NATOM,NTOTAL,NSYOP,IHOLD)
C
C STRIP DUMMY ATOMS OUT OF PERMUTATION ARRAY
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPTR(NTOTAL,NSYOP),ATMASS(*),IMAPZ2R(150)
      DIMENSION IMAPR2Z(150),IHOLD(NTOTAL,NSYOP)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      CALL DGETREC(20,'JOBARC','ATOMMASS',NTOTAL,ATMASS)
c      CALL IZERO(IMAPZ2R,NTOTAL*NSYOP)
      do i = 1, NTOTAL*NSYOP
         IMAPZ2R(i) = 0
         IMAPR2Z(i) = 0
      enddo

c      CALL IZERO(IMAPR2Z,NTOTAL*NSYOP)
      CALL ICOPY (NTOTAL*NSYOP, IPTR, 1, IHOLD, 1)

      IREAL=0
      DO 10 IATOM=1,NTOTAL
       IF(ATMASS(IATOM).NE.0.0D0)THEN
        IREAL=IREAL+1
        IMAPR2Z(IREAL)=IATOM
        IMAPZ2R(IATOM)=IREAL
       ENDIF
10    CONTINUE
C
      DO 20 IOP=1,NSYOP
       DO 30 IATOM=1,NATOM
        JATOM=IMAPR2Z(IATOM)
        KATOM=IHOLD(JATOM,IOP)
        IPTR(IATOM,IOP)=IMAPZ2R(KATOM)
30     CONTINUE
20    CONTINUE
C
      RETURN
      END        
