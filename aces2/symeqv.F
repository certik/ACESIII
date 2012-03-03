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
      Subroutine SymEqv (NAtms,GenBy)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Prepares DB's GENBY list in a much less cumbersome way, using
C  Findif data structures.  Works for all point groups.
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ORIENT
#include "mxatms.par"
      Common /Orient/ Orient(3,3)
      DIMENSION MEMBER(MXATMS),ORBPOP(MXATMS),GENBY(NATMS)
C
C INITIALIZE GENBY TO 999.  THIS WILL TAKE CARE OF DUMMY ATOMS.
C
      DO 20 I=1,NATMS
       GENBY(I)=999
20    CONTINUE
C
C GET MEMBER AND ORBPOP VECTORS FROM JOBARC.
C
      ISIZE=1
      CALL IGETREC(20,'JOBARC','COMPNORB',ISIZE,IORBIT)
      CALL IGETREC(20,'JOBARC','COMPMEMB',NATMS,MEMBER)
      CALL IGETREC(20,'JOBARC','COMPPOPV',IORBIT,ORBPOP)
C
C LOOP OVER ORBITS AND ZERO OUT POSITION OF FIRST MEMBER OF EACH ORBIT
C   IN GENBY LIST.
C
      IOFF=1
      DO 50 I=1,IORBIT
       GENBY(MEMBER(IOFF))=0
       DO 51 J=1,ORBPOP(I)-1
        GENBY(MEMBER(IOFF+J))=MEMBER(IOFF)
51     CONTINUE
       IOFF=IOFF+ORBPOP(I)
50    CONTINUE
      RETURN
      END
