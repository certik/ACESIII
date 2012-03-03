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
      SUBROUTINE SORTXYZ(XX,Y,NORD,NATOMS)
C
C     SORTS VECTOR OF NUCLEAR COORDINATES - TO CHECK FOR EQUIVALENCE
C     OF TWO ORIENTATIONS - NEEDS Q VECTOR AND ATOMIC MASS VECTOR 
C     FOR MASS WEIGHTING.  THE SORT IS FIRST DONE ON THE X COORDINATE,
C     BUT IF IDENTICAL X VALUES ARE ENCOUNTERED, THE Y COORDINATE
C     IS THEN SORTED AS WELL, AND ALSO THE Z IF NECESSARY.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     Maximum number of atoms currently allowed
#include "mxatms.par"
#include "coord_4sortxyz.com"
C
      COMMON /TOLERS/ SYMTOL,DEGTOL
C
C The purpose of the above coord common block was to get the value of the
C number of atoms, mass, and charge. However, it was a serious mistake
C to obtain the number of atoms from this common block.
C NATOMS in the COORD common block corresponds to the total number
C of atoms including the dummy atoms. There were calls to this
C routine from IRREPS() and COMVEC() with the assumption
C that NATOMS corresponds to the number of real atoms. This
C caused memory errors for sufficiently large systems (>20 atoms), but
C was hidden for smaller systems. The correct thing
C to do is to let the caller pass in the correct value as an
C argument and ignore the last entry of the COORD common block.
C 01/2006, Ajith Perera.
C
      DIMENSION X(3*MXATMS),XX(3*NATOMS),Y(3*NATOMS),NORD(MXATMS*2)
      DIMENSION SCR(3*MXATMS)
C
      CALL ZERO(X,3*MXATMS)
C
C SORT THE X COORDINATES AND MASS WEIGHT THEM
C
      CALL XDCOPY(3*NATOMS,XX,1,SCR,1)
      IOFF=1
      DO 10 I=1,NATOMS
       NORD(I+MXATMS)=I
       IOFF=IOFF+3
10    CONTINUE
    
      CALL XDCOPY(NATOMS,SCR,3,X,1)
C
C GIVE DUMMY ATOMS RIDICULOUS X COORDINATES
C
      DO 5 I=1,NATOMS
       IF(ATMASS(I).LT.1.D-8)X(I)=-99999.
5     CONTINUE
      CALL PIKSR3(NATOMS,X,NORD(MXATMS+1))
      DO 11 ITARGET=1,NATOMS
       ISOURCE=NORD(ITARGET+MXATMS)
       IOFFSRC=3*(ISOURCE-1)+1
       IOFFTAR=3*(ITARGET-1)+1
       CALL XDCOPY(3,SCR(IOFFSRC),1,Y(IOFFTAR),1)
11    CONTINUE
C
C NOW WE MUST SEARCH FOR CLUSTERS OF VALUES IN THE X VECTOR
C
      IFINDX=1
1     NLEFT=NATOMS-IFINDX+1
      ILOC=IFINDNE(NLEFT,X(IFINDX),1,X(IFINDX),SYMTOL)+IFINDX-1
      ICLSIZX=ILOC-IFINDX
      IF(ILOC.NE.IFINDX+1)THEN
C
C CLUSTER HAS BEEN LOCATED.  SORT THE ASSOCIATED Y VALUES.
C
       CALL XDCOPY(ICLSIZX,Y(2+3*(IFINDX-1)),3,X,1)
       CALL PIKSR3(ICLSIZX,X,NORD(MXATMS+IFINDX))
       DO 13 ITARGET0=1,ICLSIZX
        ISOURCE=NORD(ITARGET0+MXATMS+IFINDX-1)
        ITARGET=ITARGET0+IFINDX-1
        IOFFSRC=3*(ISOURCE-1)+1
        IOFFTAR=3*(ITARGET-1)+1
        CALL XDCOPY(3,SCR(IOFFSRC),1,Y(IOFFTAR),1)
13     CONTINUE
C
C NOW WE MUST SEARCH FOR CLUSTERS OF VALUES IN THE Y VECTOR WITHIN
C  THIS X VECTOR CLUSTER
C
       IFINDY=1
2      NLEFTY=ICLSIZX-IFINDY+1
       ILOC=IFINDNE(NLEFTY,X(IFINDY),1,X(IFINDY),SYMTOL)+IFINDY-1
       ICLSIZY=ILOC-IFINDY
       IF(ICLSIZY.NE.1)THEN
C
C CLUSTER HAS BEEN LOCATED.  SORT THE Z VALUES.
C
        IPOS=IFINDX+IFINDY-1
        CALL XDCOPY(ICLSIZY,Y(3+3*(IPOS-1)),3,X,1)
        CALL PIKSR3(ICLSIZY,X,NORD(MXATMS+IPOS))
        DO 14 ITARGET0=1,ICLSIZY
         ISOURCE=NORD(ITARGET0+MXATMS+IPOS-1)
         ITARGET=ITARGET0+IPOS-1
         IOFFSRC=3*(ISOURCE-1)+1
         IOFFTAR=3*(ITARGET-1)+1
         CALL XDCOPY(3,SCR(IOFFSRC),1,Y(IOFFTAR),1)
14      CONTINUE
       ENDIF
       IFINDY=IFINDY+ICLSIZY
       IF(IFINDY.LT.ICLSIZX)GOTO 2
      ENDIF
      IFINDX=IFINDX+ICLSIZX
      IF(IFINDX.LT.NATOMS)GOTO 1
C
      IOFF=1
      DO 100 I=1,NATOMS
       Z=FLOAT(IATNUM(NORD(MXATMS+I)))
       CALL DSCAL(3,Z,Y(IOFF),1)
       IOFF=IOFF+3
100   CONTINUE
      CALL ICOPY(MXATMS,NORD(MXATMS+1),1,NORD,1)
      RETURN
      END
