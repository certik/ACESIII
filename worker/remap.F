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
      SUBROUTINE REMAPX(CMPPGP,IORDGP,NATOMS,NORBCMP,COORD,
     &                 SCR,IMEMCMP,IPOPCMP,IMAP,itype)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: remap.F,v 1.3 2010/06/30 15:58:42 ponton Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     remap -- determine the VMol -> ZMAT reordering vector
C
C SYNOPSIS
      CHARACTER*4 CMPPGP
      Integer IOrdGp, NAtoms, NOrbCmp
      Double precision COORD(3,NATOMS), SCR(3*(IOrdGp-1))
      Integer IMEMCMP(NATOMS),IPOPCMP(NORBCMP), IMAP(NATOMS)
C
C DESCRIPTION
C     THIS ROUTINE FORMS A POINTER VECTOR RELATING THE ATOMIC ORDERING
C     USED IMPLICITY BY VMOL TO THAT SPECIFIED IN THE Z-MATRIX.  THIS
C     IS RATHER EASILY CALCULATED FROM INFORMATION RESIDENT IN THE JOBARC
C     FILE.  UPON RETURN, THE INTEGER VECTOR IMAP CONTAINS THE 
C     Z-MATRIX POSITION OF THE ATOMS, AS INDEXED BY THE VMOL ORDERING.
C     FOR EXAMPLE, IF IMAP(3)=11, THIS MEANS THAT THE 3RD ATOM USED BY
C     VMOL IS ACTUALLY THE ELEVENTH ATOM IN THE Z-MATRIX.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Integer IOff, INumAtm, IAtom, IOrbit, IRefAtm, NGen, ICent,
     $   ILoc, IOff1, I, IPos
      Integer itype
      LOGICAL QSAME

      do i = 1, natoms
         imap(i) = 0
      enddo

C LOOP OVER ORBITS IN COMPUTATIONAL POINT GROUP
C
      IOFF=1
      IATOM=1
      DO 10 IORBIT=1,NORBCMP
C
C DETERMINE THE REFERENCE ATOM FOR THE ORBIT (THE ONE USED BY VMOL
C   TO GENERATE THE REDUNDANT CENTERS)
C 
       INUMATM=IPOPCMP(IORBIT)
       IREFATM=IMEMCMP(IOFF)
       IMAP(IATOM)=IREFATM
       IATOM=IATOM+1
C
C NOW USE THE ALGORITHM OF VMOL TO GENERATE THE REDUNDANT CENTERS
C   IN THE COMPUTATIONAL POINT GROUP.  THE VALUE OF NGEN IS THE
C   NUMBER OF DISTINCT REDUNDANT CENTERS RETURNED
C
       CALL VMLGENX(CMPPGP,IORDGP,COORD(1,IREFATM),SCR,NGEN, itype)
C
C NOW LOOP OVER DISTINCT REDUNDANT CENTERS AND FIND THEIR 
C  POSITION IN THE Z-MATRIX ORDER
C
       DO 11 ICENT=1,NGEN
        ILOC=-1
        IOFF1=IOFF+1
        DO 20 I=1,INUMATM-1
         IPOS=IMEMCMP(IOFF1)
         IF(QSAME(SCR(1+(ICENT-1)*3),COORD(1,IPOS)))ILOC=IPOS
         IOFF1=IOFF1+1
20      CONTINUE
        IMAP(IATOM)=ILOC
        IATOM=IATOM+1
11     CONTINUE
       IOFF=IOFF+INUMATM
10    CONTINUE
      RETURN
      END
