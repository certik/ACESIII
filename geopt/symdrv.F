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

c This is a switch for debugging. Apparently, the CRAYs do not like
c storing characters in integer scratch arrays. Well, porting this
c to Sun/SPARC causes problems. If anyone decides to spend more time
c on this than me, then switching this define on should activate all
c the right code.
#define _USE_OPRTYP

      SUBROUTINE SYMDRV(ICORE,ZCORE,MAXINT,MAXREA,DUMMY,LABEL)
      IMPLICIT INTEGER (A-Z)

#include "mxatms.par"
#include "coord.com"
C
      DOUBLE PRECISION ZCORE
      CHARACTER*4 LABEL
      CHARACTER*4 PTGRP,DUMMY
      CHARACTER*1 PMOD

#ifdef _USE_OPRTYP
      CHARACTER*1 OPRTYP(100)
#endif /* _USE_OPRTYP */

      DIMENSION ICORE(MAXINT), ZCORE(MAXREA)
      COMMON /FLAGS/ IFLAGS(100), IFLAGS2(500)
      CHARACTER*8 STSYM(MXATMS)
      COMMON /LOCAL/ STSYM
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4     FPGrp, BPGrp, PGrp
      Common /PtGp_com/ FPGrp, BPGrp, PGrp

      PTGRP=DUMMY
      IF(DUMMY(1:2).EQ.'DX')PTGRP='D8h '
      IF(DUMMY(1:2).EQ.'CX')PTGRP='C8v '
      IORDER=IORGRP(PTGRP)

C GET ORDER OF HIGHEST AXIS
      IL=linblnk(PTGRP)
      PMOD=PTGRP(IL:IL)
      IF ((PMOD.EQ.'d').OR.(PMOD.EQ.'h').OR.(PMOD.EQ.'v')) THEN
         NORDER=ATOI(PTGRP(2:IL-1))
      ELSE
         NORDER=ATOI(PTGRP(2:IL))
      END IF

C CONVERT TO CANONICAL ORIENTATION IF THIS IS THE CALL FOR THE
C FULL POINT GROUP
      IF ((LABEL.EQ.'FULL').OR.(FPGRP.EQ.PTGRP))
     &   CALL STDORT(Q,ICORE,3*NATOMS,0)

C ALLOCATE CORE FOR CALL TO CHRTAB
      IC1=1
      IC2=IC1+120*9
      IC3=IC2+9*NATOMS
      I1=1
      I2=I1+IORDER
#ifdef _USE_OPRTYP
c YAU - OPRTYP() replaces ICORE(I1)
      I2=1
#endif /* _USE_OPRTYP */
      I3=I2+IORDER
      I4=I3+IORDER
      I5=I4+IORDER
      I6=I5+IORDER
      I8=I6+IORDER
      I9=I8+NATOMS*IORDER
      I10=I9+NATOMS
C 4/8/97
C I10 needs 2*MXATMS. See what happens to array NORD in TSTOPS, passed
C to SORTXYZ. See also COMPARE2.
Cold  I11=I10+NATOMS
      I11=I10+MXATMS*2
      I12=I11+NATOMS
      I13=I12+NATOMS

#ifdef _USE_OPRTYP
      CALL CHRTAB(PTGRP,IORDER,NORDER,OPRTYP,ICORE(I2),
     &            ICORE(I3),ICORE(I4),ICORE(I5),NCLASS,
     &            ICORE(I6),ZCORE)
#else
      CALL CHRTAB(PTGRP,IORDER,NORDER,ICORE(I1),ICORE(I2),
     &            ICORE(I3),ICORE(I4),ICORE(I5),NCLASS,
     &            ICORE(I6),ZCORE)
#endif /* _USE_OPRTYP */
      CALL TSTOPS(IORDER,Q,ZCORE,ICORE(I6),ICORE(I8),
     &            ZCORE(IC2),ICORE(I10),NATOMS)
      CALL SYMUNQ(NATOMS,IORDER,IATNUM,ZCORE,ICORE(I8),
     &            ICORE(I9),ICORE(I10),ICORE(I11),ICORE(I12),
     &            IORBIT,PTGRP)
      if (label.eq.'COMP')
     &   CALL IPUTREC(20,'JOBARC','ZMAT2MOL',NATOMS,ICORE(I9))
      CALL FLUSHS(IORDER,IORBIT,ZCORE,ICORE(I8),ICORE(I10),
     &            ICORE(I12),ICORE(I6),NATOMS,LABEL,PTGRP,NCLASS)

csb
      if ((label.eq.'FULL').and.(iflags(46).gt.0)) then
         j=1
         write(*,*) '*** COORDINATES'
         do i=1,natoms
            write(*,'(x,i3,3(2x,f15.8))') i,q(j+0),q(j+1),q(j+2)
            j=j+3
         end do
         write(*,*) '*** END'
      end if

C CONVERT BACK TO JODA ORIENTATION IF THIS IS THE CALL FOR THE
C FULL POINT GROUP
      IF ((LABEL.EQ.'FULL').OR.(PTGRP.EQ.FPGRP))
     &   CALL STDORT(Q,ICORE,3*NATOMS,1)

      RETURN
      END

