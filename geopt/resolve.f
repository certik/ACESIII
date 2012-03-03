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
      SUBROUTINE RESOLVE(NDEG,V,TMP,SYMOPS,IPTR,TMP2,BIGX,NX,NATOM,
     &                   NTOTATOM,IORDABEL)
C
C RESOLVES DEGENERATE NORMAL MODE INTO PURE IRREDUCIBLE REPRESENTATIONS
C OF THE ABELIAN SUBGROUP
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IPRINT,IPRINT2
      DIMENSION V(NX,NDEG),TMP(NX,NDEG),SYMOPS(9*IORDABEL),
     &          IPTR(NATOM*IORDABEL),BIGX(NX,NX),TMP2(NX)
      DIMENSION CHAR2(2,2),CHAR4(4,4),CHAR8(8,8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/ORIENT/ORIENT(9)
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     &   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     &   XYZTol
      DATA ONE,ONEM,ZILCH/1.0D0,-1.0d0,0.0D0/
      DATA CHAR2/ 1.0D0, 1.0D0,
     &           -1.0D0, 1.0D0/
      DATA CHAR4/ 1.D0, 1.D0, 1.D0, 1.D0,
     &           -1.D0,-1.D0, 1.D0, 1.D0,
     &           -1.D0, 1.D0,-1.D0, 1.D0,
     &            1.D0,-1.D0,-1.D0, 1.D0/
      DATA CHAR8/ 1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0,
     &            1.D0,-1.D0, 1.D0, 1.D0,-1.D0,-1.D0,-1.D0, 1.D0,
     &            1.D0, 1.D0,-1.D0,-1.D0,-1.D0,-1.D0, 1.D0, 1.D0,
     &           -1.D0, 1.D0,-1.D0, 1.D0, 1.D0,-1.D0,-1.D0, 1.D0,
     &           -1.D0, 1.D0, 1.D0,-1.D0,-1.D0, 1.D0,-1.D0, 1.D0,
     &           -1.D0,-1.D0, 1.D0,-1.D0, 1.D0,-1.D0, 1.D0, 1.D0,
     &           -1.D0,-1.D0,-1.D0, 1.D0,-1.D0, 1.D0, 1.D0, 1.D0,
     &            1.D0,-1.D0,-1.D0,-1.D0, 1.D0, 1.D0,-1.D0, 1.D0/
C
      iprnt=1
      IPRINT2=IPRNT.GT.100
      IPRINT =IPRNT.GT.1
C
C CALCULATE MASS-WEIGHTED HESSIAN IN BASIS OF DEGENERATE EIGENVECTORS
C
      IF(IPRINT)THEN
       rewind(78)
       read(78)bigx
       call xgemm('n','n',nx,ndeg,nx,one,bigx,nx,v,nx,zilch,tmp,nx)
       call xgemm('t','n',ndeg,ndeg,nx,one,v,nx,tmp,nx,zilch,tmp2,ndeg)
       write(6,*)' Mass-weighted Hessian in unresolved coordinates:'
       call output(tmp2,1,ndeg,1,ndeg,ndeg,ndeg,1)
      ENDIF
      IF(IPRINT)THEN
       write(6,*)' Degenerate modes before resolution: '
       do 901 i=1,ndeg
        write(6,'((3f20.14))')(v(j,i),j=1,nx)
        write(6,*)' normalization ',xdnrm2(nx,v(1,i),1)
        write(6,*)
901    continue
      ENDIF
c
      CALL DGETREC(20,'JOBARC','COMPSYOP',9*IORDABEL,SYMOPS)
      CALL MTRAN2(ORIENT,3)
      CALL TRNOPS(SYMOPS,ORIENT,IORDABEL)
      CALL MTRAN2(ORIENT,3)
      CALL IGETREC(20,'JOBARC','COMPPERM',NTOTATOM*IORDABEL,IPTR)
      CALL STRIPDUM(IPTR,TMP,NATOM,NTOTATOM,IORDABEL,BIGX)
C
      CALL XDCOPY(NX*NDEG,V,1,TMP,1)
      DO 5 I=2,NDEG
       CALL XDAXPY(NX,ONE,TMP(1,I),1,TMP,1)
5     CONTINUE
      IF(IPRINT2)THEN
       write(6,*)' reducible representation '
       write(6,'((3f20.14))')(tmp(j,1),j=1,nx)
      ENDIF
C
C LOOP OVER SYMMETRY OPERATIONS AND PROJECT ONTO IRREPS
C
      IDONE=1
      DO 10 IRREP=1,IORDABEL
       IOFFS=1
       IOFFI=1
       CALL ZERO(TMP2,NX)
       DO 11 IOP=1,IORDABEL
        CALL GENREP(BIGX,SYMOPS(IOFFS),IPTR(IOFFI),NATOM)
        IF(IORDABEL.EQ.2.AND.IRREP.EQ.2)THEN
         X=SSUM(9,SYMOPS(1),1)
         FACT=ONE*CHAR2(IOP,IRREP)
         IF(DABS(X-3.0D0).LT.1.D-3)FACT=ONEM*CHAR2(IOP,IRREP)
        ELSEIF(IORDABEL.EQ.2.AND.IRREP.EQ.1)THEN
         FACT=CHAR2(IOP,IRREP)
        ELSEIF(IORDABEL.EQ.4)THEN
         FACT=CHAR4(IOP,IRREP)
        ELSEIF(IORDABEL.EQ.8)THEN
         FACT=CHAR8(IOP,IRREP)
        ENDIF
        CALL XGEMM ('N','N',NX,1,NX,FACT,BIGX,NX,TMP,NX,ONE,TMP2,NX) 
        IOFFS=IOFFS+9
        IOFFI=IOFFI+NTOTATOM
11     CONTINUE
       XLEN=XDNRM2(NX,TMP2,1)
      IF(IPRINT2)THEN
       write(6,*)' after projection onto irrep ',irrep
       write(6,'((3f20.10))')tmp2
       write(6,*)' length of projection for irrep ',irrep,' is ',xlen
      ENDIF
       if(xlen.gt.1.d-2)then
        call xdscal(nx,1.0d0/dfloat(iordabel),tmp2,1)
        call xdcopy(nx,tmp2,1,v(1,idone),1)
        idone=idone+1
       endif
10    CONTINUE
      if(idone.ne.ndeg+1)then
       write(6,*)' Resolution is impossible ! '
       call dgetrec(20,'JOBARC','HOLDDEGQ',NDEG*NX,V)
      endif
      CALL ZERO(TMP,NX)
      IF(IPRINT)WRITE(6,*)' Degenerate modes after resolution and',
     &                    ' renormalization: '
       do 900 i=1,ndeg
        call xdaxpy(nx,one,v(1,i),1,tmp,1)
        IF(IPRINT2)write(6,*)' normalization ',xdnrm2(nx,v(1,i),1)
        call dscal(nx,1.D0/xdnrm2(nx,v(1,i),1),v(1,i),1)
        IF(IPRINT)write(6,'((3f20.14))')(v(j,i),j=1,nx)
        IF(IPRINT2)write(6,*)' renormalized norm ',xdnrm2(nx,v(1,i),1)
        IF(IPRINT)write(6,*)
900    continue
      IF(IPRINT2)write(6,*)' reducible representation on exit '
      IF(IPRINT2)write(6,'((3f20.14))')(tmp(j,1),j=1,nx)
C
C CALCULATE MASS-WEIGHTED HESSIAN IN BASIS OF DEGENERATE EIGENVECTORS
C
      IF(IPRINT)THEN
       rewind(78)
       read(78)bigx
       call xgemm('n','n',nx,ndeg,nx,one,bigx,nx,v,nx,zilch,tmp,nx)
       call xgemm('t','n',ndeg,ndeg,nx,one,v,nx,tmp,nx,zilch,tmp2,ndeg)
       call output(tmp2,1,ndeg,1,ndeg,ndeg,ndeg,1)
       write(6,*)' Mass-weighted Hessian in resolved coordinates:'
       call output(tmp2,1,ndeg,1,ndeg,ndeg,ndeg,1)
      ENDIF
      RETURN
      END
