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
      SUBROUTINE READGH(IGRD,IHES,NATOM,G,H,IMAP,SCR,IAVGRD,IAVHES)
C
C READS IN THE GRADIENT AND HESSIAN MATRICES
C
C INPUT :
C
C   IGRD = -1 DON'T READ GRADIENT, RETURN G UNTOUCHED
C        =  0 READ GRADIENT, ASSUMING VMOL ORDER
C        =  1 READ GRADIENT, ASSUMING ZMAT ORDER (DUMMIES INCLUDED)
C   IHES -  HESSIAN READ CONTROL FLAG.  VALUES HAVE SAME MEANING AS 
C           ABOVE
C  NATOM -  THE NUMBER OF ATOMS
C
C OUTPUT :
C
C      G -  THE GRADIENT VECTOR
C IAVGRD -  =1 GRADIENT FOUND; =0 GRADIENT NOT FOUND
C      H -  THE HESSIAN MATRIX
C IAVHES -  =1 HESSIAN FOUND;  =0 HESSIAN NOT FOUND
C
C SCRATCH :
C
C    SCR - LENGTH 9*NATOM*NATOM
C   IMAP - LENGTH NATOM
C
C The vmol and seward both generates symm. redundant atoms in the same
C order. So the MAP2ZMAT remains unchanged. See vmol2ja (v2ja.f) for
C detailed notes. Ajith Perera 07/2000
C
CEND       
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL YESNO
      DIMENSION H(3*NATOM,3*NATOM),G(3,NATOM),SCR(9*NATOM*NATOM)
      DIMENSION IMAP(NATOM)
C
      PARAMETER (TOL = 1.D-8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
C
      IONE=1
      IAVGRD=0
      IAVHES=0
      NSIZE=3*NATOM
      CALL IGETREC(20,'JOBARC','NREALATM',IONE,NREAL)
      CALL IGETREC(20,'JOBARC','MAP2ZMAT',NATOM,IMAP)
C
      IF(IGRD.GE.0)THEN
C
C READ IN GRADIENT
C 
       CALL ZERO(G,NSIZE)
       NREAD=NREAL
       IF(IGRD.EQ.1)NREAD=NATOM
       CALL DGETREC(-1,'JOBARC','GRADIENT',3*NREAD,SCR)
       X=XDNRM2(3*NREAD,SCR,1)
       IF(X.GT.TOL)IAVGRD=1
C
       IF(IGRD.EQ.0)THEN
C
C ASSUME VMOL ORDERING
C
         IOFF=1
          DO 10 IATMVML=1,NATOM
           IATMZMAT=IMAP(IATMVML)
           IF(IATMZMAT.NE.0)THEN
           CALL XDCOPY(3,SCR(IOFF),1,G(1,IATMZMAT),1) 
c           G(1,IATMZMAT) = SCR(IOFF)
c           G(2,IATMZMAT) = SCR(IOFF+1)
c           G(3,IATMZMAT) = SCR(IOFF+2)
           IOFF=IOFF+3
           ENDIF
10        CONTINUE
       ELSE
C
C ASSUME ZMAT ORDERING
C
        CALL XDCOPY(NSIZE,SCR,1,G,1)
       ENDIF
      ENDIF
C
      IF(IHES.GE.0)THEN
C
C READ IN HESSIAN
C
       NGET=3*NREAL
       IF(IHES.EQ.1)NGET=NSIZE
       CALL ZERO(H,NSIZE*NSIZE)
       CALL DGETREC(-1,'JOBARC','HESSIANM',NGET*NGET,SCR)
       X=XDNRM2(9*NREAL*NREAL,SCR,1)
       IF(X.GT.TOL)IAVHES=1
       IF(IFLAGS2(3).EQ.1)THEN
        INQUIRE(FILE='FCMFINAL',EXIST=YESNO)
        IF(.NOT.YESNO)THEN
         WRITE(6,200)
         CALL ERREX
        ENDIF
        OPEN(UNIT=15,FILE='FCMFINAL')
        READ(15,*)
        READ(15,'((3F20.10))')SCR
        CLOSE(UNIT=15, STATUS="KEEP")
       ENDIF
C
       IF(IHES.EQ.0)THEN
C
C ASSUME VMOL ORDERING
C
        IO=1
        DO 110 IATVML=1,NATOM
         IATMZMAT=IMAP(IATVML)
         DO 115 IXYZ=1,3
          ICOL =IXYZ+(IATMZMAT-1)*3
          DO 120 JATVML=1,NATOM
           JATMZMAT=IMAP(JATVML)
           IF(IATMZMAT.NE.0.AND.JATMZMAT.NE.0)THEN
            IROW=1+(JATMZMAT-1)*3
            CALL BLKCPY(SCR(IO),3,1,H,NSIZE,NSIZE,IROW,ICOL)
            IO=IO+3
           ENDIF
120       CONTINUE
115      CONTINUE
110     CONTINUE
C
C ASSUME ZMAT ORDERING
C
       ELSE
        CALL XDCOPY(NSIZE*NSIZE,SCR,1,H,1)
       ENDIF
      ENDIF
C
      OPEN(UNIT=13,FILE='FCMFINAL',FORM='FORMATTED')
      WRITE(13,'(2I5)')NATOM,3*NATOM
      WRITE(13,'((3F20.10))')H
      CLOSE(UNIT=13,STATUS='KEEP')
C
200   FORMAT(T3,'@READGH-F, Resonance Raman intensities require ',
     &       'reference state FCMFINAL file.',/,
     &       T3,'Save JOBARC and JAINDX file and run ',
     &       'xjoda twice with FCMFINAL file in place.')
C
      RETURN
      END
