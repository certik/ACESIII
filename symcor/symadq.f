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
c char*4  TYPE
c integer NPASS

c OUTPUT
c double  VGEN(NATOM*3)
c double  VREF(NATOM*3)
c double  VIMAGE(NATOM*3)
c double  CHAR(IORDER,NIRREP)
c double  SYOP(9*IORDER)
c char*8  LABEL(NIRREP)
c integer ICENSUS(NATOM)
c integer IPTR(NATOM,IORDER)
c integer NBFATM(NATOM)
c integer ILCATM(NATOM)
c double  SCR(*)
c double  SYMQ(*)
c integer ISYMTYP(*)
c integer INVOP(3*NATOM)

c RECORDS
c get TYPE//'MEMB'
c get TYPE//'CHAR'
c get TYPE//'SYOP'
c get TYPE//'PERM'
c get TYPE//'POPV'
c get TYPE//'LABL'
c get 'ORIENTMT'
c get 'NUMVIBRT'
c get 'COMPSYMQ'
c get 'COMPSYQT'
c put 'SBGRPSYM'
c put 'ORDERREF'
c put 'OPERSREF'
c put 'NVIBSYMF'
c put TYPE//'NSYQ'
c put TYPE//'SYQT'
c put TYPE//'SYMQ'

      SUBROUTINE SYMADQ(NATOM,NIRREP,IORDER,NORBIT,VGEN,VREF,
     &                  VIMAGE,CHAR,SYOP,LABEL,ICENSUS,IPTR,
     &                  NBFATM,ILCATM,SCR,SYMQ,ISYMTYP,
     &                  INVOP,TYPE,NPASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION VGEN(NATOM*3),VREF(NATOM*3),VIMAGE(NATOM*3)
      DIMENSION CHAR(IORDER,NIRREP),SYOP(9*IORDER)
      CHARACTER*8 LABEL(NIRREP)
      DIMENSION ICENSUS(NATOM),IPTR(NATOM,IORDER),NBFATM(NATOM)
      DIMENSION ILCATM(NATOM),SCR(*),SYMQ(*),ISYMTYP(*)
      DIMENSION INVOP(3*NATOM)
      CHARACTER*4 TYPE

      DIMENSION IPOP(255),buf(1000),NVIBSYM(100)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)

      DATA TOL /1.D-5/
      DATA ONE /1.0/

      IPOS=1
      NSYMOLD=0
      NSIZE=3*NATOM

      CALL ZERO(VREF,NSIZE)
      CALL ZERO(SYMQ,NSIZE*NSIZE)
      CALL IZERO(INVOP,NSIZE)

      CALL IGETREC(20,'JOBARC',TYPE//'MEMB',NATOM,ICENSUS)
      CALL DGETREC(20,'JOBARC',TYPE//'CHAR',NIRREP*IORDER,CHAR)
      CALL DGETREC(20,'JOBARC',TYPE//'SYOP',9*IORDER,SYOP)
      CALL IGETREC(20,'JOBARC',TYPE//'PERM',NATOM*IORDER,IPTR)
      CALL IGETREC(20,'JOBARC',TYPE//'POPV',NORBIT,IPOP)
      CALL DGETREC(20,'JOBARC',TYPE//'LABL',NIRREP,LABEL)
      CALL DGETREC(20,'JOBARC','ORIENTMT',9,SCR)

c   o fill basis vector skipping dummy atoms
      IF (NPASS.EQ.2) THEN
         IF (TYPE.EQ.'FULL') CALL TRNOPS(SYOP,SCR,IORDER)
         CALL IGETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
         CALL DGETREC(20,'JOBARC','COMPSYMQ',NSIZE*NSIZE,SCR)
         CALL IGETREC(20,'JOBARC','COMPSYQT',NSIZE,ISYMTYP(NSIZE+1))
         NTOP=NMODE
      ELSE
         IOFF=1
         CALL ZERO(SCR,NSIZE*NSIZE)
         IOFFS=0
         DO IORBIT=1,NORBIT
            NPOP=IPOP(IORBIT)
            DO I=1,NPOP
               J=ICENSUS(IOFF)
               IOFF=IOFF+1
               DO IXYZ=1,3
                  SCR(IOFFS+IXYZ+(J-1)*3)=ONE
                  IOFFS=IOFFS+NSIZE
               END DO
            END DO
         END DO
         NTOP=NSIZE
      END IF

c   o fill out length and offset vectors for image()
      I = 1
      DO IATOM=1,NATOM
         NBFATM(IATOM)=3
         ILCATM(IATOM)=I
         I = I+3
      END DO

      ZNORM=ONE/DFLOAT(IORDER)

c   o loop over the symmetry blocks
      DO IRREP=1,NIRREP

c      o loop over basis functions
         IOFF=1
         NVIBSYM(IRREP)=0
         DO IBAS=1,NTOP
            CALL XDCOPY(NSIZE,SCR(IOFF),1,VREF,1)
            CALL ZERO(VGEN,3*NATOM)

C NOW LOOP OVER ALL SYMMETRY OPERATIONS AND COMPUTE IMAGE OF
C REFERENCE DIRECTION AND MULTIPLY RESULT BY CHARACTER OF OPERATION

            DO IOP=1,IORDER
               CALL IMAGE(NATOM,3*NATOM,1,IOP,IPTR,NBFATM,ILCATM,
     &                    VREF,VIMAGE,SCR(NSIZE*NSIZE+1),1,3*NATOM,
     &                    SYOP,1)
               CALL XDSCAL(3*NATOM,CHAR(IOP,IRREP),VIMAGE,1)
               CALL XDAXPY(3*NATOM,ONE,VIMAGE,1,VGEN,1)
            END DO
            ZLEN=DNRM2(3*NATOM,VGEN,1)
            IF (ZLEN.GT.TOL) THEN

c            o normalize coordinates
               XX=ONE/ZLEN
               CALL XDSCAL(3*NATOM,xx,VGEN,1)

c            o schmidt orthogonalize to all previous vectors
               CALL GSCHMIDT(VGEN,SYMQ,3*NATOM,NSYMOLD,BUF,ZJUNK,1.D-8)

               ZLEN2=DNRM2(3*NATOM,VGEN,1)
               IF (ZLEN2.GT.TOL) THEN
                  NSYMOLD=NSYMOLD+1
                  ISYMTYP(NSYMOLD)=IRREP
                  ZNORM=1.D0/ZLEN
                  CALL XDCOPY(3*NATOM,VGEN,1,SYMQ(IPOS),1)
                  IPOS=IPOS+3*NATOM
                  NVIBSYM(IRREP)=NVIBSYM(IRREP)+1

C FOR FULL POINT GROUP, WE NEED TO KNOW THE SYMMETRY OF THE GENERATOR
C IN THE ABELIAN SUBGROUP.  WE ALSO WANT TO KNOW WHICH OPERATION IN THE
C FULL POINT GROUP MAPS THIS COORDINATE INTO ITS NEGATIVE (NONSYMMETRIC
C DISPLACEMENTS ONLY)

                  IF (NPASS.EQ.2) THEN
                     ISYMTYP(2*NSIZE+NSYMOLD)=ISYMTYP(NSIZE+IBAS)
                  END IF
                  IF (IFLAGS(1).GT.100) THEN
                     write(6,*)'  species ',label(irrep)(1:4)
                     ioffx=0
                     do i=1,natom
                        write(6,'(i5,3f20.10)')i,(vgen(ioffx+ix),ix=1,3)
                        ioffx=ioffx+3
                     end do
                  END IF

               END IF
            END IF
            IOFF=IOFF+NSIZE
         END DO
      END DO

      IF (NPASS.EQ.2) THEN
         CALL XDCOPY(NTOP*NSIZE,SYMQ,1,SCR,1)
         CALL XDCOPY(NSIZE*NSIZE,SCR,1,SYMQ,1)
         CALL IPUTREC(20,'JOBARC','SBGRPSYM',NSIZE,ISYMTYP(2*NSIZE+1))
         CALL IPUTREC(20,'JOBARC','ORDERREF',1,IORDER)
         CALL DPUTREC(20,'JOBARC','OPERSREF',IORDER*9,SYOP)
         CALL IPUTREC(20,'JOBARC','NVIBSYMF',NIRREP,NVIBSYM)
      END IF
      CALL IPUTREC(20,'JOBARC',TYPE//'NSYQ',1,NSYMOLD)
      CALL IPUTREC(20,'JOBARC',TYPE//'SYQT',NSIZE,ISYMTYP)
      CALL DPUTREC(20,'JOBARC',TYPE//'SYMQ',NSIZE*NSIZE,SYMQ)

      RETURN
      END

