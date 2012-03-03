
C THIS ROUTINE SERVES TWO FUNCTIONS.  FIRST, IT HANDLES THE VERY
C COMPLICATED SITUATIONS IN WHICH A DEGENERATE SET OF COORDINATES
C IS NOT RESOLVED INTO SEPARATE IRREPS OF THE ABELIAN SUBGROUP.

C THE ROUTINE ALSO DETERMINES IF ANY OPERATION IN THE POINT GROUP
C MAPS THE COORDINATE INTO ITS INVERSE.

c INPUT
c integer NATOM
c integer NIRREP
c integer IORDER
c integer NORBIT
c char*4  DOIT

c OUTPUT
c integer ISYMTYP(3*NATOM)
c double  SYMQ(9*NATOM*NATOM)
c double  SCR(9*NATOM)
c double  CHAR(IORDER,NIRREP)
c double  SYOP(9*IORDER)
c char*8  LABEL(NIRREP)
c integer ICENSUS(NATOM)
c integer IPTR(NATOM,IORDER)
c integer NBFATM(NATOM)
c integer ILCATM(NATOM)
c integer INVOP(3*NATOM)

c RECORDS
c get DOIT//'MEMB'
c get DOIT//'CHAR'
c get DOIT//'SYOP'
c get DOIT//'PERM'
c get DOIT//'LABL'
c get DOIT//'DEGN'
c get DOIT//'SYMQ'
c get 'ORIENTMT'
c get 'NVIBSYMF'
c get 'SBGRPSYM'
c put 'FULLSYMQ'
c put 'INVPSMAT'
c put 'SBGRPSYM'

      SUBROUTINE PRDEGEN(NATOM,NIRREP,IORDER,NORBIT,ISYMTYP,SYMQ,
     &                   SCR,CHAR,SYOP,LABEL,ICENSUS,IPTR,
     &                   NBFATM,ILCATM,INVOP,DOIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ISYMTYP(3*NATOM),SYMQ(9*NATOM*NATOM)
      DIMENSION SCR(9*NATOM),CHAR(IORDER,NIRREP),SYOP(9*IORDER)
      CHARACTER*8 LABEL(NIRREP)
      DIMENSION ICENSUS(NATOM),IPTR(NATOM,IORDER),NBFATM(NATOM)
      DIMENSION ILCATM(NATOM),INVOP(3*NATOM)
      CHARACTER*4 DOIT

      DIMENSION NVIBSYM(20),IDEGEN(20)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      DATA ONE /1.0/
      DATA TOL /1.D-8/

      NSIZE=3*NATOM

      CALL GETREC(9,'JOBARC',DOIT//'MEMB',NATOM,ICENSUS)
      CALL GETREC(9,'JOBARC',DOIT//'CHAR',NIRREP*IORDER*IINTFP,CHAR)
      CALL GETREC(9,'JOBARC',DOIT//'SYOP',9*IORDER*IINTFP,SYOP)
      CALL GETREC(9,'JOBARC',DOIT//'PERM',NATOM*IORDER,IPTR)
      CALL GETREC(9,'JOBARC',DOIT//'LABL',NIRREP*IINTFP,LABEL)
      CALL GETREC(9,'JOBARC',DOIT//'DEGN',NIRREP,IDEGEN)
      CALL GETREC(9,'JOBARC','ORIENTMT',9*IINTFP,SCR)
      CALL GETREC(9,'JOBARC','NVIBSYMF',NIRREP,NVIBSYM)
      CALL GETREC(9,'JOBARC',DOIT//'SYMQ',NSIZE*NSIZE*IINTFP,SYMQ)
      CALL GETREC(9,'JOBARC','SBGRPSYM',NSIZE,ISYMTYP)

      IF (DOIT.EQ.'FULL') CALL TRNOPS(SYOP,SCR,IORDER)

c   o fill out length and offset vectors for image()
      I = 1
      DO IATOM=1,NATOM
         NBFATM(IATOM)=3
         ILCATM(IATOM)=I
         I = I+3
      END DO

c   o loop over irreps
      IPOSABS=1
      DO 10 IRREP=1,NIRREP
         NVIBTOT=NVIBSYM(IRREP)

c      o continue for degenerate irreps
         IDEG=IDEGEN(IRREP)
         IF (IDEG.GT.1) THEN
            IPOSTAR=1

c         o continue if there is more than one vibration of this type
            NVIBUNQ=NVIBSYM(IRREP)/IDEG
            IF (NVIBUNQ.GT.1) THEN

c            o make sure the first mode of this symmetry is an irrep by itself
               IF (ISYMTYP(IPOSABS).EQ.ISYMTYP(IPOSABS+NVIBUNQ)) THEN
                  IPOS0=1+(IPOSABS-1)*NSIZE

C FOR DEGENERACY LEVELS >2, CHECK TO SEE IF THERE IS AT LEAST ONE
C COMPONENT WHICH BELONGS TO A UNIQUE SUBGROUP IRREP

                  IF (IDEG.GT.2) THEN
                     IUSE=0
                     IMODE=IPOSABS+2*NVIBUNQ
                     DO ICOUNT=3,IDEG
                        IF (IUSE.EQ.0) THEN
                           IF (ICOUNT.EQ.IDEG) THEN
                              IF (ISYMTYP(IMODE).NE.
     &                            ISYMTYP(IMODE-NVIBUNQ)) IUSE=IMODE
                           ELSE
                              IF (ISYMTYP(IMODE).NE.
     &                            ISYMTYP(IMODE-NVIBUNQ).AND.
     &                            ISYMTYP(IMODE).NE.
     &                            ISYMTYP(IMODE+NVIBUNQ)) IUSE=IMODE
                           END IF
                           IF (IUSE.NE.0) THEN
c                           o we are in luck. simply swap modes.
                              IPOSABS1=IMODE
                              IPOS1=1+(IPOSABS1-1)*NSIZE
                              CALL DSWAP(NVIBUNQ*NSIZE,SYMQ(IPOS0),1,
     &                                                 SYMQ(IPOS1),1)
                              DO I=0,NVIBUNQ-1
                                 ITMP=ISYMTYP(IPOSABS+I)
                                 ISYMTYP(IPOSABS+I)=ISYMTYP(IPOSABS1+I)
                                 ISYMTYP(IPOSABS1+I)=ITMP
                              END DO
                           END IF
                           IMODE=IMODE+NVIBUNQ
                        END IF
                     END DO
                     IPOSABS=IPOSABS+NVIBTOT
                  END IF

                  IF (IUSE.NE.0) GOTO 10

                  write(6,501)label(irrep)
501               format(t3,'@PRDEGEN-I, Degenerate modes of symmetry ',
     &                   A,' map to same irrep of subgroup.')

c               o find a self-adjoint symmetry operation with a character of 0
                  DO ISYMOP=1,IORDER
                  IF (ABS(CHAR(ISYMOP,IRREP)).LT.TOL) THEN
                     CALL XGEMM('N','T',3,3,3,
     &                          ONE, SYOP(1+(ISYMOP-1)*9),3,
     &                               SYOP(1+(ISYMOP-1)*9),3,
     &                          0.d0,SCR,                 3)
                     TRACE=DSUM(3,SCR,4)
                     IF (ABS(TRACE-3.d0).LT.TOL) IDO=ISYMOP
                  END IF
                  END DO

c               o project this manifold onto -1 eigenstates of sym op ISYMOP
                  IOFFSCR=3*NSIZE+1
                  IOFFSYQ=IPOS0
                  ITHRU=0
                  IOFF=IOFFSYQ
                  CHRUSE=-1.d0
                  DO IBAS=1,NVIBSYM(IRREP)
                     CALL ZERO(SCR,3*NATOM)
                     CALL IMAGE(NATOM,3*NATOM,1,IDO,IPTR,NBFATM,ILCATM,
     &                          SYMQ(IOFF),SCR,SCR(NSIZE+1),1,3*NATOM,
     &                          SYOP,1)
c                    write(6,*)' character is ',ddot(nsize,symq(ioff),1,scr,1)
                     CALL DSCAL(NSIZE,CHRUSE,SCR,1)
                     CALL VADD(SCR(2*NSIZE+1),SYMQ(IOFF),SCR,NSIZE,ONE)
                     IOFF=IOFF+NSIZE
                     ZLEN=DNRM2(3*NATOM,SCR(2*NSIZE+1),1)
                     IF (ZLEN.GT.TOL) THEN
                        XX=1.d0/ZLEN
                        CALL DCOPY(NSIZE,SCR(2*NSIZE+1),1,SCR,1)
                        CALL DSCAL(3*NATOM,XX,SCR,1)
                        CALL GSCHMIDT(SCR,SCR(3*NSIZE+1),3*NATOM,ITHRU,
     &                                SCR(2*NSIZE+1),RESID,TOL)
                        IF (RESID.GT.TOL) THEN
                           ITHRU=ITHRU+1
                           CALL DCOPY(NSIZE,SCR,1,SCR(IOFFSCR),1)
                           IOFFSCR=IOFFSCR+NSIZE
                        END IF
                     END IF
                  END DO

c               o form complementary space
                  IOFF=IOFFSYQ
                  DO IBAS=1,NVIBSYM(IRREP)
                     CALL GSCHMIDT(SYMQ(IOFF),SCR(3*NSIZE+1),3*NATOM,
     &                             ITHRU,SCR(2*NSIZE+1),RESID,TOL)
                     IF (RESID.GT.TOL) THEN
                        CALL DCOPY(NSIZE,SYMQ(IOFF),1,SCR(IOFFSCR),1)
                        IOFFSCR=IOFFSCR+NSIZE
                        ITHRU=ITHRU+1
                     END IF
                     IOFF=IOFF+NSIZE
                  END DO

                  CALL XCOPY(NSIZE*NVIBSYM(IRREP),SCR(3*NSIZE+1),1,
     &                                            SYMQ(IOFFSYQ),1)

               END IF
            END IF
         END IF
         IPOSABS=IPOSABS+NVIBTOT
10    CONTINUE

c   o determine which op of the full pt grp (if any)
c     maps this coord into its negative
      IOFF=1
      DO ICOORD=1,NSIZE
         INVOP(ICOORD)=0
         DO IOP=1,IORDER
            CALL IMAGE(NATOM,3*NATOM,1,IOP,IPTR,NBFATM,ILCATM,
     &                 SYMQ(IOFF),SCR,SCR(NSIZE+1),1,3*NATOM,SYOP,1)
            CALL DAXPY(NSIZE,ONE,SYMQ(IOFF),1,SCR,1)
            X=DNRM2(NSIZE,SCR,1)
            IF (X*X.LT.TOL) INVOP(ICOORD)=IOP
              
         END DO
         IOFF=IOFF+NSIZE
      END DO

      CALL PUTREC(20,'JOBARC','FULLSYMQ',NSIZE*NSIZE*IINTFP,SYMQ)
      CALL PUTREC(20,'JOBARC','INVPSMAT',NSIZE,INVOP)
      CALL PUTREC(20,'JOBARC','SBGRPSYM',NSIZE,ISYMTYP)

      RETURN
      END

