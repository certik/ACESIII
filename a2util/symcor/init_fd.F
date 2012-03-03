
c This routine performs all the symmetry analysis needed to generate the grid
c of displacements. It is essentially a method that uses and updates JOBARC
c records. NENER is returned solely to tell symcor that there are points to
c start working on (or not).

c INPUT
c char*4  DOIT
c integer NATOM
c integer ICRSIZ
c logical BPRINT

c OUTPUT
c integer ICORE(ICRSIZ)
c integer NENER

c RECORDS
c get DOIT//'ORDR'
c get DOIT//'NORB'
c get DOIT//'NIRX'
c get 'COMPORDR'
c get 'COMPNORB'
c get 'COMPNIRX'

#include "flags.h"

      subroutine init_fd(doit,natom,icore,icrsiz,nener,bprint)
      implicit none

      character*4 doit
      integer natom, icore(*), icrsiz, nener
      logical bprint

      integer iorder, iorderf, iorderc
      integer norbit, norbitf, norbitc
      integer nirrep, nirrepf, nirrepc
      integer I000, I010, I020, I030, I040, I050, I060, I070, I080, I090
      integer I100, I110, I120, I130, I140, I150, I160, I170, I180, I190
      integer I200, I0
      character*8 irrnm(32), label(32)

#include "control.com" /* for ROTPROJ */
#include "machsp.com"

      CALL GETREC(20,'JOBARC',DOIT//'ORDR',1,IORDERF)
      CALL GETREC(20,'JOBARC',DOIT//'NORB',1,NORBITF)
      CALL GETREC(20,'JOBARC','COMPORDR',1,IORDERC)
      CALL GETREC(20,'JOBARC','COMPNORB',1,NORBITC)
      IORDER=MAX(IORDERF,IORDERC)
      NORBIT=MAX(NORBITF,NORBITC)

      I0 = 1

      I000=I0
      I010=I000 + IINTFP*IORDER*9
      I020=I010 + IINTFP*IORDER*IORDER*9
      I030=I020 + IINTFP*IORDER*IORDER
      I040=I030 + IINTFP*IORDER*IORDER
      I050=I040 + 3*IORDER
      I060=I050 + IORDER*IORDER
      I080=I060 + IORDER + IAND(IORDER,1)
      I090=I080 + 6*IORDER
      I100=I090 + IINTFP*3*IORDER

      CALL CHRTABLE(IORDERF,ICORE(I000),ICORE(I010),ICORE(I020),
     &              ICORE(I030),ICORE(I040),ICORE(I050),ICORE(I060),
     &              IRRNM,ICORE(I080),DOIT,ICORE(I090))
      CALL CHRTABLE(IORDERC,ICORE(I000),ICORE(I010),ICORE(I020),
     &              ICORE(I030),ICORE(I040),ICORE(I050),ICORE(I060),
     &              IRRNM,ICORE(I080),'COMP',ICORE(I090))

      CALL GETREC(-1,'JOBARC',DOIT//'NIRX',1,NIRREPF)
      CALL GETREC(-1,'JOBARC','COMPNIRX',1,NIRREPC)
      NIRREP=MAX(NIRREPF,NIRREPC)

      I010=I000 + IINTFP*NATOM*NATOM*9
      I020=I010 + IINTFP*NATOM*NATOM*9
      I030=I020 + IINTFP*NATOM*NATOM*9
      I040=I030 + IINTFP*MAX(3*NATOM,NIRREPF*IORDER)
      I050=I040 + IINTFP*MAX(6*NATOM,9*IORDER)
      I060=I050 + IINTFP*MAX(6*NATOM,NIRREPF)
      I060=I050
      I070=I060 + IINTFP*9*NATOM*NATOM*3*NATOM
      I080=I070 + MAX(MAX(2,IORDER)*NATOM,9*NATOM*NATOM)
      I090=I080 + NATOM
      I100=I090 + NATOM + IAND(NATOM,1)
      I110=I100 + IINTFP*NATOM*NATOM*MAX(18,IORDER)
      I120=I110 + IINTFP*NATOM*NATOM*18
      I130=I120 + IINTFP*NATOM*NATOM*18
      BPRINT =.TRUE.

      CALL VIBINF(NATOM,NIRREPF,IORDERF,ICORE(I000),
     &            ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &            LABEL,ICORE(I060),ICORE(I070),ICORE(I080),
     &            ICORE(I090),ICORE(I100),DOIT)
      CALL SYMADQ(NATOM,NIRREPC,IORDERC,NORBITC,ICORE(I000),
     &            ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &            LABEL,ICORE(I060),ICORE(I070),ICORE(I080),
     &            ICORE(I090),ICORE(I100),ICORE(I110),ICORE(I120),
     &            ICORE(I130),'COMP',1)

      CALL TRAPRJ(NATOM,ICORE(I000),ICORE(I010),ICORE(I020),
     &            ICORE(I030))

      IF (ROTPROJ) THEN
         CALL ROTPRJ(NATOM,ICORE(I000),ICORE(I010),ICORE(I020),
     &              ICORE(I030))
      END IF

      CALL SCHMIDT(NATOM,ICORE(I010),ICORE(I020))

      CALL COLLECT(NATOM,NIRREPC,ICORE(I000),
     &             ICORE(I010),ICORE(I020),ICORE(I030),LABEL,
     &             'COMP')

      IF (BPRINT) THEN
         CALL PRTCOORD(NATOM,NIRREPC,0,ICORE(I000),ICORE(I010),
     &                 LABEL,'COMP',.FALSE.)
      END IF

      CALL SYMADQ(NATOM,NIRREPF,IORDERF,NORBITF,ICORE(I000),
     &            ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &            LABEL,ICORE(I060),ICORE(I070),ICORE(I080),
     &            ICORE(I090),ICORE(I100),ICORE(I110),ICORE(I120),
     &            ICORE(I130),DOIT,2)
      CALL PRDEGEN(NATOM,NIRREPF,IORDERF,NORBITF,ICORE(I000),
     &             ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &             LABEL,ICORE(I060),ICORE(I070),ICORE(I080),
     &             ICORE(I090),ICORE(I100),DOIT)

      IF (BPRINT) THEN
         CALL PRTCOORD(NATOM,NIRREPF,NIRREPC,ICORE(I000),ICORE(I010),
     &                 LABEL,DOIT,.TRUE.)
      END IF

      I010=I000
      I020=I010 + 3*NATOM
      I030=I020 + 9*NATOM*NATOM
      I040=I030 + 3*NATOM + IAND(NATOM,1)
      I200=I040 + IINTFP*9*NATOM*NATOM*3*NATOM
#ifdef _DEBUG_LVL0
      Print*, "Entering set points"
#endif
      CALL SETPTS(NATOM,NIRREPF,DOIT,
     &            NENER,IRRNM,ICORE(I010),ICORE(I020),ICORE(I030),
     &            ICORE(I040),ICORE(I200),(ICRSIZ-I200+I0)/IINTFP)

      return
c     end subroutine init_fd
      end

