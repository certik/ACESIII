      SUBROUTINE A2_PREP4_MOL_READ(ICORE, MAXCOR, NATOMS, NBASP, 
     &                             NBAS, IUCATMS, IUATMS, COORD, 
     &                             IATMCHRG, Norbits_fullG, 
     &                             Norbits_compG,MAXSHELL, MAXANG,
     &                             ITFCT, NUNQSHL, NTOTSHL, NTANGM, 
     &                             LNP1, INext)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION COORD(3*NATOMS), Norbits_fullG(IUATMS), 
     &          Norbits_compG(IUCATMS), ICORE(MAXCOR), 
     &          IATMCHRG(NATOMS)
C     
#include "flags.h"
#include "machsp.com"
#include "jodaflags.com"
#include "syminf.com"
#include "offsets_2getden.com"
    
      DATA IONE /1/
#ifdef _DEBUG_LVLM2
      print*, "input variables@ prep4_dens_plots"
      print*, NATOMS, IUCATMS, IUATMS,  MAXSHELL
      Write(*,*)
#endif
      IMEMBC       = IONE 
      INUC         = IMEMBC       + NATOMS
      NUFCT        = INUC         + NATOMS
      NUMOM        = NUFCT        + IUCATMS 
      NFCT         = NUMOM        + IUCATMS 
      NANGMOM      = NFCT         + NATOMS
      IATMNAM      = NANGMOM      + NATOMS
      NAOUATM      = IATMNAM      + IUCATMS*IINTFP
      NAOATM       = NAOUATM      + IUCATMS 
      IANGMOMSHL   = NAOATM       + NATOMS
      ICONFUNSHL   = IANGMOMSHL   + MAXSHELL*IUCATMS
      IPRMFUNSHL   = ICONFUNSHL   + MAXSHELL*IUCATMS
      INPRIMSHL    = IPRMFUNSHL   + MAXSHELL*IUCATMS
      IANGMOMTSHL  = INPRIMSHL    + NATOMS 
      ICONFUNTSHL  = IANGMOMTSHL  + NATOMS*MAXSHELL
      IOFFSETPRM   = ICONFUNTSHL  + NATOMS*MAXSHELL
      IOFFSETCON   = IOFFSETPRM   + IUCATMS*MAXSHELL
      IOFFSETSHL   = IOFFSETCON   + IUCATMS*MAXSHELL
      IMAP_SHL2CNT = IOFFSETSHL   + NATOMS*MAXSHELL
      IPRMFUNTSHL  = IMAP_SHL2CNT + NATOMS*MAXSHELL 
      ISHLOFF_TMP  = IPRMFUNTSHL  + NATOMS*MAXSHELL  
      ISHLOFF      = ISHLOFF_TMP  + NATOMS
      IEND1        = ISHLOFF      + NATOMS
C    
      IF (IEND1 .GE. MAXCOR) CALL INSMEM('P4_DEN_PLOTS', IEND1,
     &                                    MAXCOR)
C
      CALL GETREC(20,'JOBARC', 'COMPMEMB', NATOMS, ICORE(IMEMBC))
      IPRINT = iflags(h_IFLAGS_print)
C     
#ifdef _DEBUG_LVL0
      print*, "Entering A2RD_BASIS_4INTGRT: Orbits and Cart. Coord."
      Print*, (Norbits_compG(i), i=1, iucatms)
      Write(*,*)
      Write(6, "(3F10.5)"), (coord(i), i=1,3*natoms)
      Write(*,*)
#endif
      CALL A2RD_BASIS_4INTGRT(IUCATMS,NATOMS,ITFCT,LNP1,LNPO,NTANGM,
     &                        ICORE(IMEMBC),ICORE(INUC),ICORE(NFCT),
     &                        ICORE(NUFCT),ICORE(NANGMOM), 
     &                        ICORE(NUMOM),ICORE(IATMNAM),
     &                        COORD,Norbits_compG,
     &                        ICORE(NAOATM),ICORE(NAOUATM),
     &                        ICORE(IANGMOMSHL),ICORE(ICONFUNSHL),
     &                        ICORE(IPRMFUNSHL),IPRINT,ISHL, 
     &                        ICORE(INPRIMSHL),ICORE(IANGMOMTSHL),
     &                        ICORE(ICONFUNTSHL),ICORE(IOFFSETPRM),
     &                        ICORE(IOFFSETCON),ICORE(IOFFSETSHL),
     &                        ICORE(ISHLOFF_TMP),ICORE(ISHLOFF),
     &                        MAXSHELL,NUNQSHL,NTOTSHL,IATMCHRG)
C
      NMOMFCT  = IEND1 
      NMOMAO   = NMOMFCT  + NATOMS*NTANGM 
      IALPHA   = NMOMAO   + NATOMS*NTANGM 
      ICOEFFA  = IALPHA   + ITFCT*IINTFP
      ICOEFFB  = ICOEFFA  + NBASP*NBAS*IINTFP
      IDENSA   = ICOEFFB  + NBASP*NBAS*IINTFP
      IDENSB   = IDENSA   + NBAS*NBAS*IINTFP
      IPCOEFFA = IDENSB   + NBAS*NBAS*IINTFP
      IPCOEFFB = IPCOEFFA + NBASP*ITFCT*IINTFP
      IPCOEFF  = IPCOEFFB + NBASP*ITFCT*IINTFP
      ISCR1    = IPCOEFF  + ITFCT*IINTFP*NBAS
      ISCR2    = ISCR1    + LNP1*IINTFP
      ISTOP    = ISCR2    + LNPO*IINTFP
      INEXT    = ISCR1
C     
      IF(ISTOP.GE.MAXCOR) CALL INSMEM('P4_DEN_PLOTS',ISTOP,MAXCOR)
#ifdef _DEBUG_LVL0
      Write(*,*)
      print*, "Entering A2RD_BASIS_PRIMINF", ITFCT, NTANGM, IALPHA,
     &          LNP1,NUNQSHL

#endif
C
#ifdef _DEBUG_LVLM2
      Print*, "Offsets for various arrays initilized in Pre4_Den_plo."
      Write(*, '(31I5)'), IMEMBC,INUC,NUFCT,NUMOM,NFCT,NANGMOM,
     &                   IATMNAM,NAOUATM,NAOATM,IANGMOMSHL,
     &                   ICONFUNSHL,IPRMFUNSHL,INPRIMSHL,
     &                   IANGMOMTSHL,ICONFUNTSHL,IOFFSETPRM,
     &                   IOFFSETCON,IOFFSETSHL,IMAP_SHL2CNT,
     &                   IPRMFUNTSHL,NMOMAO,IALPHA,ICOEFFA,
     &                   ICOEFFB,IDENSA,IDENSB,IPCOEFFA,
     &                   IPCOEFFB, IPCOEFF,ISCR1,ISCR2
#endif
      CALL A2RD_BASIS_PRIMINF(NATOMS,IUCATMS,ITFCT,NBAS,LNP1,LNPO,
     &                        NTANGM,Norbits_compG,ICORE(NFCT),
     &                        ICORE(NANGMOM),ICORE(NMOMFCT),
     &                        ICORE(NMOMAO),ICORE(IMEMBC),
     &                        ICORE(IANGMOMSHL),ICORE(IPRMFUNSHL),
     &                        ICORE(ICONFUNSHL),ICORE(ICONFUNTSHL), 
     &                        ICORE(INPRIMSHL),ICORE(IMAP_SHL2CNT),
     &                        ICORE(IPRMFUNTSHL),ICORE(IALPHA),
     &                        ICORE(IPCOEFF),ICORE(NAOATM),
     &                        ICORE(ISCR1),ICORE(ISCR2),
     &                        ICORE(ISHLOFF),MAXSHELL,ISHL,MAXANG)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The number of contracted function per shell"
      Write(6,*) (ICORE(ICONFUNTSHL - 1 + I), I=1, NTOTSHL)
      Write(6,*) "The shell angular momentum"
      Write(6,*) (ICORE(IANGMOMTSHL - 1 + I), I=1, NTOTSHL)
      Write(6,*) "The number of primitive functions per shell"
      Write(6,*) (ICORE(IPRMFUNTSHL - 1 + I), I=1, NTOTSHL)
#endif
#ifdef _DEBUG_LVL0
      Print*, "Offsets for various arrays initilized in Pre4_Den_plo."
      Write(*, '(31I5)'), IMEMBC,INUC,NUFCT,NUMOM,NFCT,NANGMOM,
     &                   IATMNAM,NAOUATM,NAOATM,IANGMOMSHL,
     &                   ICONFUNSHL,IPRMFUNSHL,INPRIMSHL,
     &                   IANGMOMTSHL,ICONFUNTSHL,IOFFSETPRM,
     &                   IOFFSETCON,IOFFSETSHL,IMAP_SHL2CNT,
     &                   IPRMFUNTSHL,NMOMAO,IALPHA,ICOEFFA,
     &                   ICOEFFB,IDENSA,IDENSB,IPCOEFFA,
     &                   IPCOEFFB, IPCOEFF,ISCR1,ISCR2
      Print*, "The exponents"
      Call output(icore(ialpha), 1, nfct, 1, 1, 1, nfct, 1, 1, 1)
      Call output(icore(ipcoeff), 1, nfct, 1, nbasp, nfct, nbasp, 1)
      Write(*,*)
      print*, "Exit Prep4_density_plots"
      Write(*,*)
#endif
C     
      RETURN
      END
