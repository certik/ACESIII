
c This routine gets all of the standard information from the JOBARC file
c for the molecule.

      subroutine sb_com_mol
      implicit none

#include "aces.h"
#include "machsp.com"
#include "kscore.com"
#include "symm2_ks.com"
#include "mol.com"

      integer setptr,i
      character*8 szInt8
      integer    zccoeff,pprims,ppriml
      call callstack_push('SB_COM_MOL')

      call getrec(1,'JOBARC','NATOMS',1,natoms)
      natomsx=natoms
      if (memknown.ne.1) then
         call getrec( 1,'JOBARC','NREALATM',1,         nrealatm)
         call getrec( 1,'JOBARC','NUCREP  ',iintfp,    nucrep)
         call getrec(-1,'JOBARC','LINEAR  ',1,         linear)
         call getrec(-1,'JOBARC','ORIENTMT',3*3*iintfp,orientmt)
         call getrec(-1,'JOBARC','NMPROTON',1,         nmproton)
         call getrec(-1,'JOBARC','NOCCORB ',2,         noccorb)
         call getrec(-1,'JOBARC','NVRTORB ',2,         nvrtorb)
      end if

      patomchrg=setptr(1,F_INTEGER,natoms)
#ifdef _DEBUG
      Print*, "Pointer for atom charge:", patomchrg
#endif  
      zatommass=setptr(1,F_REAL,   natoms)
#ifdef _DEBUG
      Print*,  "Pointer for atom mass:", zatommass
#endif 
 
      zcoord   =setptr(1,F_REAL,   3*natoms)

      if (memknown.ne.0) then
         call getrec(1,'JOBARC','ATOMCHRG',natoms,kscore(patomchrg))
         call getrec(1,'JOBARC','ATOMMASS',natoms*iintfp,
     &               dcore(zatommass))
         call getrec(1,'JOBARC','COORD',3*iintfp*natoms,dcore(zcoord))
      end if

      if (memknown.ne.1) then
c        Reading strings!!!
         call getrec(1,'JOBARC','COMPPTGP',iintfp,szInt8)
         compptgp(1:4)=szInt8(1:4)
         call getrec(1,'JOBARC','FULLPTGP',iintfp,szInt8)
         fullptgp(1:4)=szInt8(1:4)
         c1symmet=0
         if (compptgp.eq.'C1  ') c1symmet=1

         compnirr=nirrep
         call getrec( 1,'JOBARC','FULLNIRR',1,fullnirr)
         call getrec( 1,'JOBARC','COMPNORB',1,compnorb)
         call getrec( 1,'JOBARC','FULLNORB',1,fullnorb)
         call getrec(-1,'JOBARC','COMPORDR',1,compordr)
         call getrec(-1,'JOBARC','FULLORDR',1,fullordr)

         spinc(1)='A'
         spinc(2)='B'
      end if

      pfullclss=setptr(1,F_INTEGER,fullordr)
      pcompclss=setptr(1,F_INTEGER,compordr)
      pfullpopv=setptr(1,F_INTEGER,natoms)
      pcomppopv=setptr(1,F_INTEGER,natoms)
      pfullmemb=setptr(1,F_INTEGER,natoms)
      pcompmemb=setptr(1,F_INTEGER,natoms)
      if (memknown.ne.0) then
         call getrec(-1,'JOBARC','FULLCLSS',fullordr,kscore(pfullclss))
         call getrec(-1,'JOBARC','COMPCLSS',compordr,kscore(pcompclss))
         call getrec( 1,'JOBARC','FULLPOPV',natoms,  kscore(pfullpopv))
         call getrec( 1,'JOBARC','COMPPOPV',natoms,  kscore(pcomppopv))
         call getrec( 1,'JOBARC','FULLMEMB',natoms,  kscore(pfullmemb))
         call getrec( 1,'JOBARC','COMPMEMB',natoms,  kscore(pcompmemb))
      end if

c   o remove dummy atoms
      if (memknown.ne.0) then
         call sb_removedummy(dcore(zatommass),dcore(zcoord),
     &                       kscore(patomchrg),kscore(pfullmemb),
     &                       kscore(pcompmemb),natomsx)
      end if
      natoms=nrealatm

c ----------------------------------------------------------------------
c sb_getbas()

      if (memknown.ne.1) then
         call getrec(1,'JOBARC','NAOBASFN',1,naobasfn)
         call getrec(1,'JOBARC','NBASTOT ',1,nbastot)
      end if

      pnprimatom =setptr(1,F_INTEGER,natoms)
      pnshellatom=setptr(1,F_INTEGER,natoms)
      pnangatom  =setptr(1,F_INTEGER,natoms)
      pnaoatom   =setptr(1,F_INTEGER,natoms)


      call sb_basis(kscore(pcomppopv),kscore(pcompmemb),
     &              kscore(pnangatom),kscore(pnshellatom),
     &              kscore(pnprimatom),kscore(pnaoatom))

      pnshellprim=setptr(1,F_INTEGER,maxshell*natoms)
      pnshellao  =setptr(1,F_INTEGER,maxshell*natoms)
      pprimoff   =setptr(1,F_INTEGER,maxshell*natoms)
      paooff     =setptr(1,F_INTEGER,maxshell*natoms)

      zalpha     =setptr(1,F_REAL,totprim)
      zprimcoef  =setptr(1,F_REAL,totprim*naobasfn)

      zccoeff   =setptr(1,F_REAL,totprim*naobasfn)
      pprims   =setptr(1,F_INTEGER,naobasfn)
      ppriml   =setptr(1,F_INTEGER,naobasfn)

      if (memknown.ne.0) then
         call sb_prim(kscore(pcomppopv),kscore(pcompmemb),
     &                kscore(pnshellatom),
     &                kscore(pnshellprim),kscore(pnshellao),
     &                kscore(pprimoff),kscore(paooff),
     &                dcore(zalpha),dcore(zprimcoef),
     &     dcore(zccoeff),kscore(pprims),kscore(ppriml))
      end if

c ----------------------------------------------------------------------

      call callstack_pop
      return
      end

