#include "flags.h"
      subroutine a2rate_main
      implicit double precision (a-h,o-z)
C
#include "icore.com"
#include "machsp.com"
#include "flags.com"
#include "flags2.com"
#include "istart.com"
C-----------------------------------------------------------------------
      if (iflags(h_IFLAGS_ref).eq.0) iuhf = 0
      maxcor=icrsiz
C
      Call Getrec(20,'JOBARC','NREALATM',1,Nreals)
      Call Getrec(20,'JOBARC','ZMATATMS',1,Natoms)

      iocc   = i0
      iatmMass = iocc + 16
      icoord   = iatmMass + 3*Natoms*iintfp 
      ivgrad   = icoord   + 3*NAtoms*iintfp
      izgrad   = ivgrad   + 3*Nreals*iintfp
      ihess    = izgrad   + 3*Nreals*iintfp
      ifreq    = ihess    + 9*Natoms*Natoms*iintfp
      imap     = ifreq    + 3*Nreals*iintfp
      iord     = imap     + Nreals
      inext    = iord     + Nreals 

      if(inext-i0.gt.maxcor)call insmem('a2rate_main',inext-i0,maxcor)

      call a2rate_wrte(.true.,.true.,.true.,.true.,natoms,nreals,
     &                 icore(iocc),icore(iatmMass),icore(icoord),
     &                 icore(ivgrad),icore(izgrad),icore(ihess),
     &                 icore(ifreq),icore(imap),icore(iord),iuhf)
c ----------------------------------------------------------------------

      return
      end

