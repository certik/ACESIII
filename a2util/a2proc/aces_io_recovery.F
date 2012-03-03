#include "flags.h"
C
      subroutine aces_io_recovery
      
      Implicit None
      Integer maxcor, length, iuhf, ioff
C
C Read list IO files to extract T, L or C amplitudes. 
C
#include "icore.com"
#include "machsp.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "istart.com"
c
      iuhf = 0
      if (iflags(h_IFLAGS_ref).ne.0) iuhf = 1
      maxcor=icrsiz
      ioff = 0
c
      call read_tlc_amps(icore(i0), maxcor, ioff, length, iuhf)
      call write_tlc_amps(icore(i0), length, "TGUESS  ")
      
      Return
      End 
