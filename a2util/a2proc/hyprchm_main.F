
c This program writes a drives creation of HYPERCHEM  format file containing 
c the coordinates, frequencies, normal mode vectors and intensities.
c
c for frequencies NORMCO
c 
c Using Ken Wilson's March 1998, Marshall Cory & Ajith Perera, 10/04.

      subroutine hyprchm_main
      implicit double precision (a-h,o-z)
      logical tbohr,hfcrap,linear
C
#include "icore.com"
#include "machsp.com"
#include "flags.com"
#include "istart.com"
#include "flags.h"
C
      maxcor=icrsiz
C
      call getrec(20,'JOBARC','NATOMS',1,natoms)

      iatchrg  = i0
      icoord   = iatchrg  + natoms+mod(natoms,2)
      ifuchrg  = icoord   + 3*natoms*iintfp
      ifucoord = ifuchrg  + natoms+mod(natoms,2)
      ifreqco  = ifucoord + 3*natoms*iintfp
      ifreq    = ifreqco  + 3*natoms*iintfp
      inormmd  = ifreq    + 3*natoms*iintfp
      intnsty  = inormmd  + 9*natoms*natoms*iintfp
      iord     = intnsty  + 3*natoms*iintfp
      itag     = iord     + 3*natoms
      inext    = itag     + 10*natoms
C
      if(inext-i0.gt.maxcor)call insmem('hyprchm-main',inext-i0,
     &                                   maxcor)

      call hyprchm_rd_geomvib(natoms,icore(iatchrg),icore(icoord),
     &                        icore(ifuchrg),icore(ifucoord),
     &                        icore(ifreqco),icore(ifreq),
     &                        icore(inormmd),icore(intnsty),
     &                        icore(iord),hfcrap,linear)


      tbohr = .true.
      if (iflags(h_IFLAGS_vib) .ne.0) hfcrap = .true.

      if(natoms .gt. 1) call wrt_hyprchm( natoms, icore(icoord), 
     .                  icore(ifreq), icore(inormmd), icore(intnsty),
     .                  icore(iatchrg),tbohr,hfcrap,linear,icore(iord),
     .                  icore(itag))

      print *, '@HYPERCHEM: successfully created interface files'

      return
      end

