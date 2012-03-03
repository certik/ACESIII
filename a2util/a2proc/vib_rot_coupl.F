#include "flags.h"
      subroutine vib_rot_coupl
C
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
      Call Getrec(20, 'JOBARC', 'LINEAR', 1, Iamlinear)
C
      If (Iamlinear .EQ. 1) Then
         Nxm6 = 3*Nreals - 5
      Else
         Nxm6 = 3*Nreals - 6
      Endif
      Nx  = 3*Nreals
C
      iatmMass = i0
      icoord   = iatmMass + 3*Natoms*iintfp 
      iBmat    = icoord   + 3*NAtoms*iintfp
      iGmat    = iBmat    + Nx*NXm6*iintfp
      ihess    = iGmat    + Nxm6*Nxm6*iintfp
      ieig_vec = ihess    + Nx*Nx*iintfp
      ivib_vec = ieig_vec + Nx*Nx*iintfp
      ibtmp    = ivib_vec + Nx*Nxm6*iintfp
      itmp1    = ibtmp    + Nx*Nxm6*iintfp 
      itmp2    = itmp1    + Nx*Nx*iintfp 
      immatx   = Itmp2    + Nx*Nx*iintfp
      immaty   = Immatx   + Nx*Nx*iintfp
      immatz   = Immaty   + Nx*Nx*iintfp
      imap     = Immatz   + Nx*Nx*iintfp
      inext    = Imap     + Natoms
C
      if(inext-i0.gt.maxcor)call insmem('vib_rot_coupl',inext-i0,
     &                      maxcor)

      call coriolis_cc(natoms,nreals,Nx,Nxm6,icore(iatmMass),
     &                 icore(icoord),icore(ibmat),icore(igmat),
     &                 icore(ihess),icore(ieig_vec),icore(ivib_vec),
     &                 icore(ibtmp),icore(itmp1),icore(itmp2),
     &                 icore(immatx),icore(immaty),icore(immatz),
     &                 icore(imap),iuhf)
c ----------------------------------------------------------------------

      return
      end

