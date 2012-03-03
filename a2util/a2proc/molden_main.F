
c This program writes a molden format file containing the coordinates,
c basis functions for each center, molecular orbitals, frequencies,
c frequency coordinates, and normal mode vectors.
c
c Basis functions are limited to gaussian which is a constraint from aces2.
c
c the following files need to be present:
c	JOBARC
c	JAINDX
c	ZMAT
c	GENBAS
c       for frequencies NORMCO
c 
c Ken Wilson March 1998

#include "flags.h"
      subroutine molden_main
      implicit double precision (a-h,o-z)
      parameter (mxangmom=7)
      parameter (mxcoef=30)
#include "icore.com"
#include "machsp.com"
#include "flags.com"
#include "flags2.com"
#include "istart.com"
      character*32 szFile
      logical bExist
      dimension iangmom(mxangmom)
c
      iuhf = 1
      if (iflags(h_IFLAGS_ref).eq.0) iuhf = 0
      maxcor=icrsiz
      iunit =10
      if (iunit.ne.6) then
         szFile = 'MOLDEN.INPUT'
         inquire(file=szFile,exist=bExist)
         if (bExist) call f_remove(szFile)
         open(unit=iunit,file=szFile,form='formatted',err=2)
 1       goto 3
 2       print *, '@MOLDEN_MAIN: could not create/open ',szFile
         call errex
 3       continue
      end if
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      iexx=0
c      iexx=1
c      call nl_init('tomolden',ierr,.true.)
c      if(ierr.eq.1)then
c         write(*,*)' *tomolden namelist not found in ZMAT file'
c         call errex
c      end if
c      call nl_str('EXX',0,iexx)
c      call nl_term
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c geometry
c
      write(iunit,'(a)')'[Molden Format]'
      write(iunit,'(a)')'[Atoms] Angs'

      call getrec(20,'JOBARC','NATOMS',1,natoms)

      iatchrg=i0
      icoord=iatchrg+natoms+mod(natoms,2)
      ifuchrg=icoord+3*natoms*iintfp
      ifucoord=ifuchrg+natoms+mod(natoms,2)
      inext=ifucoord+3*natoms*iintfp

      if(inext-i0.gt.maxcor)call insmem('rdgeom-f',inext-i0,maxcor)

      call A2rd_geom(natoms,icore(iatchrg),icore(icoord),
     &               icore(ifuchrg),icore(ifucoord),iunit,.true.)
c--------------------------------------------------------------------------
c basis set
c
      write(iunit,'(a)')'[GTO]'

      iexp=inext 
      icoef=iexp+mxangmom*mxcoef*iintfp 
      inext=icoef+mxangmom*mxcoef**2*iintfp 

      if(inext-i0.gt.maxcor)call insmem('rdbasis-f',inext-i0,maxcor)

      call A2rd_basis(natoms,icore(iatchrg),icore(iexp),
     &                icore(icoef),nshells,iangmom,iunit,.true.)
c--------------------------------------------------------------------------
c orbitals from exact calculations
c
      if((iflags(h_IFLAGS_vib).ne.3).and.
     &   (iflags2(h_IFLAGS2_grad_calc).ne.3)) then
   
        call getrec(0,'JOBARC','SCFENEG ',Ilength,Dtmp)
        if (ilength .gt. 0) then
           write(iunit,*)
           write(iunit,'(a)')'[MO]'
        else
           if (iunit.ne.6) then
              close(iunit)
              print *, '@MOLDEN_MAIN: successfully created ',szFile
           end if
            Return
        Endif
C
        call getrec(-1,'JOBARC','NAOBASFN',1,nao)
        call getrec(-1,'JOBARC','NBASTOT',1,nmo)

        iener=inext
        iocc=iener+nmo*iintfp
        iorb=iocc+nmo+mod(nmo,2)
        inext=iorb+nao*nmo*iintfp

        if(inext-i0.gt.maxcor)call insmem('rdorb-f',inext-i0,maxcor)

        call molden_rdorb(icore(iener),icore(iocc),icore(iorb),
     &    icore(inext),nao,nmo,maxcor-(inext-i0),iuhf,iexx,iunit)
      endif
c-------------------------------------------------------------------------
c vibrational frequencies
c
      if(iflags(h_IFLAGS_vib).ne.0) then
        write(iunit,*)  
        write(iunit,'(a)')'[FREQ]'

        call getrec(-1,'JOBARC','NREALATM',1,nreal_atoms)
        nvib=3*nreal_atoms
 
        ifreq=inext
        ifreqco=ifreq+nvib*iintfp
        inormmd=ifreqco+3*nreal_atoms*iintfp
        inext=inormmd+nvib*3*nreal_atoms*iintfp

        if(inext-i0.gt.maxcor)call insmem('rdvib-f',inext-i0,maxcor)
        call A2rd_vib(nreal_atoms,icore(iatchrg),icore(ifreq),
     &    icore(icoord),icore(inormmd),Nimag,nvib,iunit,.true.)

      endif

c ----------------------------------------------------------------------

      if (iunit.ne.6) then
         close(iunit)
         print *, '@MOLDEN_MAIN: successfully created ',szFile
      end if

      return
      end

