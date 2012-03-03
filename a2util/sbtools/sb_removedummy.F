
c This routine removes dummy atoms from the molecule.  It is necessary
c since JODA stores the molecular information with dummy atoms (it needs
c that information each time), but no other module should ever need
c this.  So, by default, all dummy atoms are removed.

      subroutine sb_removedummy(atommass,coord,atomchrg,
     &                          fullmemb,compmemb,natomsx)
      implicit none

      double precision atommass(*),coord(3,*)
      integer atomchrg(*),fullmemb(*),compmemb(*),natomsx

      integer iatom,iatm,ndum

      iatm=0
      ndum=0
      do iatom=1,natomsx
         if (atomchrg(iatom).eq.0) then
            ndum=ndum+1
         else
            iatm=iatm+1
            atommass(iatm)=atommass(iatom)
            coord(1,iatm)=coord(1,iatom)
            coord(2,iatm)=coord(2,iatom)
            coord(3,iatm)=coord(3,iatom)
            atomchrg(iatm)=atomchrg(iatom)
            fullmemb(iatm)=fullmemb(iatm)-ndum
            compmemb(iatm)=compmemb(iatm)-ndum
c            fullpopv(iatm)=fullpopv(iatom)
c            comppopv(iatm)=comppopv(iatom)
         end if
      end do

      return
      end

