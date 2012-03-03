
c This routine reads the MOL file and initializes data in /MOL/.
c It should be called only by sb_com_mol(), but its main purpose is
c to recast icore() addresses into convenient variables.
c
c Determines totprim, maxshlprim, maxshlorb, and maxshell.  It also
c fills the arrays nprimatom, nshellatom, and naoatom.

      subroutine sb_basis(comppopv,compmemb,nangatom,nshellatom,
     &                    nprimatom,naoatom)
      implicit none

#include "maxangshell.par"
#include "kscore.com"
#include "mol.com"

      integer comppopv(*),compmemb(*),nangatom(*),nshellatom(*),
     &        nprimatom(*),naoatom(*)

      integer
     &    iorbit,nshell,ishell,nsubshell(maxangshell),isub,natomprim,
     &    natomao,natomang,ishellang,ishellprim,ishellao,isubprim,
     &    isubao,nline,icount,ieqatm,j,iatom,iunit
      character*80 line,molfil

      integer maxorbit
      parameter (maxorbit=1000)
      integer orbitpop(maxorbit)

      call callstack_push('SB_BASIS')

c See mol.com for an explanation of this.  It will disappear soon.
      if (memknown.ne.1) then
         call getrec(0,'JOBARC','COMPPOPV',j,orbitpop)
         if (j.gt.maxorbit) then
            print *, '@SB_BASIS: hard limit exceeded; increase maxorbit'
            print *, '           dim(COMPPOPV) = ',j
            call errex
         end if
         call getrec(1,'JOBARC','COMPPOPV',natomsx,orbitpop)
      end if

      iunit=3
      call gfname('MOL',molfil,j)
      open(unit=iunit,file=molfil(1:j),form='formatted',status='old')
      rewind(iunit)
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line

      icount=0
      totprim=0
      maxshlprim=0
      maxshlao=0
      maxshlorb=0
      maxshell=0
      maxangmom=0
      do iorbit=1,compnorb
         read(iunit,'(a25,12i5)') line,nshell,
     &                            (nsubshell(ishell),ishell=1,nshell)
         read(iunit,'(a)') line
         natomprim=0
         natomao=0
         natomang=0
         do ishell=1,nshell
c           ishellang= 1,3,6,10,...   (1s, 3p, 6d, 10f, ...)
            ishellang=rshift(ishell*(ishell+1),1)
            ishellprim=0
            ishellao=0
            do isub=1,nsubshell(ishell)
               read(iunit,'(2i5)') isubprim,isubao
               ishellprim=ishellprim+isubprim
               ishellao=ishellao+isubao

               nline=(isubao-3)/4
               if ((isubao-3).gt.(nline*4)) nline=nline+1
               nline=(nline+1)*isubprim
               do j=1,nline
                  read(iunit,'(a)') line
               end do
            end do
            natomao=natomao+ishellang*ishellao
            natomprim=natomprim+ishellang*ishellprim
            natomang=natomang+ishellang
            if (ishellprim.gt.maxshlprim) maxshlprim=ishellprim
            if (ishellao.gt.maxshlao) maxshlao=ishellao
            if ((ishellprim*ishellao).gt.maxshlorb)
     &         maxshlorb=ishellprim*ishellao
c        end do ishell=1,nshell
         end do
         if (natomang.gt.maxshell) maxshell=natomang
         if (nshell.gt.maxangmom) maxangmom=nshell
         if (memknown.ne.0) then
            do ieqatm=1,comppopv(iorbit)
               icount=icount+1
               iatom=compmemb(icount)
               nprimatom(iatom)=natomprim
               nshellatom(iatom)=natomang
               naoatom(iatom)=natomao
               nangatom(iatom)=nshell
            end do
         end if
         totprim=totprim+natomprim*orbitpop(iorbit)
c     end do iorbit=1,compnorb
      end do
      
      close(iunit)
      call callstack_pop
      return
      end

