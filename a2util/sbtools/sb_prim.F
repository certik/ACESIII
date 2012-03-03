
c This routine reads the MOL file and initializes data in /MOL/.
c It should be called only by sb_com_mol(), but its main purpose is
c to recast icore() and dcore() addresses into convenient variables.
c
c Loads primitive coefficients of the basis.

      subroutine sb_prim(comppopv,compmemb,nshellatom,nshellprim,
     &                   nshellao,primoff,aooff,alpha,primcoef,
     &          dpmcoef,nprimao,nprimaol)
      implicit none

#include "maxangshell.par"
#include "mol.com"
#include <machsp.com>
      integer comppopv(*),compmemb(*),nshellatom(*),
     &        nshellprim(maxshell,*),nshellao(maxshell,natoms),
     &        primoff(maxshell,*),aooff(maxshell,*)
      double precision alpha(*),primcoef(totprim,*)
       integer iip,ikk,ill,iloo,ioo,ijj,
     &  nprimao(naobasfn),nprimaol(naobasfn),incre

      double precision nthres,
     &   dpmcoef(totprim,naobasfn)
      double precision pi,s,tmp,xnorm,ai,aj,dpow
      integer
     &    iorbit,ishell,ishelloff,ishellang,isub,iang,
     &    iequiv,nsubshell(maxangshell),nshell,iatom,iprimoff,iaooff,
     &    isubprim,isubao,i,nline,iprim0,iprim1,iprim,iao0,iao1,iao,
     &    jprim,jprim0,jao0,jprim1,jao1,jao,icount,iequivatom,iunit
      character*80 line,molfil

      call callstack_push('SB_PRIM')

      pi = acos(-1.d0)
      nthres=10.d0**(-12.d0)
 9000 format(a25,12i5)
 9010 format(2i5)
 9020 format(4f18.10)

      iunit=3
      call gfname('MOL',molfil,i)
      open(unit=iunit,file=molfil,form='formatted',status='old')

c Set nshellprim and nshellao
      rewind(iunit)
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line

c iorbit     the orbit number
c icount     the pointer to compmemb for the current atom of the current orbit
c iatom      the atom at that location in compmemb
c iequivatom the current equivalent atom

      icount=1
      do iorbit=1,compnorb

         iatom=compmemb(icount)
         icount=icount+1
         read(iunit,9000) line,nshell,
     &                    (nsubshell(ishell),ishell=1,nshell)

         do iip=1,nshell
         end do
         read(iunit,'(a)') line
         ishellang=0
         ishelloff=1
         do ishell=1,nshell
c ishelloff=1,2,5,11,21,...  (where the first shell of ishellang starts)
c ishellang=1,3,6,10,15,...  (1s, 3p, 6d, 10f, etc)
            ishelloff=ishelloff+ishellang
       
            ishellang=ishell*(ishell+1)/2
             
     
            do isub=1,nsubshell(ishell)
               read(iunit,9010) isubprim,isubao
                 

               do iang=ishelloff,ishelloff+ishellang-1
                  nshellprim(iang,iatom)=nshellprim(iang,iatom)+isubprim

                  nshellao(iang,iatom)=nshellao(iang,iatom)+isubao

               end do
               nline=(isubao-3)/4
                
               if ((isubao-3).gt.(nline*4)) nline=nline+1
               nline=(nline+1)*isubprim
               do i=1,nline
                  read(iunit,'(a)') line
               end do
            end do
         end do

c Copy primitive and ao info to equivalent atoms
         if (comppopv(iorbit).gt.1) then
            do iequiv=2,comppopv(iorbit)
               iequivatom=compmemb(icount)
               icount=icount+1
               do ishell=1,nshellatom(iatom)
                  nshellprim(ishell,iequivatom)=nshellprim(ishell,iatom)
                  nshellao(ishell,iequivatom)=nshellao(ishell,iatom)
               end do
            end do
         end if

c     end do iorbit=1,compnorb
      end do

c Set up primoff and aooff
      iprimoff=1
      iaooff=1
      do iatom=1,natoms
         do ishell=1,nshellatom(iatom)
            primoff(ishell,iatom)=iprimoff
            iprimoff=iprimoff+nshellprim(ishell,iatom)
            aooff(ishell,iatom)=iaooff
            iaooff=iaooff+nshellao(ishell,iatom)
         end do
      end do

c Read in the exponents and primitive coefficients (renormalizing them)
      rewind(iunit)
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line

      icount=1
      do iorbit=1,compnorb
         iatom=compmemb(icount)
         read(iunit,9000) line,nshell,
     &                    (nsubshell(ishell),ishell=1,nshell)
         read(iunit,'(a)') line
         ishellang=0
         ishelloff=1
         do ishell=1,nshell
            ishelloff=ishelloff+ishellang
            ishellang=ishell*(ishell+1)/2
c iprim0,iprim1 : the first and last primitives in the current set
c iao0,iao1     : the first and last AOs in the current set
            iprim0=primoff(ishelloff,iatom)
            iao0=aooff(ishelloff,iatom)
            do isub=1,nsubshell(ishell)
               read(iunit,9010) isubprim,isubao
               iprim1=iprim0+isubprim-1
               iao1=iao0+isubao-1
c Read all exponents and coefficients

               do iprim=iprim0,iprim1
                  read(iunit,9020) alpha(iprim),
     &                             (primcoef(iprim,iao),iao=iao0,iao1)
                 do iao=iao0,iao1
               end do
               end do
c Renormalize

               do iao=iao0,iao1
                  s=0.d0
                  dpow=dble(ishell)+0.5d0
                  do iprim=iprim0,iprim1
                     ai=alpha(iprim)
                     do jprim=1,iprim
                        aj=alpha(jprim)
                        tmp=primcoef(iprim,iao)
     &                     *primcoef(jprim,iao)
     &                     *(2.d0*sqrt(ai*aj)/(ai+aj))**dpow
                        s=s+tmp
                        if (iprim.ne.jprim) s=s+tmp
                     end do
                  end do
                  xnorm=((0.5d0/pi)**(3.d0*0.25d0))/sqrt(s)
                  dpow=0.5d0*dpow
                  do iprim=iprim0,iprim1
                     primcoef(iprim,iao)=primcoef(iprim,iao)
     &                                  *xnorm
     &                                  *(4.d0*alpha(iprim))**dpow
                  end do
               end do
               iao0=iao0+isubao
               iprim0=iprim0+isubprim
            end do
         end do

c Copy to equivalent blocks (given iorbit,ishell)
c This time we have to loop over ALL equivalent atoms since we'll also
c need to fill in the remaining shells of the first atom.
         do iequiv=1,comppopv(iorbit)
            iequivatom=compmemb(icount)
            icount=icount+1
            ishellang=0
            ishelloff=1
            do ishell=1,nshell
               ishelloff=ishelloff+ishellang
               ishellang=ishell*(ishell+1)/2
               iprim0=primoff(ishelloff,iatom)
               iao0=aooff(ishelloff,iatom)
               do iang=ishelloff,ishelloff+ishellang-1
                  jprim0=primoff(iang,iequivatom)
                  jprim1=jprim0+nshellprim(iang,iequivatom)-1
                  jao0=aooff(iang,iequivatom)
                  jao1=jao0+nshellao(iang,iequivatom)-1
                  do jprim=jprim0,jprim1
                     iprim=(jprim-jprim0)+iprim0
                     alpha(jprim)=alpha(iprim)
                     do jao=jao0,jao1
                        iao=(jao-jao0)+iao0
                        primcoef(jprim,jao)=primcoef(iprim,iao)
                     end do
c                 end do jprim=jprim0,jprim1
                  end do
c              end do iang=ishelloff,ishelloff+ishellang-1
               end do
c           end do ishell=1,nshell
            end do
c        end do iequiv=1,comppopv(iorbit)
         end do
c     end do iorbit=1,compnorb
      end do

      close(iunit)

c          iloo=1
c          nprimao(1)=1
c           nprimaol(1)=0 
        do ikk=1,naobasfn
           incre=0
        do ill=1,totprim
           if(abs(primcoef(ill,ikk)) .gt. nthres) then
               if(incre .eq.0) then
                   nprimao(ikk)=ill
                   nprimaol(ikk)=ill
                   incre=incre+1
                else
                   nprimaol(ikk)=ill
                end if
            end if
         end do
         end do
                  













c              if(ill .le. nprimaol(ikk-1)) then
                
c                 if(incre .eq.0) then
c                    nprimao(ikk)=ill
c                  end if
c                  incre=incre+1 
c                  iloo=ill
c               dpmcoef(iloo,ikk)=primcoef(ill,ikk)
c              iloo=iloo+1
c              else
c                dpmcoef(iloo,ikk)=primcoef(ill,ikk)
c                iloo=iloo+1
c              end if
c          end if
c        end do
c             nprimaol(ikk)=iloo-1
c             write(*,*) ikk,nprimaol(ikk) 
c             if(ikk .lt. naobasfn) then
c               nprimao(ikk+1)=iloo
c            end if
c        end do

c       call putrec(1,'JOBRAC','CCOEFFC',totprim*naobasfn*iintfp
c     &             ,dpmcoef)
     
       call putrec(1,'JOBARC','NPRIMST',naobasfn,nprimao)
       call putrec(1,'JOBARC','NPRIMLT',naobasfn,nprimaol) 
      
      call callstack_pop
      return
      end

