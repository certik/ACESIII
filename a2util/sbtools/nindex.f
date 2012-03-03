      function nindex(str,substr,ind)
      implicit none
      integer nindex
      character*(*) str,substr
      integer ind,strlen,slen,sublen,indmax,i,ncindex
      slen=strlen(str)
      sublen=strlen(substr)
      indmax=slen-sublen+1
      if (sublen.gt.slen .or. ind.ge.indmax .or. ind.lt.0) then
         nindex=0
         return
      end if
      i=ind
 10   continue
      nindex=ncindex(str,substr(1:1),i)
      if (nindex.eq.0) return
      if (nindex.gt.indmax) then
         nindex=0
         return
      end if
      if (str(nindex:nindex+sublen-1).eq.substr(1:sublen)) return
      i=nindex
      goto 10
      end
