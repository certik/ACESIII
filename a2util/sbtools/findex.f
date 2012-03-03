      function findex(str,substr)
      implicit none
      integer findex
      character*(*) str,substr
      integer strlen,slen,sublen
      slen=strlen(str)
      sublen=strlen(substr)
      findex=index(str(1:slen),substr(1:sublen))
      end
