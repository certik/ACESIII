#ifndef _SYM_COM_
#define _SYM_COM_
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
#endif /* _SYM_COM_ */
