#ifndef _ISTART_COM_
#define _ISTART_COM_
c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
#endif /* _ISTART_COM_ */
