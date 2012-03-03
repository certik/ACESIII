#ifndef _SYMPOP_COM_
#define _SYMPOP_COM_
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
#endif /* _SYMPOP_COM_ */
