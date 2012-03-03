#ifndef _INFO_COM_
#define _INFO_COM_
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
#endif /* _INFO_COM_ */
