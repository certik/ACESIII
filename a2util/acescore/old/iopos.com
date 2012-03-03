#ifndef _IOPOS_COM_
#define _IOPOS_COM_
c iopos.com : begin
      integer        icrsiz, ichcsz, ioff(2), lenrec
      common /iopos/ icrsiz, ichcsz, ioff,    lenrec
c iopos.com : end
#endif /* _IOPOS_COM_ */
