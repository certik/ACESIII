#ifndef _FILSPC_COM_
#define _FILSPC_COM_
c filspc.com : begin

c This common block contains the dimensions of the physical records used by the
c MOIO storage files.

c iprcln : the    byte-length of a physical record
c iprcwd : the integer-length of a physical record
c iprcrl : the    recl-length of a physical record

      integer         iprcln, iprcwd, iprcrl
      common /filspc/ iprcln, iprcwd, iprcrl
      save   /filspc/

c filspc.com : end
#endif /* _FILSPC_COM_ */
