
#ifndef _VMOL_COM_
#define _VMOL_COM_

c ilnbuf : the length of a buffer to read in chunks of integral files
c (There is no practical reason why this would be anything but 600.)

      integer           ilnbuf
      common /vmol_com/ ilnbuf
      save   /vmol_com/

#ifndef NO_EXTERNAL
      external sb_bd_vmol
#endif

#endif /* _VMOL_COM_ */

