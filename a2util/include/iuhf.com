
#ifndef _IUHF_COM_
#define _IUHF_COM_

      integer           iuhf
      common /iuhf_com/ iuhf
      save   /iuhf_com/

#endif /* _IUHF_COM_ */

