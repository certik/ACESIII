#ifndef _METHOD_COM_
#define _METHOD_COM_
c method.com : begin
      integer iuhf
      logical scf, nonhf
      common /method/ iuhf, scf, nonhf
c method.com : end
#endif /* _METHOD_COM_ */
