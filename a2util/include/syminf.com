#ifndef _SYMINF_COM_
#define _SYMINF_COM_
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
#endif /* _SYMINF_COM_ */
