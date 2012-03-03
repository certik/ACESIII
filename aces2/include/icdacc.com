#ifndef _ICDACC_COM_
#define _ICDACC_COM_
c icdacc.com : begin
c Nevin 8-30-95 added record length common to facilitate change from
c Bytes to Words for SGI and DecAlpha
      integer         idaccm
      common /icdacc/ idaccm
c icdacc.com : end
#endif /* _ICDACC_COM_ */
