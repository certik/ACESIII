#ifndef _FLAGS2_COM_
#define _FLAGS2_COM_
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
#endif /* _FLAGS2_COM_ */
