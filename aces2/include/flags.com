#ifndef _FLAGS_COM_
#define _FLAGS_COM_
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
#endif /* _FLAGS_COM_ */
