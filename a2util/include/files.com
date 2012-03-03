#ifndef _FILES_COM_
#define _FILES_COM_
c files.com : begin
      integer        luout, moints
      common /files/ luout, moints
c files.com : end
#endif /* _FILES_COM_ */
