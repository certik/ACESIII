
#ifndef _LIFO_COM_
#define _LIFO_COM_

#define _MAX_STACKS 100
      integer           lifo(4,_MAX_STACKS), nStacks
      common /lifo_com/ lifo, nStacks
      save   /lifo_com/

#endif /* _LIFO_COM_ */

