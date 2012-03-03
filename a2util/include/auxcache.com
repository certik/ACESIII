#ifndef _AUXCACHE_COM_
#define _AUXCACHE_COM_
c auxcache.com : begin

c The auxiliary cache is a programmer-controlled list cache. The dimensions are
c the same as those of the MOIO arrays. The programmer may load lists into icore
c memory and set quikget(?,?) to their icore addresses. When getlst and putlst
c operate on ANY list, quikget is checked to see if the list lives in icore.
c If so, the operation is performed on the in-core data instead of hitting the
c storage file(s).

c WARNING
c    There is no automatic updating of storage files with data in the auxiliary
c cache. If the memory-resident data is altered and must be stored to disk, then
c the quikget value must be destroyed (zeroed) and the actual location of the
c data must then be passed to putlst. The aces_auxcache_flush routine
c systematically stores the quikget values, calls putlst, and restores the
c quikget values.

#include "lists.h" /* for MOIO dimensions */

      integer           quikget(_MAX_IO_GRPS,_MAX_IO_FAMS)
      common /auxcache/ quikget
      save   /auxcache/

c auxcache.com : end
#endif /* _AUXCACHE_COM_ */
