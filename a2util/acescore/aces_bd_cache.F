
c This routine defines the default values for the common blocks in
c cache.com.

      blockdata aces_bd_cache
      implicit none
#define NO_EXTERNAL
#include "cache.com"
      data cachemiss, cacheskip,
     &     cacheread, cachewrite, cachewriteback /5*0/
      data bCacheUp /.false./
      end

