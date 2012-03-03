#ifndef _OFFSETS_4DENS_
#define _OFFSET_4DENS_

      COMMON /OFFSETS_4DENS/ISCF_TDEN, ICOR_TDEN, ISCF_DDEN, 
     &                      ICOR_DDEN, IBEGIN_P_DENS
C
#endif /* _OFFSET_4DENS_ */


