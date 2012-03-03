
c This routine initializes the chemical system in order that the executable
c understands symmetry, orbital, and list information.

      subroutine aces_init_chemsys
      implicit none

c COMMON BLOCKS
#include "syminf.com" /* for nirrep */

c ----------------------------------------------------------------------

c   o nocco and nvrto
      call aces_com_info

c   o initialize the point group and direct product table
      call aces_com_syminf
      if (nirrep.lt.1) then
         print '(/)'
         print *, '@ACES_INIT_CHEMSYS: Information regarding the ',
     &            'chemical system has not been'
         print *, '                    created. Common blocks relying ',
     &            'on this information will'
         print *, '                    not be initialized.'
         print '(/)'
         return
      end if

c   o irpdpd (needs /syminf/)
      call aces_com_sympop

c   o pop and vrt (needs /syminf/)
      call aces_com_sym

c   o isymoff (needs /syminf/ and /sym/)
      call aces_com_symloc

      return
c     end subroutine aces_init_chemsys
      end

