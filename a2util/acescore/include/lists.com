#ifndef _LISTS_COM_
#define _LISTS_COM_
c lists.com : begin

c These common blocks contain global information about the arrays in storage.
c Elements prepended with "bw" are for storing the file metadata while working
c on multiple references.

#include "lists.h" /* for MOIO dimensions */
#include "bwcc.com" /* for BWCC maxref */

#ifndef NO_EXTERNAL
      external aces_bd_lists
#endif

c moio  (iGrp,iFam) : the physical record that contains the first element
c                     of the array (iGrp,iFam)
c moiowd(iGrp,iFam) : the integer-word index of the first element
c moiods(iGrp,iFam) : the number of columns in the array
c moiosz(iGrp,iFam) : the number of rows    in the array
c moiofl(iGrp,iFam) : the external file unit that contains the array

      integer          moio  (_MAX_IO_GRPS,_MAX_IO_FAMS),
     &                 moiowd(_MAX_IO_GRPS,_MAX_IO_FAMS),
     &                 moiosz(_MAX_IO_GRPS,_MAX_IO_FAMS),
     &                 moiods(_MAX_IO_GRPS,_MAX_IO_FAMS),
     &                 moiofl(_MAX_IO_GRPS,_MAX_IO_FAMS),
     &               bwmoio  (_MAX_IO_GRPS,_MAX_IO_FAMS,maxref),
     &               bwmoiowd(_MAX_IO_GRPS,_MAX_IO_FAMS,maxref),
     &               bwmoiosz(_MAX_IO_GRPS,_MAX_IO_FAMS,maxref),
     &               bwmoiods(_MAX_IO_GRPS,_MAX_IO_FAMS,maxref),
     &               bwmoiofl(_MAX_IO_GRPS,_MAX_IO_FAMS,maxref)
      common /lists/   moio,   moiowd,   moiosz,   moiods,   moiofl,
     &               bwmoio, bwmoiowd, bwmoiosz, bwmoiods, bwmoiofl
      save   /lists/

c moiomxsz(iGrp,iFam) : the original length of a one-dimensional array
c                       (This is shameful. Arrays should not be re-dimensioned
c                        at will during a job.)

      integer               moiomxsz(_MAX_IO_GRPS,_MAX_IO_FAMS),
     &                    bwmoiomxsz(_MAX_IO_GRPS,_MAX_IO_FAMS,maxref)
      common /lists_mxsz/   moiomxsz,
     &                    bwmoiomxsz
      save   /lists_mxsz/

c pRec(i)    : the index of the physical record in file i containing free space
c              (i is the internal unit number of the storage file.)
c iIntOff(i) : the integer offset from the beginning of the physical record
c              needed to address the free space

      integer            pRec   (_MAX_IO_LUNS),
     &                   iIntOff(_MAX_IO_LUNS),
     &                 bwpRec   (_MAX_IO_LUNS,maxref),
     &                 bwiIntOff(_MAX_IO_LUNS,maxref)
      common /io_ptrs/   pRec,   iIntOff,
     &                 bwpRec, bwiIntOff
      save   /io_ptrs/

c bIOUp  : a flag for bombing in get/putlst if aces_io_init has not been called
c bIOMod : a flag for updating the records in aces_io_fin

      logical           bIOUp, bIOMod
      common /io_flags/ bIOUp, bIOMod
      save   /io_flags/

c lists.com : end
#endif /* _LISTS_COM_ */
