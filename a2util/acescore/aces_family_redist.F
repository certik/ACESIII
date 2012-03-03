
c This routine redistributes the array sizes for an MOIO family of arrays
c based on two-particle symmetry processing. It DOES NOT perform the resort,
c and it requires the groups to be 'packed' in the sense that aces_list_touch
c was not told to start each group on a physical record boundary (besides the
c first one). Although it updates the auxiliary cache addresses, every
c array in the family MUST be in memory with no spaces between them (just
c like the disk storage requirement).
c    This is really a shameful routine. The main idea is that a diagram
c like <Ab|Ij> is to be resorted and stored as <AI|bj>. The programmer
c does the resort with sstgen, calls this routine, and then calls putlst
c to update the disk records. People should consider better algorithms...

c INPUT
c int IRP_DIAG : the total symmetry of the diagram (normally 1)
c int IFAMNDX  : the family index of the arrays to redistribute
c int DISTTYPE_BRA : the distribution type (cf. disttype.h) of the bra (rows)
c int DISTTYPE_KET : the distribution type (cf. disttype.h) of the ket (cols)
c logical BUPDATE : a behavior flag
c                   = T;        update /sympop/isytyp(*,IFAMNDX)
c                   = F; do not update ...
c                   Note: This flag is so important that it is ignored
c                         and assumed to be true.

#include "lists.h" /* for MOIO dimensions */

c#define _DEBUG_ACES_FAMILY_REDIST

      subroutine aces_family_redist(irp_diag,iFamNdx,
     &                              DistType_bra,DistType_ket,
     &                              bUpdate)
      implicit none

c ARGUMENTS
      integer irp_diag, iFamNdx, DistType_bra, DistType_ket
      logical bUpdate

c INTERNAL VARIABLES
      integer irp_bra, irp_ket
      integer nRows, nCols, nWords
      integer iNumNew, iNumMax
      integer iRec, iWrdNdx
      integer quiktmp(_MAX_IO_GRPS)
      integer iTmp

c COMMON BLOCKS
#include "lists.com"
#include "filspc.com" /* for iprcwd */
#include "sympop.com" /* for irpdpd and isytyp */
#include "syminf.com" /* for nirrep and dirprd */
#include "auxcache.com" /* for quikget */
#include "machsp.com" /* for iintfp */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iTmp = 1
      end if
c   o assert irp_diag is in [1,nirrep]
      if ((irp_diag.lt.1).or.(nirrep.lt.irp_diag)) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   irp_diag = ',irp_diag
         print *, '   nirrep   = ',nirrep
         iTmp = 1
      end if
c   o assert iFamNdx is properly bound
      if ((iFamNdx.lt.1).or.(_MAX_IO_FAMS.lt.iFamNdx)) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   iFamNdx = ',iFamNdx
         print *, '   max ndx = ',_MAX_IO_FAMS
         iTmp = 1
      end if
c   o assert distribution types are in [1,22]
      if ((DistType_bra.lt.1).or.(22.lt.DistType_bra).or.
     &    (DistType_ket.lt.1).or.(22.lt.DistType_ket)    ) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   DistType_bra = ',DistType_bra
         print *, '   DistType_ket = ',DistType_ket
         print *, '   max type     = ',22
         iTmp = 1
      end if
c   o assert bUpdate is .true. (even though it's ignored)
      if (.not.bUpdate) then
         print *, '@ACES_FAMILY_REDIST: Assertion failed.'
         print *, '   bUpdate = ',bUpdate
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o process the first row/group/irrep (peeled off the do loop)
      irp_bra = dirprd(1,irp_diag)
      nRows   = irpdpd(irp_bra,DistType_bra)
      nCols   = irpdpd(1,      DistType_ket)
      iNumNew = nRows*nCols
      iNumMax = moiomxsz(1,iFamNdx)
      moiosz(1,iFamNdx) = nRows
      moiods(1,iFamNdx) = nCols

c   o process 'the others'
      if (nirrep.ne.1) then
         nWords = nRows*nCols*iintfp
         quiktmp(1) = nWords
         iWrdNdx = moiowd(1,iFamNdx) + nWords
         iTmp    = (iWrdNdx-1)/iprcwd
         iRec    = moio(1,iFamNdx) + iTmp
         iWrdNdx = iWrdNdx         - iTmp*iprcwd
         do irp_ket = 2, nirrep
            irp_bra = dirprd(irp_ket,irp_diag)
            nRows   = irpdpd(irp_bra,DistType_bra)
            nCols   = irpdpd(irp_ket,DistType_ket)
            iTmp    = nRows*nCols
            iNumNew = iNumNew + iTmp
            nWords  = iTmp*iintfp
            iNumMax = iNumMax + moiomxsz(irp_ket,iFamNdx)
            moiosz(irp_ket,iFamNdx) = nRows
            moiods(irp_ket,iFamNdx) = nCols
            moio  (irp_ket,iFamNdx) = iRec
            moiowd(irp_ket,iFamNdx) = iWrdNdx
            quiktmp(irp_ket) = nWords
            iWrdNdx = iWrdNdx + nWords
            iTmp    = (iWrdNdx-1)/iprcwd
            iRec    = iRec    + iTmp
            iWrdNdx = iWrdNdx - iTmp*iprcwd
         end do
         if (quikget(1,iFamNdx).ne.0) then
            do irp_ket = 2, nirrep
               quikget(irp_ket,iFamNdx) =   quikget(irp_ket-1,iFamNdx)
     &                                    + quiktmp(irp_ket-1)
            end do
         end if
c     end if (nirrep.ne.1)
      end if

cYAU - I don't think we should modify moiomxsz. Those elements are set in
c      aces_list_touch which also moves the free-space pointers.

c   o make sure programmers are not creating new data
      if (iNumNew.gt.iNumMax) then
#ifdef _DEBUG_ACES_FAMILY_REDIST
         print *, '@ACES_FAMILY_REDIST: The number of elements has ',
     &            'increased!'
         print *, '                     The original NEWTYP2 was ',
     &            'gracious enough to remove this'
         print *, '                     trap, but since we do not ',
     &            'trust our own programmers,'
         print *, '                     we have to complain.'
         print *, '                     The main problem is that the ',
     &            'new data overwrites the'
         print *, '                     next list family or that the ',
     &            'next list_touch will'
         print *, '                     overwrite the tail of this ',
     &            'family.'
         print *, '   array family index  = ',iFamNdx
         print *, '   maximum family size = ',iNumMax
         print *, '   new family size     = ',iNumNew
#endif
      end if

c   o update the type lookup table
      isytyp(1,iFamNdx) = DistType_bra
      isytyp(2,iFamNdx) = DistType_ket

c   o mark MOIO as modified
      bIOMod = .true.

      return
c     end subroutine aces_family_redist
      end

