
c This routine redimensions the array sizes for an MOIO family of arrays
c based on two arrays: one of the number of rows and the other of the number
c of columns. It requires the groups to be 'packed' in the sense that
c aces_list_touch was not told to start each group on a physical record
c boundary (besides the first one). Although it updates the auxiliary cache
c addresses, every array in the family MUST be in memory with no spaces between
c them (just like the disk storage requirement).
c    This is really a shameful routine. Only one routine uses it: vea(updlsts2).
c NO ONE should use this routine besides sloppy programmers. (The reason why
c it lives here is that only libacescore.a knows about lists.com.)

c INPUT
c int IO_ROWS : the number of arrays in the family
c int IFAMNDX : the family index of the arrays to redimension
c int ROWS_PER_IRP(IO_ROWS) : the one-dimensional array of rows by irrep
c int COLS_PER_IRP(IO_ROWS) : the one-dimensional array of cols by irrep
c int ISYTYP1 : the value (distribution type) to enter in to isytyp(1,IFAMNDX)
c int ISYTYP2 : the value (distribution type) to enter in to isytyp(2,IFAMNDX)

#include "lists.h" /* for MOIO dimensions */

      subroutine aces_family_danger(IO_Rows,     iFamNdx,
     &                              Rows_Per_Irp,Cols_Per_Irp,
     &                              iSyTyp1,     iSyTyp2)
      implicit none

c ARGUMENTS
      integer IO_Rows,         iFamNdx
      integer Rows_Per_Irp(*), Cols_Per_Irp(*)
      integer iSyTyp1,         iSyTyp2

c INTERNAL VARIABLES
      integer nRows, nCols, nWords
      integer irrep
      integer iRec, iWrdNdx
      integer quiktmp(_MAX_IO_GRPS)
      integer iTmp

c COMMON BLOCKS
#include "lists.com"
#include "iopos.com" /* for lenrec */
#include "sympop.com" /* for irpdpd and isytyp */
#include "syminf.com" /* for nirrep and dirprd */
#include "auxcache.com" /* for quikget */
#include "machsp.com" /* for iintfp */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_FAMILY_DANGER: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iTmp = 1
      end if
c   o assert IO_Rows is properly bound
      if ((IO_Rows.lt.1).or.(_MAX_IO_GRPS.lt.IO_Rows)) then
         print *, '@ACES_FAMILY_DANGER: Assertion failed.'
         print *, '   IO_Rows = ',IO_Rows
         print *, '   max ndx = ',_MAX_IO_GRPS
         iTmp = 1
      end if
c   o assert iFamNdx is properly bound
      if ((iFamNdx.lt.1).or.(_MAX_IO_FAMS.lt.iFamNdx)) then
         print *, '@ACES_FAMILY_DANGER: Assertion failed.'
         print *, '   iFamNdx = ',iFamNdx
         print *, '   max ndx = ',_MAX_IO_FAMS
         iTmp = 1
      end if
c   o assert dimensions are whole
      do irrep = 1, IO_Rows
         if ((Rows_Per_Irp(irrep).lt.0).or.
     &       (Cols_Per_Irp(irrep).lt.0)    ) then
            print *, '@ACES_FAMILY_DANGER: Assertion failed.'
            print *, '   Rows = ',Rows_Per_Irp(irrep)
            print *, '   Cols = ',Cols_Per_Irp(irrep)
            iTmp = 1
         end if
      end do
c   o assert distribution types are in [1,22]
      if ((iSyTyp1.lt.1).or.(22.lt.iSyTyp1).or.
     &    (iSyTyp2.lt.1).or.(22.lt.iSyTyp2)    ) then
         print *, '@ACES_FAMILY_DANGER: Assertion failed.'
         print *, '   iSyTyp1  = ',iSyTyp1
         print *, '   iSyTyp2  = ',iSyTyp2
         print *, '   max type = ',22
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

      if (IO_Rows.lt.1) return

c ----------------------------------------------------------------------

c   o process the first row/group/irrep (peeled off the do loop)
      nRows = Rows_Per_Irp(1)
      nCols = Cols_Per_Irp(1)
      moiosz(1,iFamNdx) = nRows
      moiods(1,iFamNdx) = nCols

c   o process 'the others'
      if (nirrep.ne.1) then
         nWords = nRows*nCols*iintfp
         quiktmp(1) = nWords
         iWrdNdx = moiowd(1,iFamNdx) + nWords
         iTmp    = (iWrdNdx-1)/lenrec
         iRec    = moio(1,iFamNdx) + iTmp
         iWrdNdx = iWrdNdx         - iTmp*lenrec
         do irrep = 2, nirrep
            nRows  = Rows_Per_Irp(irrep)
            nCols  = Cols_Per_Irp(irrep)
            nWords = nRows*nCols*iintfp
            moiosz(irrep,iFamNdx) = nRows
            moiods(irrep,iFamNdx) = nCols
            moio  (irrep,iFamNdx) = iRec
            moiowd(irrep,iFamNdx) = iWrdNdx
            quiktmp(irrep) = nWords
            iWrdNdx = iWrdNdx + nWords
            iTmp    = (iWrdNdx-1)/lenrec
            iRec    = iRec    + iTmp
            iWrdNdx = iWrdNdx - iTmp*lenrec
         end do
         if (quikget(1,iFamNdx).ne.0) then
            do irrep = 2, nirrep
               quikget(irrep,iFamNdx) =   quikget(irrep-1,iFamNdx)
     &                                  + quiktmp(irrep-1)
            end do
         end if
c     end if (nirrep.ne.1)
      end if

c   o update the type lookup table
      isytyp(1,iFamNdx) = iSyTyp1
      isytyp(2,iFamNdx) = iSyTyp2

c   o mark MOIO as modified
      bIOMod = .true.

      return
c     end subroutine aces_family_danger
      end

