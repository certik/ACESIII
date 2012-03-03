
c This routine takes a matrix or a block of a matrix and transforms it from
c one form to another.  This is for a matrix of the same general form as the
c fock matrix (i.e. a block diagonal matrix where each block is symmetric,
c there is one block for each irrep, and the size of each block is kept in
c the variables in the sym.com include file).
c
c Matrices may be stored in one of three forms:
c   triangular  : the matrix is stored in blocks of upper triangular matrices
c   square      : each block of the matrix is stored fully (but none of the
c                 zero's outside each block)
c   full        : the full matrix (including zero's) is stored

      subroutine mat_trans(intype,outtype,in,out,irr)

#include <aces.h>
      implicit none
#include <symm2_ks.com>
#include <callstack.com>
#include <mol.com>

      integer irr,intype,outtype
      M_REAL in(1),out(1)

      integer irrep
      if (intype.eq.outtype) return

      callstack_prev=callstack_curr
      callstack_curr='MAT_TRANS'

      irrep=irr

c Full matrix transformations

      if (intype.eq.MAT_SQR) then
        if (outtype.eq.MAT_TRI) then
          do irrep=1,nirrep
            call mat_trans_sqr_tri(in(irrsqroff(irrep)),
     &          out(irrtrioff(irrep)),numbasir(irrep),
     &          irrtrilen(irrep),1)
          end do
        else if (outtype.eq.MAT_FULL) then
          call dzero(out,nbastot*nbastot)
          do irrep=1,nirrep
            call mat_trans_sqr_sqr(in(irrsqroff(irrep)),out,
     &          numbasir(irrep),nbastot,1,irrorboff(irrep))
          end do
        else
          goto 900
        endif

      else if (intype.eq.MAT_TRI) then
        if (outtype.eq.MAT_SQR) then
          do irrep=1,nirrep
            call mat_trans_tri_sqr(in(irrtrioff(irrep)),
     &          out(irrsqroff(irrep)),irrtrilen(irrep),
     &          numbasir(irrep),1)
          end do
        else if (outtype.eq.MAT_FULL) then
          call dzero(out,nbastot*nbastot)
          do irrep=1,nirrep
            call mat_trans_tri_sqr(in(irrtrioff(irrep)),out,
     &          irrtrilen(irrep),nbastot,irrorboff(irrep))
          end do
        else
          goto 900
        endif

      else if (intype.eq.MAT_FULL) then
        if (outtype.eq.MAT_SQR) then
          do irrep=1,nirrep
            call mat_trans_sqr_sqr(in,out(irrsqroff(irrep)),
     &          nbastot,numbasir(irrep),irrorboff(irrep),1)
          end do
        else if (outtype.eq.MAT_TRI) then
          do irrep=1,nirrep
            call mat_trans_sqr_tri(in,out(irrtrioff(irrep)),nbastot,
     &          irrtrilen(irrep),irrorboff(irrep))
          end do
        else
          goto 900
        endif

c Block transformations

      else if (intype.eq.MAT_FULLSEG) then
        if (outtype.eq.MAT_TRISEG) then
          call mat_trans_sqr_tri(in,out(irrtrioff(irrep)),nbastot,
     &        irrtrilen(irrep),irrorboff(irrep))
        else if (outtype.eq.MAT_TRIBLK) then
          call dzero(out,maxirrtri)
          call mat_trans_sqr_tri(in,out,nbastot,irrtrilen(irrep),
     &        irrorboff(irrep))
        else if (outtype.eq.MAT_SQRSEG) then
          call mat_trans_sqr_sqr(in,out(irrsqroff(irrep)),
     &        nbastot,numbasir(irrep),irrorboff(irrep),1)
        else if (outtype.eq.MAT_SQRBLK) then
          call dzero(out,maxirrsqr)
          call mat_trans_sqr_sqr(in,out,nbastot,numbasir(irrep),
     &        irrorboff(irrep),1)
        else
          goto 900
        endif

      else if (intype.eq.MAT_TRISEG) then
        if (outtype.eq.MAT_FULLSEG) then
          call mat_trans_tri_sqr(in(irrtrioff(irrep)),out,
     &        irrtrilen(irrep),nbastot,irrorboff(irrep))
        else if (outtype.eq.MAT_TRIBLK) then
          call dzero(out,maxirrtri)
          call B_COPY(irrtrilen(irrep),in(irrtrioff(irrep)),1,out,1)
        else if (outtype.eq.MAT_SQRSEG) then
          call mat_trans_tri_sqr(in(irrtrioff(irrep)),
     &        out(irrsqroff(irrep)),irrtrilen(irrep),
     &        numbasir(irrep),1)
        else if (outtype.eq.MAT_SQRBLK) then
          call dzero(out,maxirrsqr)
          call mat_trans_tri_sqr(in(irrtrioff(irrep)),out,
     &        irrtrilen(irrep),numbasir(irrep),1)
        else
          goto 900
        endif

      else if (intype.eq.MAT_TRIBLK) then
        if (outtype.eq.MAT_FULLSEG) then
          call mat_trans_tri_sqr(in,out,irrtrilen(irrep),nbastot,
     &        irrorboff(irrep))
        else if (outtype.eq.MAT_TRISEG) then
          call B_COPY(irrtrilen(irrep),in,1,out(irrtrioff(irrep)),1)
        else if (outtype.eq.MAT_SQRSEG) then
          call mat_trans_tri_sqr(in,out(irrsqroff(irrep)),
     &        irrtrilen(irrep),numbasir(irrep),1)
        else if (outtype.eq.MAT_SQRBLK) then
          call dzero(out,maxirrsqr)
          call mat_trans_tri_sqr(in,out,irrtrilen(irrep),
     &        numbasir(irrep),1)
        else
          goto 900
        endif

      else if (intype.eq.MAT_SQRSEG) then
        if (outtype.eq.MAT_FULLSEG) then
          call mat_trans_sqr_sqr(in(irrsqroff(irrep)),out,
     &        numbasir(irrep),nbastot,1,irrorboff(irrep))
        else if (outtype.eq.MAT_TRISEG) then
          call mat_trans_sqr_tri(in(irrsqroff(irrep)),
     &        out(irrtrioff(irrep)),numbasir(irrep),
     &        irrtrilen(irrep),1)
        else if (outtype.eq.MAT_TRIBLK) then
          call dzero(out,maxirrtri)
          call mat_trans_sqr_tri(in(irrsqroff(irrep)),out,
     &        numbasir(irrep),irrtrilen(irrep),1)
        else if (outtype.eq.MAT_SQRBLK) then
          call dzero(out,maxirrsqr)
          call B_COPY(irrsqrlen(irrep),in(irrsqroff(irrep)),1,out,1)
        else
          goto 900
        endif

      else if (intype.eq.MAT_SQRBLK) then
        if (outtype.eq.MAT_FULLSEG) then
          call mat_trans_sqr_sqr(in,out,numbasir(irrep),nbastot,
     &        1,irrorboff(irrep))
        else if (outtype.eq.MAT_TRISEG) then
          call mat_trans_sqr_tri(in,out(irrtrioff(irrep)),
     &        numbasir(irrep),irrtrilen(irrep),1)
        else if (outtype.eq.MAT_TRIBLK) then
          call dzero(out,maxirrtri)
          call mat_trans_sqr_tri(in,out,numbasir(irrep),
     &        irrtrilen(irrep),1)
        else if (outtype.eq.MAT_SQRSEG) then
          call B_COPY(irrsqrlen(irrep),in,1,out(irrsqroff(irrep)),1)
        else
          goto 900
        endif

      else
        goto 900
      endif

      callstack_curr=callstack_prev
      return

  900 continue
      write(*,9000)
 9000 format('@MAT_TRANS-F: illegal transformation')
      stop
      end

