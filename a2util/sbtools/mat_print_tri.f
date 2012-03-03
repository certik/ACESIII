
c The following routines allow you to print part of the upper or lower
c diagonal part of a matrix, with or without the diagonal.  The matrix
c is stored as the full square matrix.

      subroutine mat_print_tri(which,diag,mat,nrow,title,collab,rowlab,
     &                         numwid,ndec)
      implicit none
      integer nrow,numwid,ndec
      character*1 which,diag
      double precision mat(nrow,nrow)
      character*(*) title,collab,rowlab
      call mat_print_tri_sub(which,diag,mat,nrow,title,collab,rowlab,
     &                       numwid,ndec,1,nrow,1,nrow)
      return
      end

