
c The following routines print out a matrix.  Because of the many ways that
c matrices are stored, this is a complicated set of tasks.
c
c Variables are:
c   mat    : an nrow x ncol rectangular matrix OR and nrow x nrow square
c            matrix
c   title  : a string containing the title
c   collab : a (short) label for columns
c   rowlab : a (short) label for rows
c   numwid : the width of each number when printed (W in fW.D)
c   ndec   : the number of decimal places to print (D in fW.D)
c   r0,r1  : the first/last row to print
c   c0,c1  : the first/last column to print
c   which  : 'U' or 'L' for upper/lower diagonal
c   diag   : 'D' or 'N' for inclusion or not of the diagonal

c The simplest sort of matrix is a rectangular matrix with all numbers
c stored and printing a rectangular part of the matrix.  The following
c two routines print all or part of a matrix.

      subroutine mat_print(mat,nrow,ncol,title,collab,rowlab,
     &                     numwid,ndec)
      implicit none
      integer nrow,ncol,numwid,ndec
      double precision mat(nrow,ncol)
      character*(*) title,collab,rowlab
      call mat_print_sub(mat,nrow,ncol,title,collab,rowlab,
     &                   numwid,ndec,1,nrow,1,ncol)
      return
      end

