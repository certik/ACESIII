C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine ISctr ( NZ, X, INDX, Y)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: isctr.f,v 1.1.1.1 2003/04/02 19:21:47 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     isctr -- Scatter elements of a vector into another
C
C SYNOPSIS
      Integer NZ, INDX (NZ), Y(*), X(NZ)
C
C ARGUMENTS
C     NZ      Number of elements to be gathered (input)
C     X       Source vector (input)
C     Indx    Vector of indices of Y to put elements of X into (input)
C     Y       Source vector (input)
C
C DESCRIPTION
C     Performs the operation Y( Indx(i) ) = X(i) for i = 1, NZ.
C     Only elements of Y listed in Indx are accessed.
C
C NOTES
C     This is the integer analog to [DS]Sctr of the Dodson, Lewis, and
C     Grimes Sparse-BLAS1 proposal.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     LOCAL VARIABLES
C
      INTEGER             I
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
         Y(INDX(I)) = X(I)
   10 CONTINUE
C
      RETURN
      END
