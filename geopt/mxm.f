C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.

C    Routine on VAX to emulate CRAY SCILIB routine, MXM p. 4-22
C         LIBRARY manual, implemented in VAX single precision.
C
C This routine replaces the brain dead matrix multiply used here with
C   a call to the blas routine.
C SG 4/96

      SUBROUTINE MXM(A,NAR,B,NAC,C,NBC)
      INTEGER NAR,NAC,NBC
      DOUBLE PRECISION A(NAR,*),B(NAC,*),C(NAR,*)
      DOUBLE PRECISION ONE,ZILCH
      PARAMETER (ONE=1.0D0, ZILCH=0.0D0)
      CALL XGEMM('N','N',NAR,NBC,NAC,
     &           ONE,  A,NAR,
     &                 B,NAC,
     &           ZILCH,C,NAR)
      RETURN
      END

c NAME
c      MXM - Computes matrix-times-matrix product (unit increments)
c
c SYNOPSIS
c      CALL MXM(a,nra,b,nca,c,ncb)
c
c DESCRIPTION
c      MXM computes the nra-by-ncb matrix product C = AB of the nra-by-nca
c      matrix A and the nca-by-ncb matrix B.
c
c      This routine has the following arguments:
c
c      a      Real array of dimension (nra,nca).  (input)
c             Matrix A, the first factor.
c
c      nra    Integer.  (input)
c             Number of rows in A (same as number of rows in C).
c
c      b      Real array of dimension (nca,ncb).  (input)
c             Matrix B, the second factor.
c
c      nca    Integer.  (input)
c             Number of columns in A (same as number of rows in B).
c
c      c      Real array of dimension (nra,ncb).  (output)
c             Matrix C, the product AB.
c
c      ncb    Integer.  (input)
c             Number of columns in B (same as number of columns in C).
c
c NOTES
c      You should use the Level 3 Basic Linear Algebra Subprogram (Level 3
c      BLAS) SGEMM(3S) rather than MXM.  BLAS routines are preferred because
c      they are the de facto standard linear algebra interface.  Using Level
c      3 BLAS routines will enhance your program's portability, and also
c      should enhance its performance portability.
c
c      For example,
c
c           CALL MXM(A,NRA,B,NCA,C,NCB)
c
c      is equivalent to,
c
c           CALL SGEMM('N','N',NRA,NCB,NCA,1.0,A,NRA,B,NCA,0.0,C,NRA)
c
c      MXM is restricted to multiplying matrices that have elements stored by
c      columns in successive memory locations.  MXMA(3S) is a general
c      subroutine for multiplying matrices that can be used to multiply
c      matrices that do not satisfy the requirements of MXM (although SGEMM
c      also supersedes MXMA).  If B and C have only one column, MXV(3S) or
c      MXVA(3S) (both superseded by Level 2 BLAS routine SGEMV, see
c      SGEMV(3S)) are similar subroutines, each of which computes the product
c      of a matrix and a vector.
c
c CAUTIONS
c      The product must not overwrite either factor.  For example, the
c      following call will not (in general) assign the product AB to A:
c
c           CALL MXM(A,NRA,B,NCA,A,NCA)

