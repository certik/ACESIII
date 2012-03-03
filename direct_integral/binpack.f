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
c This program implements a worst fit bin packing algorithm. It searches
c over the number of possible bins and the range of maximum bin sizes and
c minimizes the group mass variance.

c isort and dependencies come from http://www.netlib.org/slatec/src/isort.f

      subroutine binpack(nShells,mass,ShlMax, ishells)
      implicit none

      integer nShells, mass(*), ishells(*), ShlMax

      integer nCells, MasLim
      integer TotMass, lrg, sml

      integer GrpMax, nGrps, pGrp0, SmlGrp

      integer m(nshells), GrpNdx(nshells), GrpMass(nshells)
      integer FinGrps, FinNdx(nshells), FinMass(nshells)

      double precision avg, FinAvg, var, FinVar, dTmp
      logical bValid
      logical bAccept

      integer i, j

c ----------------------------------------------------------------------

      nCells = nShells
      if (nCells.lt.1) then
         print *,'Error in binpack: nshells = ',nshells
         call abort_job()
      endif


      TotMass = 0
      lrg = 0
      sml = 2**30
      do i = 1, nCells
         if (mass(i).lt.1) then
            print *,'Error in binpack: Cells must have a mass > 0.'
            call abort_job()
         endif
         TotMass = TotMass + mass(i)
         lrg = max(lrg,mass(i))
         sml = min(sml,mass(i))
         m(i) = mass(i)
         ishells(i) = i
      end do

      MasLim = ShlMax
      if (MasLim.lt.lrg) then
         print *,'Error in binpack: ',
     *           'The group limit is ',MasLim,
     *           ' which is less than the largest shell size (',
     *           lrg,').'
         call abort_job()
      endif

c   o sort the masses in descending order
      call isort(m,ishells,nCells,-2)

      do MasLim = lrg, ShlMax 
c   o initialize the final groups to the individual cells
      do i = 1, nCells
         FinNdx(i)  = i
         FinMass(i) = m(i)
      end do
      FinGrps = nCells
      call statx(FinMass,FinGrps,FinAvg,FinVar)

      do GrpMax = lrg, min(MasLim,TotMass)
         do nGrps = (TotMass+GrpMax-1)/GrpMax, nCells

c   o fill in the largest cells (pGrp0 points to the first group to scan from)
      pGrp0 = 1
      do i = 1, nGrps
         GrpNdx(i) = i
         GrpMass(i) = m(i)
         if (m(i).gt.GrpMax-sml) pGrp0 = pGrp0 + 1
      end do
      bValid=(pGrp0.le.nGrps.or.nGrps.eq.nCells)

c   o assign the remaining cells
      do i = nGrps+1, nCells
c      o find the smallest group
         SmlGrp = pGrp0
         do j = pGrp0+1, nGrps
            if (GrpMass(j).lt.GrpMass(SmlGrp)) SmlGrp = j
         end do
         GrpNdx(i) = SmlGrp
         GrpMass(SmlGrp) = GrpMass(SmlGrp) + m(i)
         if (GrpMass(SmlGrp).gt.GrpMax) bValid = .false.
      end do

      if (bValid) then

c      o find the variance of the group masses
         call statx(GrpMass,nGrps,avg,var)

         bAccept = .false.
         if (nGrps .eq. FinGrps .and.
     *       var .lt. FinVar) bAccept = .true.
         if (nGrps .lt. FinGrps) bAccept = .true.

c      o compare this grouping with the best known grouping
         if (bAccept) then
c            print *, '> lower var at ',nGrps,' groups under ',GrpMax
            do i = 1, nCells
               FinNdx(i)  = GrpNdx(i)
               FinMass(i) = GrpMass(i)
            end do
            FinAvg = avg
            FinVar = var
            FinGrps = nGrps
         end if
      end if

         end do
      end do
      enddo

 1000 continue

      call isort(FinNdx,ishells,nCells,2)
      do i = 1, nCells
         m(i) = mass(ishells(i))
      enddo

c   o print the original scheme
      call statx(m,nCells,avg,var)
      dTmp = (1.d0*lrg)/sml

c   o print the new results
      lrg = 0
      sml = 2**30
      do i = 1, FinGrps
         lrg = max(lrg,FinMass(i))
         sml = min(sml,FinMass(i))
      end do
      dTmp = (1.d0*lrg)/sml

      return
      end

c ----------------------------------------------------------------------

c m(1:n) is considered the population, not a sample.

      subroutine statx(m,n,avg,var)
      implicit none
      double precision avg, var, a, v, dTmp
      integer m(*), n, i
      if (n.eq.1) then
         avg = m(1)
         var = 0.d0
         return
      end if
      a = m(1)
      do i = 2, n
         dTmp = 1.d0/i
         a = a + (1.d0*m(i)-a)*dTmp
      end do
      dTmp = 1.d0*m(1)-a
      v = dTmp*dTmp
      do i = 2, n
         dTmp = 1.d0*m(i)-a
         v = v*(i-1)/i
         v = v + dTmp*dTmp/i
      end do
      var = v
      avg = a
      return
      end

      SUBROUTINE ISORT (IX, IY, N, KFLAG)
C***BEGIN PROLOGUE  ISORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2A
C***TYPE      INTEGER (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Kahaner, D. K., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   ISORT sorts array IX and optionally makes the same interchanges in
C   array IY.  The array IX may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      IX - integer array of values to be sorted
C      IY - integer array to be (optionally) carried along
C      N  - number of values in integer array IX to be sorted
C      KFLAG - control parameter
C            =  2  means sort IX in increasing order and carry IY along.
C            =  1  means sort IX in increasing order (ignoring IY)
C            = -1  means sort IX in decreasing order (ignoring IY)
C            = -2  means sort IX in decreasing order and carry IY along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***REVISION HISTORY  (YYMMDD)
C   761118  DATE WRITTEN
C   810801  Modified by David K. Kahaner.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to IX,IY. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  ISORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      INTEGER IX(*), IY(*)
C     .. Local Scalars ..
      REAL R
      INTEGER I, IJ, J, K, KK, L, M, NN, T, TT, TTY, TY
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  ISORT
      NN = N
      IF (NN .LT. 1) THEN
           print *,'ERROR in ISORT: ',
     *       ' The number of values to be sorted is not positive.'
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         print *,'ERROR IN ISORT:',
     +      'The sort control parameter, K, is not 2, 1, -1, or -2.'
         RETURN
      ENDIF
C
C     Alter array IX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            IX(I) = -IX(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort IX only
C
      M = 1
      I = 1
      J = NN
      R = 0.375E0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IX(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX(I) .GT. T) THEN
         IX(IJ) = IX(I)
         IX(I) = T
         T = IX(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (IX(J) .LT. T) THEN
         IX(IJ) = IX(J)
         IX(J) = T
         T = IX(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX(I) .GT. T) THEN
            IX(IJ) = IX(I)
            IX(I) = T
            T = IX(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (IX(L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (IX(K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX(L)
         IX(L) = IX(K)
         IX(K) = TT
         GO TO 40
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = IX(I+1)
      IF (IX(I) .LE. T) GO TO 80
      K = I
C
   90 IX(K+1) = IX(K)
      K = K-1
      IF (T .LT. IX(K)) GO TO 90
      IX(K+1) = T
      GO TO 80
C
C     Sort IX and carry IY along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IX(IJ)
      TY = IY(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX(I) .GT. T) THEN
         IX(IJ) = IX(I)
         IX(I) = T
         T = IX(IJ)
         IY(IJ) = IY(I)
         IY(I) = TY
         TY = IY(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (IX(J) .LT. T) THEN
         IX(IJ) = IX(J)
         IX(J) = T
         T = IX(IJ)
         IY(IJ) = IY(J)
         IY(J) = TY
         TY = IY(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX(I) .GT. T) THEN
            IX(IJ) = IX(I)
            IX(I) = T
            T = IX(IJ)
            IY(IJ) = IY(I)
            IY(I) = TY
            TY = IY(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (IX(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (IX(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX(L)
         IX(L) = IX(K)
         IX(K) = TT
         TTY = IY(L)
         IY(L) = IY(K)
         IY(K) = TTY
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = IX(I+1)
      TY = IY(I+1)
      IF (IX(I) .LE. T) GO TO 170
      K = I
C
  180 IX(K+1) = IX(K)
      IY(K+1) = IY(K)
      K = K-1
      IF (T .LT. IX(K)) GO TO 180
      IX(K+1) = T
      IY(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            IX(I) = -IX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END
