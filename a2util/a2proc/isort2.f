C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
      Subroutine iSort2(Ndim,K,M,iKey)
C
C     -----
C     sorts the input array and return the result in the same array
C
C     -----
C     iKey=1: sort downward
C          2: sort upward
C
C     -----
C     called by
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
      Implicit double precision (a-h,o-z)
      Dimension K(*),M(*)
C
      Logical debug
      Parameter (debug=.False.)
      If (debug) write(*,*) '\n--- In the iSort2'
C
C---------------------------------------------------------------------
C
C     Check input
C     ----------------------------------------------------------------
C
*     --- print out ---
      If (debug) then
      Write(*,9920) Ndim
      Write(*,9930) (K(i),i=1,Ndim)
      Write(*,9930) (M(i),i=1,Ndim)
      Write(*,9940) iKey
      End if
*     -----------------
 9920 Format(4X,'--- ','Ndim =',I7)
 9930 format(4X,5I10)
 9940 Format(4X,'--- ','iKey =',I7)
C
C---------------------------------------------------------------------
C
C     Sort downward
C     ----------------------------------------------------------------
      If (iKey .EQ. 1) then
C
         Do i=1,Ndim
            kMax=K(i)
            jMax=i
            Do j=i+1,Ndim
               If (K(j) .GT. kMax) then
                  kMax=K(j)
                  jMax=j
               End if
            End do
            Mtmp=M(jMax)
            Do l=jMax,i+1,-1
               K(l)=K(l-1)
               M(l)=M(l-1)
            End do
            K(i)=kMax
            M(i)=Mtmp
         End do
C
*        --- print out ---
         If (debug) then
         Write(*,9810)
         Write(*,9930) (K(i),i=1,Ndim)
         Write(*,9930) (M(i),i=1,Ndim)
         End if
*        -----------------
 9810    Format(4X,'--- ','sort downward')
C
C
C     Sort upward
C     ----------------------------------------------------------------
      Else if (iKey .EQ. 2) then
C
         Do i=1,Ndim
            kMin=K(i)
            jMin=i
            Do j=i+1,Ndim
               If (K(j) .LT. kMin) then
                  kMin=K(j)
                  jMin=j
               End if
            End do
            Mtmp=M(jMin)
            Do l=jMin,i+1,-1
               K(l)=K(l-1)
               M(l)=M(l-1)
            End do
            K(i)=kMin
            M(i)=Mtmp
         End do
C
*        --- print out ---
         If (debug) then
         Write(*,9710)
         Write(*,9930) (K(i),i=1,Ndim)
         Write(*,9930) (M(i),i=1,Ndim)
         End if
*        -----------------
 9710    Format(4X,'--- ','sort upward')
C
C
C     neither
C     ----------------------------------------------------------------
      Else
         Call ErrExit('dSort',0,'Invalid iKey.#')
      End if
C
C---------------------------------------------------------------------
C
      If (debug) write(*,'(A)') '--- Out of the iSort2'
      Return
C
      End
