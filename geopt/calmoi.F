C---------------------------------------------------------------------
C
      Subroutine calMOI(Natom,Coord,AMass,MOI)
C
C calculates the moments of inertia matrix for general molecules
C
C Coord : input coordinate (center of mass as the origin)
C
C---------------------------------------------------------------------
C
      Implicit double precision (a-h,o-z)
C
      Double Precision MOI(3,3),MOImn,MOImm
C
      Dimension Coord(3,Natom)
      Dimension AMass(Natom)
      Dimension COM(3)
C
      Logical debug
      Parameter (debug=.False.)
      If (debug) write(*,*) '\n--- In the calMOI'
C
C---------------------------------------------------------------------
C
C     Check input
C---------------------------------------------------------------------
C
*     --- print out ---
      If (debug) then
      Write(*,9910) Natom
      Write(*,9920)
      Do i=1,Natom
         Write(*,9900) (Coord(j,i),j=1,3)
      End do
      Write(*,9930)
      Write(*,9905) (AMass(i),i=1,Natom)
      End if
*     -----------------
 9910 Format(4X,'--- ','Number of atoms =',I5)
 9920 Format(4X,'--- ','Input coordinate')
 9930 Format(4X,'--- ','Atomic mass')
C
C--------------------------------------------------------------------
C     Calculation of MOI
C---------------------------------------------------------------------
C
      Do m=1,3
         Do n=1,3
c           If (im .NE. m) then
            If (n .NE. m) then
               MOImn=0.0d0
               Do i=1,Natom
                  factor=AMass(i)
                  MOImn=MOImn-Coord(m,i)*Coord(n,i)*factor
               End do
               MOI(m,n)=MOImn
            End if
         End do
         MOImm=0.0d0
         Do im=1,3
            If (im .NE. m) then
                  Do i=1,Natom
                     factor=AMass(i)
                     MOImm=MOImm+factor*Coord(im,i)**2
                  End do
            End if
         End do
         MOI(m,m)=MOImm
      End do
C
*     --- print out ---
      If (debug) then
      Write(*,9710)
      Do i=1,3
         Write(*,9900) (MOI(i,j),j=1,3)
      End do
      End if
*     -----------------
 9710 Format(4X,'--- ','Moments of inertia')
C
C---------------------------------------------------------------------
C
      If (debug) write(*,'(A)') '--- Out of the calMOI'
      Return
C
 9900 format(3X,3f20.10)
 9905 format(7X,5F12.6)
C
      End
C
C---------------------------------------------------------------------
C
