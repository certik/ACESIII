      Subroutine Stpsize(Omega, Deltaq, Delx, Nvibs)
C 
      Implicit Double Precision (A-H, O-Z)
C
#include "machsp.com"
C
      Double Precision Na
C
      Dimension Deltaq(Nvibs), Omega(Nvibs)
     
      Parameter(h=6.62608D-34, c=2.99792458D+10, Na=6.02214D+23)
      
      PI = Acos(-1.0D+00)

      Do Imode = 1, Nvibs
         Deltaq(Imode) = Delx/(2.0D+00*PI*SQRT(c*Omega(Imode)/
     &                   h/Na*1.0D-23))
      Enddo 
C
      Return
      End

