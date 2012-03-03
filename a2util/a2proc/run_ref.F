      Subroutine Run_ref(Deltaq, Deltax, Norm_coords, NReals, 
     &                   Nvibs)
C 
      Implicit Double Precision (A-H, O-Z)
C
#include "machsp.com"
C
      Double Precision Norm_coords
C
      Dimension Deltaq(Nvibs), Norm_coords(3*Nreals, Nvibs), 
     &          Deltax(3*Nreals)
      
      Call Xgemm("N", "N", 3*Nreals, 1, Nvibs, 1.0D0, Norm_coords,
     &            3*Nreals, Deltaq, Nvibs, 1.0D0, Deltax, 3*Nreals)
   
      Call Putrec(20, "JOBARC", "COORD  ", Nreals*3*IINTFP, Deltax)
C
#ifdef _DEBUG_LVL0
                Write(6,*)
                Write(6, "(a)") "The refrence geometry"
                Write(6, "(3F10.5)") (Deltax(i),i=1,3*Nreals)
#endif          
C
      Call aces_ja_fin
C    
C Run the point and read the energy and the gradients
C
      Call Runit("runaces2b")
      Call aces_ja_init

      Return
      End

