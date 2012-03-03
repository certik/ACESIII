      Subroutine reord_vecs(Vector, Scratch, Iorder, Length, Print)

      Implicit Double Precision (A-H, O-Z)
C
      Logical Print
 
      Dimension Vector(Length), Scratch(Length), Iorder(Length)
C
      Call Dcopy(Length, Vector, 1, Scratch, 1)
   

      Do I = 1, Length
      
         Vector(I) = Scratch(Iorder(i))

      Enddo

#ifdef _DEBUG_LVL0
       if (print) then
       write(6,*)
       Write(6,"(a)") "The orginal array"
       write(6,"(5(1x,e14.7))") (Scratch(i),i=1, Length)
       Write(6,"(a)") "The reordered array"
       write(6,"(5(1x,e14.7))") (vector(i),i=1, Length)
       endif 
#endif

      Return
      End
