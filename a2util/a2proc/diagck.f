
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine DiagCk(M, N, Tol, A, LDA, Info)
      Integer M, N, LDA, Info
      Double precision A(LDA, N), Tol
      Integer I, J
      Info = 0
      Do 1000 J = 1, N
         Do 1010 I = 1, M
            If (Abs( A(i,j) ) .gt. Tol .AND. i .ne. j) Info = Info + 1
 1010    Continue
 1000 Continue
      Return
      End
