      
      
      SUBROUTINE HEADER(HEAD, IN, LUPRI)
C
C Printing of Headers in a neat format.
C
C Arguments:
C  
C  Head  : String to be printed
C  IN    : Controls the indentation
C          When it is greater than zero, indentations is the 
C          the value of the variable IN plus 1. Otherwise
C          it will center the string in the output page.
C
C  LUPRI : The out-put unit number
C
      CHARACTER HEAD*(*)
C
      LENGTH = LEN(HEAD)
C
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = (72 - LENGTH)/2 + 1
      END IF
C
      WRITE (LUPRI, '(/,80A)') (' ',I=1,INDENT), HEAD
      WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (LUPRI, '()')
C
      RETURN
      END
