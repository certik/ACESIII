      subroutine Read_IRCnmlist(Stride, Trans_state, Direction, 
     &                          Max_IRC_Step)
C
      Implicit Double Precision (a-h,o-z)
C
      Logical Printdef, Trans_state, Forward
      Character*80  Direction 

c Call to nl_init find the *IRC namelist if it is in the ZMAT

      printdef = .True.
      call nl_init('IRC', ierr, printdef)

C Step size control the size of the individual steps.

      call nl_real('STEP_SIZE', 0.3D0, Stride)

c Starting point for the IRC search, it can be from the
C saddle point or any other point on the IRC. The saddle
C point (the default) is preferred. 

      call nl_log('SADDLE', .TRUE., Trans_state)

C The direction of the search, Forward=.true. towards products
C and false the opposite. 
 
      call nl_str('IRC_SEARCH', "FORWARD", Direction)

C The number of steps for the entire IRC search
  
      call nl_int("MAX_IRC_STEPS", 20, Max_IRC_Step)
C
      call nl_term

      RETURN
      END
