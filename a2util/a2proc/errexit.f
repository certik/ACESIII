C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
      Subroutine ErrExit(PrgNam,iStat,Text)
C
C     -----
C     print error message given as Text and stop program.
C     the end of the message has to be insisted with '#'.
C
C     -----
C     called by 
C
C     -----
C     PrgNam : name of the program from which this is called
C     iStat  : whatever dummy integer
C     Text   : error message to be printed out
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
      Implicit double precision (a-h,o-z)
      Character*10 PrgNam
      Character*70 Text
C
      Logical debug
      debug=.False.
      If (debug) write(*,*) '\n--- In the ErrExit'
C
C---------------------------------------------------------------------
C     Pulling out the error message
C     ----------------------------------------------------------------
      np=0
      istop=0
      Do ip=1,70
         If (Text(ip:ip) .EQ. '#' .AND. istop .EQ. 0) then
            istop=1
         Else if (istop .EQ. 0) then
            np=np+1
         End if
      End do
C
C     ----------------------------------------------------------------
C     Print error messages and stop program
C     ----------------------------------------------------------------
      Write(6,900) Text(1:np)
      Write(6,910) PrgNam
      call aces_exit(1)
C
  900 Format(A)
  910 Format('Error occured in',2X,A6,'.\n',
     &       'Program terminated.')
C
C---------------------------------------------------------------------
C
      If (debug) write(*,'(A)') '--- Out of the ErrExit'
      Return
C
      End
