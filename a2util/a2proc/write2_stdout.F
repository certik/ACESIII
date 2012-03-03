      Subroutine Write2_stdout(Etran, Cvtran, Stran, Erot, Cvrot,
     &                         Srot, Evib, Cvvib, Svib, Etot,
     &                         Cvtot, Stot, H) 
      
      Implicit Double Precision(A-H, O-Z)

      write(6,*)
      write(6,955)'E','Cv','S'
      write(6,955)'kJ/mol','J/mol-K','J/mol-K'
      write(6,960)'Electronic',0.0,0.0,0.0
      write(6,960)'Translational',etran,cvtran,stran   
      write(6,960)'Rotational',erot,cvrot,srot
      write(6,960)'Vibrational',evib,cvvib,svib
      write(6,960)'TOTAL',etot,cvtot,stot
      write(6,*)
      write(6,955)'kcal/mol','cal/mol-K','cal/mol-K'
      write(6,960)'Electronic',0.0,0.0,0.0
      write(6,960)'Translational',etran/4.184,cvtran/4.184,stran/4.184
      write(6,960)'Rotational',erot/4.184,cvrot/4.184,srot/4.184
      write(6,960)'Vibrational',evib/4.184,cvvib/4.184,svib/4.184
      write(6,960)'TOTAL',etot/4.184,cvtot/4.184,stot/4.184
 955  format(20x,3(5x,a9))
 960  format(a20,3(5x,f9.4))
      write(6,*)
      write(6,965)h,'kJ/mol'
      write(6,965)h/4.184,'kcal/mol'
 965  format(' H =',f9.4,a10)

      Return
      End
