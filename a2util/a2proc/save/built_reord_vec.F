      Subroutine Built_Reord_Vec(Nreal_atoms, Nprim_4shl, Iangmom,
     &                           IShell_4atm, Itot_prim, Nshls,
     &                           Inatural_order)

      Implicit integer(A-Z)
C
#include <mxatms.par>
C
      Dimension Nprim_4shl(Nreal_atoms*NShls),
     &          Iangmom(Nreal_atoms*Nshls), IShell_4atm(Nreal_atoms),
     &          Inatural_order(Itot_prim), Ispdfgh_4atm(mxatms,7),
     &          Ifwd_4shell(mxatms,7),Ihigh(mxatms),
     &          Ibwd_4shell(mxatms,7),Isavei(mxatms),
     &          Isavep(mxatms), Inprims_4atm(mxatms,7),
     &          Ifwd_4prims(mxatms,7), Ibwd_4prims(mxatms,7),
     &          Iprims(0:6), Nprims(mxatms),Isaven(mxatms,mxatms),
     &          Ioff(8)
     
      Data (Iprims(I), I=0, 6) /1, 3, 6, 10, 15, 21, 36/
      Data (Ioff(I), I=1, 8) /0, 0, 2, 7, 19, 34, 55, 91/
C      
      Icount   = 0
      Iatm_off = 0
C
      Do Iatm = 1, Nreal_atoms
C
         Isfun  = 0
         Ipfun  = 0
         Idfun  = 0
         Iffun  = 0
         Ihfun  = 0
         Ilfun  = 0
C
C
         Isprm  = 0
         Ipprm  = 0
         Idprm  = 0
         Ifprm  = 0
         Ihprm  = 0
         Ilprm  = 0
         Iprmc  = 0
         
         Iprmc = Iatm_off
C
         Do Ishell = 1, IShell_4atm(Iatm)
C
            Iang  = Iangmom(Iatm_off+ishell)
            Iprmc = Iprmc + 1
C
            If (Iang .Eq. 1)  Then
                Isfun = Isfun + 1
                Isprm = Isprm + Nprim_4shl(Iprmc)
                Ispdfgh_4atm(Iatm, 1) = Isfun
                Inprims_4atm(Iatm, 1) = Isprm
                Ihigh(Iatm) = 1
            Elseif (Iang .Eq. 3) Then
                Ipfun = Ipfun + 1
                Ipprm = Ipprm + Nprim_4shl(Iprmc)*3
                Ispdfgh_4atm(Iatm, 2) = Ipfun
                Inprims_4atm(Iatm, 2) = Ipprm
                Ihigh(Iatm) = 2
            Elseif (Iang .Eq. 6) Then
                Idfun = Idfun + 1
                Idprm = Idprm + Nprim_4shl(Iprmc)*6
                Ispdfgh_4atm(Iatm, 3) = Idfun
                Inprims_4atm(Iatm, 3) = Idprm
                Ihigh(Iatm) = 3
            Elseif (Iang .Eq. 10) Then
                Iffun = Iffun + 1
                Ifprm = Ifprm + Nprim_4shl(Iprmc)*10
                Ispdfgh_4atm(Iatm, 4) = Iffun
                Inprims_4atm(Iatm, 4) = Ifprm
                Ihigh(Iatm) = 4
            Elseif (Iang .Eq. 15) Then
                Igfun = Igfun + 1
                Igprm = Igprm + Nprim_4shl(Iprmc)*15
                Ispdfgh_4atm(Iatm, 5) = Igfun
                Inprims_4atm(Iatm, 5) = Igprm
                Ihigh(Iatm) = 5
            Elseif (Iang .Eq. 21) Then
                Ihfun = Ihfun + 1
                Ihprm = Ihprm + Nprim_4shl(Iprmc)*21
                Ispdfgh_4atm(Iatm, 6) = Ihfun
                Inprims_4atm(Iatm, 6) = Ihprm
                Ihigh(Iatm) = 6 
            Elseif (Iang .Eq. 36) Then
                Ilfun = Ilfun + 1
                Ilprm = Ilprm + Nprim_4shl(Iprmc)*36
                Ispdfgh_4atm(Iatm, 7) = Ilfun
                Inprims_4atm(Iatm, 7) = Ilprm
                Ihigh(Iatm) = 7
            Elseif (Iang .Gt. 36) Then
                Write(6, "(a,a,a)") "The WFN file is only ",
     &                   "available for up to L angular ",
     &                   "momentum."
                Call Errex
            Endif
C
            Do Iang =1,  Iangmom(Iatm_off)
C
               Do Iprim = 1, Nprim_4shl(Ishell)

                  Icount = Icount + 1

                  Inatural_order(Icount) = Icount             
           
               Enddo
            Enddo
         Enddo
         Iatm_off = Iatm_off + IShell_4atm(Iatm)
      Enddo
 
      Nprims(1) = 0
      Do Iatm = 2, Nreal_atoms 
         Nprims(Iatm) = 0
         Do Ishell = 1, 7
            Nprims(Iatm) = Nprims(Iatm) + Inprims_4atm(iatm-1,Ishell)
         Enddo
      Enddo
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Natural order"
      Write(6,"(6(1x,i3))")(Inatural_order(i),i=1,Itot_prim)
      Write(6,*)
      Write(6,*) "The # of s,p,d..shells per atom"
      do i=1, Nreal_atoms
      Write(6,"(7(1x,i2))")(Ispdfgh_4atm(i,j), j=1,Ihigh(i))
      enddo 
      Write(6,*)
      Write(6,*) "The # of s,p,d.. primitives per atom"
      do i=1, Nreal_atoms
      Write(6,"(7(1x,i2))")(Inprims_4atm(i,j), j=1,Ihigh(i))
      enddo 
      Write(6,*)
      Write(6,"(6(1x,i2))")(Nprims(i), i=1, Nreal_atoms)
#endif
          
      Do Iatm = 1, Nreal_atoms
         Ifwd_4shell(Iatm, 1) = Ispdfgh_4atm(Iatm, 1) 
         Ifwd_4prims(Iatm, 1) = Inprims_4atm(Iatm, 1) + Nprims(Iatm)
         Do Ishell = 2, Ihigh(Iatm)
            Ifwd_4shell(Iatm,Ishell) = Ifwd_4shell(Iatm, Ishell-1) + 
     &                                 Ispdfgh_4atm(Iatm,Ishell) - 1
            Ifwd_4prims(Iatm,Ishell) = Ifwd_4prims(Iatm, Ishell-1) +
     &                                 Inprims_4atm(Iatm,Ishell) - 
     &                                 Iprims(Ishell-1) 
         Enddo                 
C
         Ibwd_4shell(Iatm,Ihigh(Iatm))=Ispdfgh_4atm(Iatm,Ihigh(Iatm)) 
         Ibwd_4prims(Iatm,Ihigh(Iatm))= Inprims_4atm(Iatm,Ihigh(Iatm))
         Do Jshell = Ihigh(Iatm)-1, 1, -1
            Ibwd_4shell(Iatm,Jshell) = Ibwd_4shell(Iatm,Jshell+1) +  
     &                                  + Ispdfgh_4atm(Iatm,Jshell)
            Ibwd_4prims(Iatm,Jshell) = Ibwd_4prims(Iatm,Jshell+1) +
     &                                  + Inprims_4atm(Iatm,Jshell)
         Enddo
C
CSSS         Ifwd_4shell(Iatm,Ihigh(Iatm))=Ifwd_4shell(Iatm,Ihigh(Iatm))+1
      Enddo
      

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,*) "The fwd and bwd offset for Ang. Moms. and Prims"
      do i=1, Nreal_atoms
      Write(6,"(7(1x,i2))")(Ifwd_4shell(i,j), j=1,Ihigh(i))
      Write(6,"(7(1x,i3))")(Ifwd_4prims(i,j), j=1,Ihigh(i))
      enddo
      Write(6,*)
      do i=1, Nreal_atoms
      Write(6,"(7(1x,i2))")(Ibwd_4shell(i,j), j=1,Ihigh(i))
      Write(6,"(7(1x,i2))")(Ibwd_4prims(i,j), j=1,Ihigh(i))
      enddo
#endif
C
      Iatm_off = 0
      Do Iatm = 1, Nreal_atoms
C
         Do Ishell = 1, Ihigh(Iatm)
C
            Isavei(Ishell)=Iangmom(Iatm_off+Ifwd_4shell(Iatm,Ishell))
            Isavep(Ishell)=Nprim_4shl(Iatm_off+Ifwd_4shell(Iatm,
     &                                Ishell))
            Write(6,*) Isavep(Ishell)
C
            Do Jshell = 1, Ibwd_4shell(Iatm, Ishell+1)-1
                Iangmom(Iatm_off+Ifwd_4shell(Iatm,Ishell)+Jshell-1)= 
     &          Iangmom(Iatm_off+Ifwd_4shell(Iatm,Ishell)+Jshell)
C
                Nprim_4shl(Iatm_off+Ifwd_4shell(Iatm,Ishell)+
     &                     Jshell-1) = Nprim_4shl(Iatm_off+Ifwd_4shell
     &                                        (Iatm,Ishell)+Jshell)
             Write(6,*) Nprim_4shl(Iatm_off+Ifwd_4shell(Iatm,Ishell)+
     &                  Jshell-1)

C
            Enddo
         Enddo
C
         Do Ishell = 1, Ihigh(Iatm)
            Iangmom(Iatm_off + Ifwd_4shell(Iatm,Ihigh(Iatm))+Ishell-1)
     &              = Isavei(Ishell)
C
            Nprim_4shl(Iatm_off + Ifwd_4shell(Iatm,Ihigh(Iatm))+
     &                 Ishell-1) = Isavep(Ishell)
C
         Enddo
         iatm_off = iatm_off + IShell_4atm(Iatm)
      Enddo
      
      Iatm_off = 0
      Do Iatm = 1, Nreal_atoms
C
         Do Ishell = 1, Ihigh(Iatm)

            Do Jshell = 1, Iprims(Ishell-1)
            Isaven(Ishell,Jshell)=Inatural_order(Iatm_off+
     &                            Ifwd_4prims(Iatm, Ishell)+Jshell-1)
            Enddo 
C
            Do Kshell = 1, Iprims(Ishell-1)
               Do Jshell = 1, Ibwd_4prims(Iatm, Ishell+1)+
     &                                    Iprims(Ishell-1)
                  Inatural_order(Iatm_off+Ifwd_4prims(Iatm,Ishell)+
     &                           Jshell-1) = Inatural_order(Iatm_off +
     &                                               Ifwd_4prims(Iatm,
     &                                               Ishell)+Jshell)
               Enddo
            Enddo
         Enddo
C
         Do Ishell = 1, Ihigh(Iatm)

            Ioffst = Ioff(Ishell)
            Do Jshell = 1, Iprims(Ishell-1)
               Inatural_order(Iatm_off+Ifwd_4prims(Iatm,Ihigh(Iatm))+
     &                        Ishell-1+Jshell-1+ioffst) = 
     &                        Isaven(Ishell,Jshell)
            Enddo 
        
         Enddo
      Enddo
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The reordered Ang. Moms. array"
      joff = 0
      Do I= 1, Nreal_atoms
      Write(6, "(4(1x,I3))") (Iangmom(joff+J), J=1,Ishell_4atm(I))
      joff = joff + IShell_4atm(I) 
      Enddo
      Write(6,*)
      Write(6,"(a)") "The reordered primitive array"
      joff = 0
      Do I= 1, Nreal_atoms
      Write(6, "(4(1x,I3))") (Nprim_4shl(joff+J), J=1,Ishell_4atm(I))
      joff = joff + IShell_4atm(I)
      Enddo 
      Write(6,*)
      Write(6,"(a)") "The reordered Natural order array"
      Write(6, "(10(1x,I3))") (Inatural_order(J), J=1,
     &                        Itot_prim)
#endif
C
      Return
      End 
