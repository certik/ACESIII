C
CJDW  6/ 6/95. Add JFS new improved version of PDSWAP.
C
      SUBROUTINE PDSWAP(EVEC,IANGBF,SCR,NAO, NBAS)
      Integer NAO, NBas
      Double precision EVEC(NAO,NBAS),SCR(NAO,NBAS)
      Integer IANGBF(NAO),junk(10000)
cAP - 500 is the maximum number of shells (not basis functions)
      dimension nshlang(500),nshlnum(500)
      Integer MaxAng, IOne, ICnt, IJump, i, j, ik, iq, INew
      PARAMETER (MAXANG = 7)
      Integer ISIZ(MAXANG)
C
      Integer LUIn, LuOut, LuErr
      Parameter (LuIn = 5, LuOut = 6, LuErr = 0)
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA IONE /1/
      DATA ISIZ /1,3,6,10,15,21,28/
C
C READ IN SHELL INFORMATION
C
      CALL GETREC(20,'JOBARC','FULSHLNM',1,NSHELL)
      if (nshell.gt.500) then
         print *, '@PDSWAP: Assertion failed.'
         print *, '         maximum number of shells = 500'
         print *, '         nshell = ',nshell
         call errex
      end if
      CALL GETREC(20,'JOBARC','FULSHLTP',NSHELL,NSHLANG)
      CALL GETREC(20,'JOBARC','FULSHLSZ',NSHELL,NSHLNUM)
      CALL ZERO(SCR,NAO*NBAS)
C
C LOOP OVER NUMBER OF SHELLS
C
      ISTART=1
      DO 10 ISHELL=1,NSHELL
       IANGMOM=NSHLANG(ISHELL)
       NFUNC  =NSHLNUM(ISHELL)
       ISIZE=ISIZ(IANGMOM+1)*NFUNC
       IF(IANGMOM.GE.1.AND.NFUNC.GT.1)THEN
C
C SWAPPING IS REQUIRED.  DO IT.
C
        DO 100 J=ISTART,ISTART+ISIZE-1
         INEW=MOD(J-ISTART,NFUNC)*ISIZ(IANGMOM+1)+((J-ISTART)/NFUNC)+
     &        ISTART
         CALL SCOPY(NBAS,EVEC(J,1),NAO,SCR(INEW,1),NAO)
100     CONTINUE
       ELSE
C
C NO SWAPPING NEEDED FOR THIS SHELL.  JUST COPY EVEC INTO SCR.
C
        IPOS=ISTART
        DO 200 I=1,ISIZE
         CALL SCOPY(NBAS,EVEC(IPOS+I-1,1),NAO,SCR(IPOS+I-1,1),NAO)
200     CONTINUE
C
       ENDIF
C
       ISTART=ISTART+ISIZE
C
10    CONTINUE
C
c YAU : old
c     CALL ICOPY(NAO*NBAS*IINTFP,SCR,1,EVEC,1)
c YAU : new
      CALL DCOPY(NAO*NBAS,SCR,1,EVEC,1)
c YAU : end
C
      RETURN
      END
