      SUBROUTINE ATPT(PCOORD,NPTS,IPT,ACOORD,NATOM,ALPHA,
     $     PCOEFFN,MOMFCT,MAXANG,NTANGM,ITFCT,NBASP,PVAL,
     $     NBAS,PVAL2,IPCNT,APDX,APDY,
     $     APDZ,APDR,PCOEFSO,FINMO,FINSALC,IDORO,NORBOUTN)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION IXPWR(15),IYPWR(15),IZPWR(15)
      DIMENSION PCOORD(3*(NPTS+1))
      DIMENSION ACOORD(3*NATOM)
      DIMENSION APDX(NATOM),APDY(NATOM),APDZ(NATOM),
     $     APDR(NATOM)
      DIMENSION ALPHA(ITFCT)
      DIMENSION PCOEFFN(ITFCT*NBASP),PVAL(ITFCT*NBASP)
      DIMENSION PVAL2(ITFCT*NBASP)
      DIMENSION MOMFCT(NATOM*NTANGM)
      DIMENSION IPCNT(NBAS)
      DIMENSION PCOEFSO(ITFCT*NBASP)
      DIMENSION FINMO(NBASP),FINSALC(NBASP)
C     
C     Calculate the distance from this point to each atom
       CALL DIS(IPT,ACOORD,NATOM,PCOORD,NPTS,APDR,
     $     APDX,APDY,APDZ)
C     
CXXXXXXXXX
C     
C     Calculate the value of each primitive,for each
C     molecular orbital and sale, at each point
C     
CXXXXXXXXX
C
       IF (IDORO.NE.0) THEN
          DO 5 J=0,(NORBOUTN-1)
             K=1
             N=0
             FINMO(J+1)=0.0
             DO 10 I=1,NATOM
C     
                DO 20 NUMSHL=1,MAXANG
C     
                   M=0
                   DO 30 IIII=1,NUMSHL
                      M=M+IIII
 30                CONTINUE
C     
                   CALL LAME(NUMSHL,IXPWR,IYPWR,IZPWR,M)
C     
                   N=N+M
C     
                   DO 50 IJ=1,M
                      TEMPVAL1=(APDX(I)**REAL(IXPWR(IJ)))*
     $                     (APDY(I)**REAL(IYPWR(IJ)))*
     $                     (APDZ(I)**REAL(IZPWR(IJ)))
                      DO 40 ML=1,MOMFCT(N-(M-1))
C     Calculate the value of each primitive,for each
C     molecular orbital, at each point
                         PVAL(J*ITFCT+K)=TEMPVAL1*PCOEFFN(J*ITFCT+K)*
     $                        DEXP(-ALPHA(K)*(APDR(I))**2)
                         FINMO(J+1)=FINMO(J+1)+PVAL(J*ITFCT+K)
C     IF (IPT.EQ.1) THEN
C     WRITE(*,*) "PVAL(",J*ITFCT+K,")=~,PVAL(J*ITFCT+K)
C     WRITE(*,*) "XPWR=",IXPWR(IJ),",","YPWR=",IYPWR(IJ),
C     & ",","ZPWR=",IZPWR(IJ)
C     END IF
                         K=K+1
 40                   CONTINUE
 50                CONTINUE
C     
 20             CONTINUE
 10          CONTINUE
 5        CONTINUE
C     
       END IF
C     
CXXXXXXX
C     
       IF (IDORO.NE.1) THEN
C     
          DO 105 J=0,(NBASP-1)
             K=1
             N=0
             FINSALC(J+1)=0.0
             DO 110 I=1,NATOM
C     
                DO 120 NUMSHL=1,MAXANG
C     
                   M=0
                   DO 130 IIII=1,NUMSHL
                      M=M+IIII
 130               CONTINUE
C     
                   CALL LAME(NUMSHL,IXPWR,IYPWR,IZPWR,M)
C     
                   N=N+M
C     
                   DO 150 IJ=1,M
                        TEMPVAL2=(APDX(I)**REAL(IXPWR(IJ)))*
     $                     (APDY(I)**REAL(IYPWR(IJ)))*
     $                     (APDZ(I)**REAL(IZPWR(IJ)))
                        DO 140 ML=1,MOMFCT(N-(M-1))
C     Calculate the value of each primitive,for each
C     SALC, at each point
                           PVAL2(J*ITFCT+K)=TEMPVAL2*PCOEFSO(J*ITFCT+K)*
     $                          DEXP(-ALPHA(K)*(APDR(I))**2)
                           FINSALC(J+1)=FINSALC(J+1)+PVAL2(J*ITFCT+K)
C     IF (IPT.EQ.1) THEN
C     WRITE(*,*) "PVAL2(",J*ITFCT+K,")=",PVAL2(J*ITFCT+K)
C     WRITE(*,*) "XPWR=",IXPWR(IJ),",","YPWR=",IYPWR(IJ),
C     $ ",","ZPWR=",IZPWR(IJ)
C     END IF
                           K=K+1
 140                    CONTINUE
 150                 CONTINUE
C     
 120              CONTINUE
 110           CONTINUE
 105        CONTINUE
         END IF
C     
C     
C     IF (IPT.EQ.1) THEN
C     WRITE(*,*) ''pval''
C     CALL PRINTMAT(PVAL,ITFCT,NBASP)
C     WRITE(*,*) "pval2~'
C     CALL PRINTMAT(PVAL2,ITFCT,NBASP)
C     END IF
C     
         RETURN
         END SUBROUTINE ATPT
         
         
         






