      SUBROUTINE BLTAOTOSO(I2,ITRAN,CUN,CONTUN,CTEMP,NU,NR,AT,KABS,
     &                  I1,NXYZ,NUC,NRC,NBAS,CORR,KA,C,JCO,PT,
     &                  TITLE,JTRAN,A,CTRAN,ICDS,NDEG,CHARGT,X,
     &                  IAD,NFUNMX,NHT2,KWD,NH4,
     &                  JtranC, ITranC, CTranC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     (Other) Parameters
C-----------------------------------------------------------------------
C
#include "baslims.par"
C
      PARAMETER (NPRIMX = 10)
      PARAMETER (NCONTS = 15, NCONTP = 8, NCONTD = 4)
      PARAMETER (NCONTF = 2, NCONTG = 1)
      PARAMETER (NCONTH = 1, NCONTI = 1)
C-----------------------------------------------------------------------
      INTEGER  NAMN(MXCBF),LAMN(MXTNSH),MAMN(MXCBF)
      INTEGER  IAD(NHT2)
      INTEGER  IFBA(100)
      INTEGER  NCONTR(7)
C
      character*6 kaset
      character*6 slask
      character*8 hblnk
      character*8 nanum(8)
      character*1 ispd(7)
      character*1 ko(3)
      character*1 kka(8,3) 
      character*1 iblank
      character*2 koko(3)
      character*4 kwo(84)
      character*4 jprx(MXCBF)
      character*4 jjprx(MXCBF)
      character*4 kktyp(84,2)
      character*6 iblnk
cch
      INTEGER   AND,OR,EOR,DSTRT
      LOGICAL IFIF,LPRT,LHARM,DELE
      DIMENSION NXBAS(8),JCNTAM(28,7),KDELB(200)
      DIMENSION I2(21), ITRAN(MXCBF,224), CUN(MXTNCC),
     1            CONTUN(MXTNCC), CTEMP(NFUNMX,NFUNMX), 
     2            NU(NFUNMX), NR(NFUNMX), AT(NFUNMX),
     3            KABS(MXSHEL,28,15,8)
      dimension iadq(8),charm(1596)
C SG 6/25/97 Changed to agree with Stanton's code to use 500 functions
      dimension last(MXCBF),iast(MXCBF)
      dimension nast(MXCBF),jaast(MXCBF)
C
C     DIMENSION I1(50),NXYZ(KWD,3),NUC(NHT2,25),NRC(NHT2,25),
      DIMENSION I1(50),NXYZ(KWD,3),NUC(NHT2,NFUNMX),NRC(NHT2,NFUNMX),
     1 NBAS(8),
     2 CORR(100,3),KA(8,3),C(MXTNCC),JCO(NHT2),PT(8),TITLE(18),
     3 JTRAN(MxCBF),A(50),CTRAN(MxCBF,224),
     4 ICDS(MXTNPR), NDEG(100), CHARGT(100), X(3)
C
C     To hold SOAO matrix pointers for cartesian basis
C
      Integer JTranC(MxCBF), ITranC(MxCBF, 224)
      Double precision CTranC(MxCBF, 224)
C
      DOUBLE PRECISION CJUNK(MxCBF*MxCBF)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD 
      COMMON /REP/ NEWIND(MXCBF) , MSTOLD(8) 
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(100),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(100,3),CHARGE(100),FMULT(8),TLA, TLC
      COMMON /SYMIND/ LABORB(MXCBF)
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA IECPU /35/
      DATA INU /5/, IF10 /0/, IF18 /0/
      DATA HBLNK/'        '/
      DATA JCNTAM /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     & 0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1/
C
      DATA LUGRAD /36/, LUT /6/, LUAUX /62/, LUSCR /13/
C
      DATA NANUM /'       1','       2','       3','       4',
     1            '       5','       6','       7','       8'/
      DATA ISPD /'s', 'p','d','f','g','h','i'/
      DATA KO /'X','Y','Z'/   , IBLANK /' '/
      DATA KOKO /'YZ','XZ','XY'/
      DATA (KWO(I),I=1,56)
     O         /'S   ', 'X   ', 'Y   ', 'Z   ',
     1          'XX  ', 'XY  ', 'XZ  ', 'YY  ',
     2          'YZ  ', 'ZZ  ', 'F300', 'F210',
     3          'F201', 'F120', 'F111', 'F102',
     4          'F030', 'F021', 'F012', 'F003',
     5          'G400', 'G310', 'G301', 'G220',
     6          'G211', 'G202', 'G130', 'G121',
     7          'G112', 'G103', 'G040', 'G031',
     8          'G022', 'G013', 'G004', 'H500',
     9          'H410', 'H401', 'H320', 'H311',
     O          'H302', 'H230', 'H221', 'H212',
     1          'H203', 'H140', 'H131', 'H122',
     2          'H113', 'H104', 'H050', 'H041',
     3          'H032', 'H023', 'H014', 'H005'/
      data (kwo(i),i=57,84) /
     $ 'I600', 'I510', 'I501', 'I420', 'I411', 'I402',
     $ 'I330', 'I321', 'I312', 'I303', 'I240', 'I231',
     $ 'I222', 'I213', 'I204', 'I150', 'I141', 'I132',
     $ 'I123', 'I114', 'I105', 'I060', 'I051', 'I042',
     $ 'I033', 'I024', 'I015', 'I006'/
      DATA NCONTR /NCONTS, NCONTP, NCONTD, NCONTF, NCONTG,
     $ NCONTH, NCONTI/
      data iadq/0,1,10,46,146,371,812,1596/
      data (charm(i),i = 1,146) /1,  1,0,0, 0,1,0, 0,0,1,
     D       -1., 0 , 0 ,-1., 0 , 2.,
     D        0 , 1., 0 , 0 , 0 , 0 ,
     D        0 , 0 , 1., 0 , 0 , 0 ,
     D        1., 0 , 0 ,-1., 0 , 0 ,
     D        0 , 0 , 0 , 0 , 1., 0 ,
     D        1., 0 , 0 , 1., 0 , 1.,
     F       -1., 0 , 0 ,-1., 0 , 4 , 0 , 0 , 0 , 0 ,
     F        0 ,-1., 0 , 0 , 0 , 0 ,-1., 0 , 4 , 0 ,
     F        0 , 0 ,-3., 0 , 0 , 0 , 0 ,-3., 0 , 2.,
     F        1., 0 , 0 ,-3., 0 , 0 , 0 , 0 , 0 , 0 ,
     F        0 , 0 , 0 , 0 , 1., 0 , 0 , 0 , 0 , 0 ,
     F        1., 0 , 0 , 1., 0 , 1., 0 , 0 , 0 , 0 ,
     F        0 , 3., 0 , 0 , 0 , 0 ,-1., 0 , 0 , 0 ,
     F        0 , 0 , 1., 0 , 0 , 0 , 0 ,-1., 0 , 0 ,
     F        0 , 1., 0 , 0 , 0 , 0 , 1., 0 , 1., 0 ,
     F        0 , 0,  1., 0 , 0 , 0 , 0 , 1., 0 , 1./
      data (charm(i),i = 147,371) /
     G  3., 0 , 0 , 6 , 0 ,-24, 0 , 0 , 0 , 0 , 3., 0 ,-24, 0 , 8 ,
     G  0 ,-1., 0 , 0 , 0 , 0 ,-1., 0 , 6 , 0 , 0 , 0 , 0 , 0 , 0 ,
     G  0 , 0 ,-3., 0 , 0 , 0 , 0 ,-3., 0 , 4 , 0 , 0 , 0 , 0 , 0 ,
     G  1., 0 , 0 ,-6 , 0 , 0 , 0 , 0 , 0 , 0 , 1., 0 , 0 , 0 , 0 ,
     G  0 , 0 , 0 , 0 , 3., 0 , 0 , 0 , 0 , 0 , 0 ,-1., 0 , 0 , 0 ,
     G -1., 0 , 0 , 0 , 0 , 6 , 0 , 0 , 0 , 0 , 1., 0 ,-6 , 0 , 0 ,
     G  0 , 1., 0 , 0 , 0 , 0 ,-1., 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     G  0 , 0 , 1., 0 , 0 , 0 , 0 ,-3., 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     G  0 , 1., 0 , 0 , 0 , 0 , 1., 0 , 1., 0 , 0 , 0 , 0 , 0 , 0 ,
     G  0 , 0 , 1., 0 , 0 , 0 , 0 , 1., 0 , 1., 0 , 0 , 0 , 0 , 0 ,
     G  1., 0 , 0 , 2., 0 , 2., 0 , 0 , 0 , 0 , 1., 0 , 2., 0 , 1.,
     G  0 , 0 , 0 , 0 ,-3., 0 , 0 , 0 , 0 , 0 , 0 ,-3., 0 , 4 , 0 ,
     G  1., 0 , 0 , 0 , 0 , 1., 0 , 0 , 0 , 0 ,-1., 0 ,-1., 0 , 0 ,
     G  0 , 0 , 0 , 0 , 1., 0 , 0 , 0 , 0 , 0 , 0 , 1., 0 , 1., 0 ,
     G -1., 0 , 0 ,-2., 0 , 1., 0 , 0 , 0 , 0 ,-1., 0 , 1., 0 , 2./
      DATA (charm(I),I = 372,539) /
     H  1.,0.,0.,2.,0.,-12.,0.,0.,0.,0.,1.,0.,
     H    -12.,0.,8.,0.,0.,0.,0.,0.,0.,
     H  0.,1.,0.,0.,0.,0.,2.,0.,-12.,0.,0.,0.,
     H     0.,0.,0.,1.,0.,-12.,0.,8.,0.,
     H  0.,0.,-1.,0.,0.,0.,0.,0.,0.,2.,0.,0.,
     H     0.,0.,0.,0.,1.,0.,-2.,0.,0.,
     H -1.,0.,0.,2.,0.,8.,0.,0.,0.,0.,3.,0.,
     H    -24.,0.,0.,0.,0.,0.,0.,0.,0.,
     H  0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,-1.,
     H     0.,0.,0.,0.,0.,0.,0.,0.,0.,
     H  1.,0.,0.,2.,0.,2.,0.,0.,0.,0.,1.,0.,
     H     2.,0.,1.,0.,0.,0.,0.,0.,0.,
     H  0.,5.,0.,0.,0.,0.,-10.,0.,0.,0.,0.,0.,
     H     0.,0.,0.,1.,0.,0.,0.,0.,0.,
     H  0.,0.,1.,0.,0.,0.,0.,-6.,0.,0.,0.,0.,
     H     0.,0.,0.,0.,1.,0.,0.,0.,0. /
      DATA (charm(I),I = 540,707) /
     H  0.,-3.,0.,0.,0.,0.,-2.,0.,24.,0.,0.,0.,
     H     0.,0.,0.,1.,0.,-8.,0.,0.,0.,
     H  0.,0.,15.,0.,0.,0.,0.,30.,0.,-40.,0.,0.,
     H     0.,0.,0.,0.,15.,0.,-40.,0.,8.,
     H  1.,0.,0,-10.,0.,0.,0.,0.,0.,0.,5.,0.,
     H     0.,0.,0.,0.,0.,0.,0.,0.,0.,
     H  0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,1.,
     H     0.,1.,0.,0.,0.,0.,0.,0.,0.,
     H  1.,0.,0.,-2.,0.,1.,0.,0.,0.,0.,-3.,0.,
     H    -3.,0.,0.,0.,0.,0.,0.,0.,0.,
     H  0.,0.,0.,0.,-1.,0.,0.,0.,0.,0.,0.,-1.,
     H     0.,2.,0.,0.,0.,0.,0.,0.,0.,
     H -1.,0.,0.,-2.,0.,3.,0.,0.,0.,0.,-1.,0.,
     H     3.,0.,4.,0.,0.,0.,0.,0.,0.,
     H  0.,1.,0.,0.,0.,0.,2.,0.,2.,0.,0.,0.,
     H     0.,0.,0.,1.,0.,2.,0.,1.,0. /
      DATA (charm(I),I = 708,812) /
     H  0.,0.,1.,0.,0.,0.,0.,2.,0.,2.,0.,0.,
     H     0.,0.,0.,0.,1.,0.,2.,0.,1.,
     H  0.,3.,0.,0.,0.,0.,2.,0.,3.,0.,0.,0.,
     H     0.,0.,0.,-1.,0.,-1.,0.,0.,0.,
     H  0.,0.,1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,
     H     0.,0.,0.,0.,-1.,0.,-1.,0.,0.,
     H  0.,-1.,0.,0.,0.,0.,-2.,0.,3.,0.,0.,0.,
     H     0.,0.,0.,-1.,0.,3.,0.,4.,0.,
     H  0.,0.,-3.,0.,0.,0.,0.,-6.,0.,-1.,0.,0.,
     H     0.,0.,0.,0.,-3.,0.,-1.,0.,2. /
      DATA (charm(I),I = 813,980) /
     I   1.,   0.,   0., -15.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  15.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   0.,   0.,   0.,   2.,   0., -16.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   1.,   0., -16.,   0.,  16.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   0.,   0.,   0., -10.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   5.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  -1.,   0.,   0.,   5.,   0.,  10.,   0.,   0.,   0.,   0.,
     I   5.,   0., -60.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,  10.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   5.,   0.,   0.,   0.,   0.,   0.,
     I   0., -10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,
     I   1.,   0.,   0.,   1.,   0., -16.,   0.,   0.,   0.,   0.,
     I  -1.,   0.,   0.,   0.,  16.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,  16.,   0., -16.,   0.,   0. /
      DATA (charm(I),I = 981,1148) /
     I   0.,   6.,   0.,   0.,   0.,   0., -20.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   6.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  -3.,   0.,   0.,   0.,   0.,   6.,   0.,   8.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   9.,   0., -24.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -4.,   0.,   0.,   0.,   0.,   0.,   0.,  40.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   4.,   0., -40.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   0.,   0.,   0.,  -2.,   0.,   1.,
     I   0.,   0.,   0.,   0.,   0.,   0.,  -3.,   0.,  -3.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  -5.,   0.,   0., -15.,   0.,  90.,   0.,   0.,   0.,   0.,
     I -15.,   0., 180.,   0.,-120.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -5.,   0.,  90.,   0.,-120.,   0.,  16.,
     I   0.,   0.,   0.,   0.,  -9.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -6.,   0.,  24.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   3.,   0.,  -8.,   0.,   0.,   0. /
      DATA (charm(I),I = 1149,1316) /
     I   1.,   0.,   0.,  -5.,   0.,   1.,   0.,   0.,   0.,   0.,
     I  -5.,   0.,  -6.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   1.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  20.,   0., -40.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  10.,   0., -40.,   0.,  16.,   0.,
     I  -1.,   0.,   0.,  -1.,   0.,   5.,   0.,   0.,   0.,   0.,
     I   1.,   0.,   0.,   0.,   6.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,  -5.,   0.,  -6.,   0.,   0.,
     I   0.,  -1.,   0.,   0.,   0.,   0.,  -2.,   0.,   5.,   0.,
     I   0.,   0.,   0.,   0.,   0.,  -1.,   0.,   5.,   0.,   6.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  10.,   0.,   0.,   0.,   0.,  20.,   0., -40.,
     I   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0., -40.,   0.,
     I  16.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   0.,   0.,   0.,   0.,   0.,   1.,   0.,
     I   0.,   0.,   0.,   0.,   0.,  -1.,   0.,  -1.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0. /
      DATA (charm(I),I = 1317,1484) /
     I   0.,   0.,  -3.,   0.,   0.,   0.,   0.,  -6.,   0.,   1.,
     I   0.,   0.,   0.,   0.,   0.,   0.,  -3.,   0.,   1.,   0.,
     I   4.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   0.,   0.,   0.,   2.,   0.,   2.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   1.,   0.,   2.,   0.,   1.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   0.,   0.,   0.,   2.,   0.,   2.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   1.,   0.,   2.,   0.,
     I   1.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  -1.,   0.,   0.,  -3.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  -3.,   0.,   0.,   0.,   3.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,   0.,   0.,   3.,   0.,   2.,
     I   0.,   0.,   0.,   0.,   3.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   2.,   0.,   3.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  -1.,   0.,  -1.,   0.,   0.,   0.,
     I   1.,   0.,   0.,   1.,   0.,   2.,   0.,   0.,   0.,   0.,
     I  -1.,   0.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,  -2.,   0.,  -1.,   0.,   0. /
      DATA (charm(I),I = 1485,1596) /
     I   0.,   0.,   0.,   0.,  -3.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -6.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  -3.,   0.,   1.,   0.,   4.,   0.,
     I   3.,   0.,   0.,   9.,   0., -21.,   0.,   0.,   0.,   0.,
     I   9.,   0., -42.,   0., -16.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   3.,   0., -21.,   0., -16.,   0.,   8.,
     I   0.,   0.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   2.,   0.,   2.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   2.,   0.,   1.,   0.,
     I   1.,   0.,   0.,   3.,   0.,   3.,   0.,   0.,   0.,   0.,
     I   3.,   0.,   6.,   0.,   3.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   3.,   0.,   3.,   0.,   1. /
      DATA (KKTYP(I,1),I = 1,84) /'S000', 'P100', 'P010', 'P001',
     1            'D200', 'D110', 'D101', 'D020', 'D011', 'D002',
     2            'F300', 'F210', 'F201', 'F120', 'F111', 'F102',
     3            'F030', 'F021', 'F012', 'F003',
     4            'G400', 'G310', 'G301', 'G220', 'G211', 'G202',
     5            'G130', 'G121', 'G112', 'G130', 'G040', 'G031',
     6            'G022', 'G031', 'G004',
     7            'H500', 'H410', 'H401', 'H320', 'H311', 'H302',
     8            'H230', 'H221', 'H212', 'H203', 'H140', 'H131',
     9            'H122', 'H113', 'H104', 'H050', 'H041', 'H032',
     O            'H023', 'H014', 'H005',
     1            'I600', 'I510', 'I501', 'I420', 'I411', 'I402',
     1            'I330', 'I321', 'I312', 'I303', 'I240', 'I231',
     1            'I222', 'I213', 'I204', 'I150', 'I141', 'I132',
     1            'I123', 'I114', 'I105', 'I060', 'I051', 'I042',
     1            'I033', 'I024', 'I015', 'I006'/
      DATA (KKTYP(I,2),I = 1,84) /'1s  ', '2px ', '2py ', '2pz ',
     1            '3d0 ', '3d2-', '3d1+', '3d2+', '3d1-', '3s  ',
     2            '4f1+', '4f1-', '4f0 ', '4f3+', '4f2-', '4px ',
     3            '4f3-', '4f2+', '4py ', '4pz ',
     4            '5g0 ', '5g2-', '5g1+', '5g4+', '5g3-', '5g2+',
     5            '5g4-', '5g3+', '5d2-', '5d1+', '5s  ', '5g1-',
     6            '5d2+', '5d1-', '5d0' ,
     H '6h1+','6h1-','6h2+','6h3+','6h4-','6px ','6h5-','6h4+','6h3-',
     H   '6h0 ','6h5+','6f2-','6f3+','6h2-','6f1+','6py ','6pz ',
     H   '6f3-','6f2+','6f1-','6f0 ',
     I            '7i6+', '7i2-', '7i5+', '7i4+', '7i5-', '7i2+',
     I            '7i6-', '7i3+', '7i4-', '7g3+', '7i0 ', '7i3-',
     I            '7g4+', '7i1-', '7g2+', '7g2-', '7i1+', '7g4-',
     I            '7g1+', '7d2-', '7d1+', '7d0 ', '7g3-', '7d2+',
     I            '7g1-', '7g0 ', '7d1-', '7s  '/
      DATA IBLNK/'      '/
      DATA NCMAX /100/
C
C...  STATEMENT FUNCTIONS
C
      PTAND(I,J) = PC(AND(I,J))
C
C  OPEN INPUT FILE MOL
C
      OPEN(UNIT=INU,FILE='MOL',FORM='FORMATTED',STATUS='OLD')
      OPEN(UNIT=LUSCR,FILE='VMSCR1',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
      IPRINT = IFLAGS(1)
      IECP=0
      ISTATS=0
      ID3=0 
      TLA=1.D-10 
      SCALE=0.0
      JUNK=1
      NONT=1 
      TLC=0.D0
      NTAP=2

      DO 3335 I = 1, 3
         ISYTYP(I) = 1
 3335 CONTINUE
C
      NM=0
      DO I = 0, NHT2
         DO IX = I, 0, -1
            DO  IY =I-IX, 0, -1
               NM = NM + 1
               NXYZ(NM,1) = IX
               NXYZ(NM,2) = IY
               NXYZ(NM,3) = I-IX-IY
            ENDDO
         ENDDO
      ENDDO
C
      PI=2.D0*DACOS(0.0D00)
      PIPPI=(0.5D0/PI)**0.75D0
C
C-----------------------------------------------------------------------
C     Begin reading the MOL file.
C-----------------------------------------------------------------------
C
      READ(INU,83383) SLASK,IGENER,IF10,IF18,IECP,IECPU,IJUNK,
     &                INUTMP,ISTATS,NHARM,IDELA,IDELB,IBORTH
83383 FORMAT(A6,4X,12I5)
C
      IF (SLASK .EQ. 'SUPMAT') THEN
         WRITE(6,88777) 
         CALL ERREX
      ENDIF
88777 FORMAT(///,'* * * *  THIS CODE DOES NOT GENERATE SUPERMATRICES',
     1         '  * * * *')
C
      LHARM = NHARM .EQ. 0
C
      IF (LHARM) THEN
        IDOSPH = 1
      ELSE
         IDOSPH = 0
         idela=1
      ENDIF
C
      IF (INUTMP .NE. 0) INU = INUTMP
      LPRT=.FALSE.
      IF(IPRINT.GE.1)LPRT = .TRUE.
C
      READ(INU, 833) TITLE
  833 FORMAT (9A8)
C
C-----------------------------------------------------------------------
C Read number of nuclear types (NONTYP), symmetry information
C NSYMOP, KKA), tolerances (TLA, TLC), and coordinate units (ID3).
C IRSTRT is not used anywhere.
C-----------------------------------------------------------------------
      READ(INU,8800) NONTYP, NSYMOP, ((KKA(J,I),I=1,3),J=1,3),
     1              TLA,  TLC, ID3, IRSTRT
8800  FORMAT(2I5,9A1,1X,E10.2,E10.2,I5,I5)
C
      IF (TLC.eq.0.0) TLC=0.1*TLA
      READ(INU,180)TIML,SAFL
 180  FORMAT (2F10.5)
C
C-----------------------------------------------------------------------
C     Process symmetry information and construct symmetry arrays.
C-----------------------------------------------------------------------
C
C MAXLOP = NUMBER OF IRREDUCIBLE REPRESENTATIONS,
C THE LOOPS ARE ALWAYS FROM 1 TO MAXLOP
C
C ISYTYP DETERMINES THE BEHAVIOUR OF X,Y,Z UNDER THE SYMMETRY
C OPERATIONS: BIT 1,2,3 ARE SET TO ONE IF X, OR Y OR Z ARE
C CHANGING UNDER THIS PARTICULAR OPERATIONS.
C !ATTENTION! ISYTYP IS INITIALIZED TO 1 AND NOT TO 0 IN VMOL.
C THIS HAS TO BE CONSIDERED BY USING THE DEFINITION OF ISYTYP
C GIVEN ABOVE!
C EXAMPLE: ISYTYP(2)=2 MEANS THE VALUE OF Y CHANGES UNDER
C          1 SYMMETRYOPERATION.
C GENERATORS:
C NSYMOP: NUMBER OF GENERATORS,
C MAXLOP: NUMBER OF SYMMETRYOPERATIONS (NUMBER OF IRREPS).
C
      MAXLOP=2**NSYMOP
C
      IF(NSYMOP.GT.0) THEN
         DO J=1,NSYMOP
            NBL=1
            DO I=1,3
               IF (KKA(J,I).EQ.IBLANK) NBL=NBL+1
               DO K=1,3
                  IF (KKA(J,I).EQ.KO(K)) ISYTYP(K)=ISYTYP(K)+2**(J-1)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
C
C Determine KA(I,J), PT(I), MULT(I) and FMULT(I)
C KA, PT: fixed
C MULT (FMULT): depend only on the value of NSYMOP
C
      DO I=1,8
         IS = 0
         DO K=1,3
            KA(I,K)=(I-1)/(2**(K-1))-2*((I-1)/(2**K))
            IS=IS+KA(I,K)
         ENDDO
         PT(I)=(-1.0)**IS
         MULT(I)=2**MAX0(0,NSYMOP-IS)
         FMULT(I)=DFLOAT(MULT(I))
      ENDDO
C
      DO I=1,8
         DO J=1,8
            IAN=1
            IOR=1
            DO K=1,3
               IAN=IAN+MIN0(KA(I,K),KA(J,K))*2**(K-1)
               IOR=IOR+MAX0(KA(I,K),KA(J,K))*2**(K-1)
            ENDDO
            AND(I,J)=IAN
            OR(I,J)=IOR
            EOR(I,J)=IOR-IAN+1
         ENDDO
      ENDDO
C
C CREATE ITYPE-ARRAY
C ITYPE(1,1)=1: BEHAVIOUR OF R^0 (NO CHANGE OF SIGN AT ALL!)
C ITYPE(2,1..3)=ISYTYP(1..3): BEHAVIOUR OF R^1 (X, Y, Z)
C AND SO ON UP TO R^6.
C
      ITYPE(1,1)=1
      NDO=0
      NFO=0
      NGO=0
      NHO=0
      NIO=0
C
      DO JQM=1,3
         ITYPE(2,JQM)=ISYTYP(JQM)
         DO KQM=JQM,3
            NDO=NDO+1
            ITYPE(3,NDO) = EOR(ISYTYP(JQM),ISYTYP(KQM))
            DO LQM=KQM,3
               NFO=NFO+1
               ITYPE(4,NFO)=EOR(ITYPE(3,NDO),ISYTYP(LQM))
               DO M2M=LQM,3
                  NGO=NGO+1
                  ITYPE(5,NGO)=EOR(ITYPE(4,NFO),ISYTYP(M2M))
                  DO M3M=M2M,3
                     NHO=NHO+1
                     ITYPE(6,NHO)=EOR(ITYPE(5,NGO),ISYTYP(M3M))
                     DO M4M=M3M,3
                        NIO=NIO+1
                        ITYPE(7,NIO)=EOR(ITYPE(6,NHO),ISYTYP(M4M))
                     ENDDO
                  ENDDO   
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C ENTARTUNG OF L-QUANTUM NUMBERS 1-7
C
      DO JJ=1,7
         KHKT(JJ)=JJ*(JJ+1)/2
      ENDDO
C
C RUN OVER TYPES OF ATOMS IN THE MOLECULE
C
      KMAX=0
      LMAX=0
      LLMAX = 0
      LLMAX = 0
      NMAX=0
      MAXLO=8
      NHTYP=2
C-----------------------------------------------------------------------
C     Begin loop over center (nuclear) types.
C-----------------------------------------------------------------------
      DO I=1,NONTYP
         DO J = 1,5
            IAD(J) = 0
         ENDDO
C
         READ(INU,11) KASET, Q,NONT,IQM,(JCO(J),J=1,IQM), IGENES
 11      FORMAT (A6,F14.8,12I5)
         LUBAS = INU
         IQM4 = MAX0(IQM,1)
         NMIN=NMAX+1
C     
C LOOP OVER SYMMETRY INDEPENDENT ATOMS IN THE MOLECULE
C
         DO N=1, NONT
C
C NMAX IS THE TOTAL NUMBER OF SYMMETRY INDEPENDENT ATOMS
C     
            NMAX = NMAX + 1
            IF (NMAX .LE. NCMAX) GO TO 113
            WRITE(6,1305) NCMAX
 1305       FORMAT ('1',5X,'SPACE LIMITATION ENCOUNTERED IN READIN'
     1         /11X,'NUMBER OF CENTERS EXCEEDS',I4)
            CALL ERREX
C
 113        IFBA(NMAX) = IQM
            READ(INU,12)NAMN(NMAX),(CORD(NMAX,J),J=1,3),IDX
 12         FORMAT (A4,3F20.12,I5)
            IF (IDX .EQ. 0) GO TO 115
            READ(INU,30)(CORD(NMAX,J),J=1,3)
 30         FORMAT (3F20.5)
 115        CONTINUE
C
            IF (ID3 .EQ. 1) THEN
               DO J=1,3
                  CORD(NMAX,J)=CORD(NMAX,J)/0.5291771
               ENDDO
            ENDIF
         ENDDO
C-----------------------------------------------------------------------
C     End loop over symmetry inequivalent centers of each type.
C-----------------------------------------------------------------------
         IMAX=0
         IIMAX = 0
         IF(IQM.EQ.0) GO TO 1272
C-----------------------------------------------------------------------
C     Loop over angular momentum for each nuclear type.
C-----------------------------------------------------------------------
         DO J=1,IQM
            KM=JCO(J) + IAD(J)
C
 13         FORMAT (10I5)
 211        FORMAT(8F10.4)
 213        FORMAT(4F18.4)
 14         FORMAT (8F9.4)
C
C
            IF (JCO(J) .EQ. 0 .AND. IAD(J) .EQ. 0) GOTO 2
C-----------------------------------------------------------------------
C     Loop over shells of each angular momentum. Read number of primit-
C     ives and contracted functions.
C-----------------------------------------------------------------------
            DO K = 1, KM
               IF (IGENER .EQ. 0 .AND. IGENES .EQ. 0) THEN
                  READ(LUSCR) NUC(J,K), NRC(J,K)
               ELSE
                  IF (K .LE. JCO(J)) THEN
                     READ(LUBAS,13) NUC(J,K), NRC(J,K)
                  ELSE
                     READ(INU,13) NUC(J,K), NRC(J,K)
                  ENDIF
               ENDIF
C
               IF(NUC(J,K) .GT. NFUNMX)THEN
               WRITE(6,*)' @READIN-F, Too many primitives in this shell'
               WRITE(6,*) '             Value is ',NUC(J,K)
               WRITE(6,*) '             Limit is ',NFUNMX
               CALL ERREX
               ENDIF
C
               IF(NRC(J,K) .GT. NFUNMX)THEN
               WRITE(6,*)' @READIN-F, Too many cont. fns. in this shell'
               WRITE(6,*) '             Value is ',NRC(J,K)
               WRITE(6,*) '             Limit is ',NFUNMX
               CALL ERREX
               ENDIF
C
               IMIN = IMAX+1
               IMAX = IMAX+NUC(J,K)
C-----------------------------------------------------------------------
C Loop over primitives in each shell, reading exponents and cont-
C raction coefficients.
C-----------------------------------------------------------------------
               DO LL = 1, NUC(J,K)
C
                  IF (IGENER .EQ. 0 .AND. IGENES .EQ. 0) THEN
                     READ(LUSCR) A(IMIN-1+LL), (CTEMP(LL,MC), 
     &                           MC = 1, NRC(J,K))
                  ELSE
                     IF (KASET .NE. '      ') THEN
                        IF (K .LE. JCO(J)) THEN
                           IF (IF10 .EQ. 1) THEN
                              READ(LUBAS,211) A(IMIN-1+LL),
     &                                     (CTEMP(LL,MC),MC=1,NRC(J,K))
                           ELSE
                              READ(LUBAS,14) A(IMIN-1+LL),
     &                                     (CTEMP(LL,MC),MC=1,NRC(J,K))
                           ENDIF
C
C SCALE S FUNCTIONS IF DESIRED
C     
                           IF (J .EQ. 1 .AND. SCALE .NE. 0.) THEN
                              DO L = IMIN, IMAX
                                 A(L) = A(L)*SCALE*SCALE
                              ENDDO
                           ENDIF
                        ELSE
                           IF (IF10 .EQ. 1) THEN
                              READ(INU,211) A(IMIN-1+LL),
     &                                     (CTEMP(LL,MC),MC=1,NRC(J,K))
                           ELSE
                              READ(INU,14) A(IMIN-1+LL), 
     &                                     (CTEMP(LL,MC),MC=1,NRC(J,K))
                           ENDIF
                        ENDIF
                     ELSE
                        IF(IF18.EQ.1)GO TO 88212
                        IF(IF10.EQ.1)GO TO 88210
C
                        READ(INU,14) A(IMIN-1+LL), 
     &                              (CTEMP(LL,MC),MC = 1,NRC(J,K))
                        GO TO 88215
C
88210                   CONTINUE
                        READ(INU,211) A(IMIN-1+LL), 
     &                               (CTEMP(LL,MC),MC = 1,NRC(J,K))
                        GO TO 88215
C
88212                   CONTINUE
                        READ(INU,213) A(IMIN-1+LL), 
     &                                (CTEMP(LL,MC),MC = 1,NRC(J,K))
C
88215                   CONTINUE
                     ENDIF
                  ENDIF
               ENDDO
C
C NOW RENORMALIZE COEFFICIENT ARRAY AND MOVE FROM CTEMP TO C.
C INPUT COEFFICIENTS (POSSIBLY UNNORMALIZED) ARE HELD IN CUN.
C
               IF (IBORTH .EQ. 0) THEN
C
C Loop over contracted functions in each shell.
C
                  DO MC = 1, NRC(J,K)
                     IIMIN = IIMAX
                     SUM = 0.
                     DO LL = 1, NUC(J,K)
                        AL = A(IMIN-1+LL)
                        DO LM = 1, LL
                           AM = A(IMIN-1+LM)
                           T = CTEMP(LL,MC)*CTEMP(LM,MC)*
     &                         (2.0*DSQRT(AL*AM)/(AL+AM))**(J+0.5)
                           SUM = SUM + T
                           IF (LL .NE. LM) SUM = SUM + T
                        ENDDO
                     ENDDO

                     XNORM = 1./DSQRT(SUM)
                     DO LL = 1, NUC(J,K)
                        AL = A(IMIN-1+LL)
                        CUN(IIMIN+LL) = CTEMP(LL,MC)
                        C(IIMIN+LL) = CTEMP(LL,MC)*PIPPI*
     &                                XNORM*(4.*AL)**(0.5*J+0.25)
                     ENDDO
                     IIMAX = IIMAX + NUC(J,K)
C-----------------------------------------------------------------------
C     End of loop over contracted functions.
C-----------------------------------------------------------------------
                  ENDDO
               ELSE
C-----------------------------------------------------------------------
C     Loop over contracted functions in each shell.
C-----------------------------------------------------------------------
                  DO MB = 1, NRC(J,K)
                     IIMIN = IIMAX
                     DO MC = 1, MB
                        SUM = 0.0D0
                        DO LL = 1, NUC(J,K)
                           AL = A(IMIN-1+LL)
                           DO LM = 1, NUC(J,K)
                              AM = A(IMIN-1+LM)
                              T = CTEMP(LL,MB)*CTEMP(LM,MC)*
     &                            (2.0*DSQRT(AL*AM)/(AL+AM))**(J+0.5)
                              SUM = SUM+T
                           ENDDO
                        ENDDO
C
                        IF (MC .eq. MB) THEN
                           DO ll=nuc(j,k),1,-1
                              if(abs(ctemp(ll,mc)).gt.0.1) goto 88202
                           ENDDO
C
                           Write(6,*)' Normalization failed. C<0.1'
                           CALL ERREX
C
88202                      ac=ctemp(ll,mc)/abs(ctemp(ll,mc))
C
                           DO LL = 1, NUC(J,K)
                              AL = A(IMIN-1+LL)
                              ctemp(ll,mc) = ac*ctemp(ll,mc)/sqrt(SUM)
                              CUN(IIMIN+LL) = CTEMP(LL,MC)
                              C(IIMIN+LL)=CTEMP(LL,MC)*PIPPI*
     &                                    (4.*AL)**(0.5*J+0.25)
                           ENDDO
                           IIMAX = IIMAX + NUC(J,K)
                        ELSE
                           DO  LL=1, NUC(J,K)
                              CTEMP(LL,MB) = CTEMP(LL,MB) -
     &                                       SUM*CTEMP(LL,MC)
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
C-----------------------------------------------------------------------
C End of loop over contracted functions.
C-----------------------------------------------------------------------
               ENDIF
C
            ENDDO
C-----------------------------------------------------------------------
C End of loop over shells of each angular momentum.
C-----------------------------------------------------------------------
            IF (IGENER .NE.0.OR.IGENES.NE.0) JCO(J) = JCO(J) + IAD(J)
 2          CONTINUE
C
         ENDDO
C-----------------------------------------------------------------------
C End of loop over angular momentum for each nuclear type.
C-----------------------------------------------------------------------
C
 1272    CONTINUE
C
C-----------------------------------------------------------------------
C Loop over centers of current type, fill arrays with exponents,
C contraction coefficients, symmetry information, and data for each
C shell.
C-----------------------------------------------------------------------
         DO N = NMIN,NMAX
C
C-----------------------------------------------------------------------
C Load up exponents of primitives.
C-----------------------------------------------------------------------
            IF(LMAX+IMAX .GT. MXTNPR)THEN
               WRITE(6,*)' @READIN-F, Too many primitives. Limit is ',
     &                   MXTNPR
               CALL ERREX
            ENDIF
C
            DO L=1,IMAX
               ALPHA(L+LMAX) = A(L)
            ENDDO
            LMAX = LMAX + IMAX
C
C-----------------------------------------------------------------------
C Load up contraction coefficients.
C-----------------------------------------------------------------------
            IF(LLMAX+IIMAX .GT. MXTNCC)THEN
               WRITE(6,*)' @READIN-F, Too many contraction coefficients'
               WRITE(6,*)' The limit is', MXTNCC
               CALL ERREX
            ENDIF
C
            DO L=1,IIMAX
               CONTUN(L+LLMAX) = CUN(L)
               CONT(L+LLMAX) = C(L)
            ENDDO
            LLMAX = LLMAX+ IIMAX
            MULK=1
            LL=1
            DO L=1,NSYMOP
               DO M=1,3
                  IF(AND(LL+1,ISYTYP(M)).NE.1 .AND. 
     &               ABS(CORD(N,M)).GE.1.E-6) GO TO 9

               ENDDO
               MULK=MULK+LL
 9             CONTINUE
               LL=2*LL
            ENDDO
C     
            MAXLO=AND(MULK,MAXLO)
            MULNUC(N)=MULK
            CHARGE(N)=Q
            NHTYP=MAX0(NHTYP,IQM)
C
            IF (.NOT.(IQM.EQ.0))THEN
               DO J=1,IQM
                  JCOJ=JCO(J)
                  DO K=1,JCOJ
                     KMAX=KMAX+1
C
                     IF(KMAX .GT. MXTNSH)THEN
                        WRITE(6,*) ' @READIN-F, Too many shells.'
                        WRITE(6,*) ' Limit is ',MXTNSH
                        CALL ERREX
                     ENDIF
C     
                     NUCO(KMAX)=NUC(J,K)
                     NRCO(KMAX) = NRC(J,K)
                     LAMN(KMAX)=NAMN(N)
                     MUL(KMAX)=MULK
                     NHKT(KMAX)=J
C
                     DO M=1,3
                        CENT(M,KMAX)=CORD(N,M)
                     ENDDO
C
                  ENDDO
               ENDDO
            ENDIF
C
         ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C End of loop over nuclear types. All of MOL file has now been read.
C-----------------------------------------------------------------------
C
      if(idelb .gt. 0) read 828, (kdelb(i), i=1,idelb)
 828  format (25i3)
C
      NPARSU(1)=0
      MST(1) = 0
      NORB=1
      norxb=0
C
C Tracks the number of deleted orbitals and where
C their symmetry transformations are stored
C
      NOrbDel = MxCBF
C     
      IBBAS=0
      DO LA=1,MAXLOP
         IF(LA.EQ.1) IABAS=0
         KABAS=0
         NBAS(LA)=0
         JSTRT(1)=0
         JRS(1) = 0
         NSTRT(1)=0
         DO IA=1,KMAX
            NHKTA=NHKT(IA)
            KHKTA=KHKT(NHKTA)
            IF(IA.NE.KMAX) THEN
               JSTRT(IA+1)=JSTRT(IA)+NUCO(IA)
               JRS(IA+1) = JRS(IA) + NRCO(IA) * NUCO(IA)
               NSTRT(IA+1)=NSTRT(IA)+KHKTA
            ENDIF
            DO NA=1,KHKTA
               MXYZ=((NHKT(IA)+1)*NHKT(IA)*(NHKT(IA)-1))/6+NA
               IFIF = AND(MUL(IA),EOR(LA,ITYPE(NHKT(IA),NA))) .EQ. 1
               DO JAA = 1, NRCO(IA)
                  JKB=0
                  DO KB=1,MAXLOP
                     IF (.NOT.(AND(KB,MUL(IA)) . NE . 1)) THEN
                        JKB=JKB+1
                        KABAS=KABAS+1
C
                        IF(NORB.GT.MXCBF)THEN
                           WRITE(6,*)' @READIN-F, NORB too large ',NORB
                           CALL ERREX
                        ELSEIF(JKB.GT.224)THEN
                           WRITE(6,*)' @READIN-F, JKB too large ',JKB
                           CALL ERREX
                        ENDIF
C
                        CTRAN(NORB,JKB)=
     $                  PC(AND(KB,EOR(LA,ITYPE(NHKT(IA),NA))))
                        ITRAN(NORB,JKB)=KABAS
                        JTranC(NOrb) = JKB
                        CTranC(NOrb, JKB) = CTran(NOrb, JKB)
                        ITranC(NOrb, JKB) = ITran(NOrb, JKB)
                     ENDIF
                  ENDDO
C
                  IF(IFIF) THEN
                     dele=(idela.eq.0).and.(jcntam(na,nhkta).eq.1)
                     if(idelb .gt. 0) then
                        do i=1,idelb
                           dele = dele .or. kdelb(i).eq.norb
                        enddo
                     endif
                     if (dele) then
                        newind(norb)= -NOrbDel
C
                        IF(NORBDEL.GT.MXCBF)THEN
                           WRITE(6,*)' @READIN-F, NORBDEL too large',
     &                                 NORBDEL
                           CALL ERREX
                        ENDIF
C
                        JTRAN(NOrbDel)=JKB
                        JPRX(NOrbDel)=KWO(MXYZ)
                        MAMN(NOrbDel)=LAMN(IA)
C     
C NOrbDel starts at a useful value, so the update is done _after_ the
C value is used.
C     
                        NOrbDel = NOrbDel - 1
                     else
C
                        IF(NORXB.GT.MXCBF)THEN
                           WRITE(6,*)' @READIN-F, NORXB too large ',
     &                                 NORXB
                           CALL ERREX
                        ENDIF
C
                        norxb=norxb+1
                        newind(norb)=norxb
                        JTRAN(NORXB)=JKB
                        JPRX(NORXB)=KWO(MXYZ)
                        MAMN(NORXB)=LAMN(IA)
                     endif
C
                     NORB=NORB+1
                     NBAS(LA)=NBAS(LA)+1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
C
         IF(LA.GT.1) THEN
            NPARSU(LA)=NPARSU(LA-1)+NPAR(LA-1)
            MST(LA) = MST(LA-1) +  NBAS(LA-1)
         ENDIF
C
         NPAR(LA)=NBAS(LA)*(NBAS(LA)+1)/2
      ENDDO
C
      if (lharm) then
         kab=0
         do ia=1,kmax
            khkta=khkt(nhkt(ia))
            khkta=nhkt(ia)*(nhkt(ia)+1)/2
            do nb=1,khkta
               do jaa=1,nrco(ia)
                  do kb=1,maxlop
                     if (and(kb,mul(ia))  .eq.  1) then
                        kab=kab+1
                        kabs(ia,nb,jaa,kb)=kab
                     endif
                  enddo
               enddo
            enddo
         enddo
C
C map 1-dim. seqno onto la,ia,na and jaa
C
         noryb=0
C
         do la=1,maxlop
            do ia=1,kmax
               khkta=nhkt(ia)*(nhkt(ia)+1)/2
               do na=1,khkta
                  do jaa=1,nrco(ia)
                     if (.not.(and(mul(ia),eor(la,itype(nhkt(ia),
     &                    na))) .ne. 1)) then
                        noryb=noryb+1
C     
                        IF(NORYB.GT.MXCBF)THEN
                           WRITE(6,*)'@READIN-F, NORYB too large ',
     &                                 NORYB
                           CALL ERREX
                        ENDIF
C     
                        last(noryb)=la
                        iast(noryb)=ia
                        nast(noryb)=na
                        jaast(noryb)=jaa
                     endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
c     
c find connection between SO's and CGTO's
C NOTE: NOrb is one too _large_ from the way it is updated.  Hence
C the limit norb-1.  The value of NOrb is corrected later on.
C
         do i=1,norb-1
            inew=newind(i)
C
            IF(INEW.GT.MXCBF)THEN
               WRITE(6,*) ' @READIN-F, INEW too large ',INEW
               CALL ERREX
            ENDIF
C
C Instead of skiping deleted functions, construct their
C transformations.  Remember they live at the very end of
C the JTran/CTran/ITran arrays, so they won't harm anything.
C
            if(inew .le. 0) inew = -inew
            jkb=0
            la=last(i)
            ia=iast(i)
            khkta=khkt(nhkt(ia))
            na=nast(i)
            jaa=jaast(i)
            iaq=iadq(nhkt(ia))+(na-1)*khkta
            inx=(nhkt(ia)+1)*(nhkt(ia)-1)*nhkt(ia)/6 + na
C
            do kb=1,maxlop
               if (.not.(and(kb,mul(ia)) .ne. 1)) then
                  p=pc(and(kb,eor(la,itype(nhkt(ia),na))))
                  do nb=1,khkta
                     ch=charm(iaq+nb)
                     if (.not.(abs(ch) .lt. 1.e-03)) then
                        kab=kabs(ia,nb,jaa,kb)
                        jkb=jkb+1
                        ctran(inew,jkb)=p*ch
                        itran(inew,jkb)=kab
                     endif
                  enddo
               endif
            enddo
C
            jtran(inew)=jkb
            jjprx(inew) = kktyp(inx,1+idosph)
C
         enddo
C
      endif
C
CJDW 6/6/95. Add Ajith's block of PUTRECs.
C
C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  &                                                               &
C  & Added by ajith 03/94. Information required to symmetry adapt  &
C  & nuclear perturbations.                                        &
C  &                                                               &
C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C      
      CALL PUTREC (20, 'JOBARC', 'MULCINP ', 100, MULNUC)
      CALL PUTREC (20, 'JOBARC', 'SYMTYPE ', 3 , ISYTYP)
      CALL PUTREC (20, 'JOBARC', 'NOSYMMOP', 1 , NSYMOP)
      CALL PUTREC (20, 'JOBARC', 'NOSYATOM', 1 , NMAX) 
      CALL PUTREC (20, 'JOBARC', 'ATOMSYMB', 300, NAMN)
C
c      IAB=0
c      KNTR=0
c      JBBAS = 0
c      DO IA=1,KMAX
c         NHKTA=NHKT(IA)
c         NAK=(NHKTA*(NHKTA+1)*(NHKTA-1))/6
c         KHKTA=KHKT(NHKTA)
c         MULA=MUL(IA)
c         LAM=LAMN(IA)
c         JRSIA = JRS(IA)
c         JSTRIA = JSTRT(IA)
c         NSTRIA = NSTRT(IA)
c         DO 1141 NA=1,KHKTA
c            NUCA=NUCO(IA)
c            NRCA = NRCO(IA)
c            NAKNA=NAK+NA
c            NAKO=KWO(NAKNA)
C
C CONSTRUCT PSEUDO-SEGMENTED BASIS FOR PROPERTIES AND ECP
C PROGRAMS BY DUPLICATING PRIMITIVES IN CONTRACTED FUNCTIONS.
C
c            DO JAA = 1,NRCA
c               DO KB=1,MAXLOP
c                  IF (.NOT. (AND(KB,MULA).NE.1)) THEN
c                     KBBAS = 0
c                     DO JA = 1,NUCA
c                        EX = ALPHA(JSTRIA+JA)
c                        CON = CONTUN(JRSIA+JA+(JAA-1)*NUCA)
c                        IF (ABS(CON) .GT. 1.E-8) THEN
c                           JBBAS = JBBAS + 1
c                           KBBAS = KBBAS + 1
c                           WRITE(LUSCR) LAM, KB, NAKO, EX, CON
c                        ENDIF
c                     ENDDO
c                  KNTR = KNTR + 1
c                  ICDS(KNTR) = KBBAS
c               ENDDO
c            ENDDO
C
c            DO KB =1, MAXLOP
c               IF(AND(KB,MULA).NE.1) GO TO 1142
c                DO 1144 JA=1,NUCA
c                  IAB=IAB+1
c                  IOVERP = ' '
c                  IF (JA .EQ. 1) IOVERP = '+'
c 1144          CONTINUE
c            ENDDO
c         ENDDO
c      ENDDO
C
      ji=0
      norxb  = 0
      do i=1,maxlop
         nxbas(i) = 0
         do j=1,nbas(i)
            ji=ji+1
            if (.not.(newind(ji) .le. 0)) Then
               nxbas(i) = nxbas(i) + 1
               norxb = norxb+1
               LABORB(NORXB) = I*1024 + NXBAS(I)
            endif
         enddo
      enddo
C
C After the actual transformation, write the rest of the
C matrix (deleted orbitals) necessary to make it square.
C Earlier we placed these at the very end of the arrays.
C    
C NOTE: NOrbDel is one too _small_ because if the way it is updated.
C Hence, the limit NOrbDel+1.
C         
      NORB = NORB - 1 
      CALL ZERO(CJUNK,NORB*NORB)
      DO I=1,NORXB
         J=JTRAN(I)
         DO K=1,J
            CJUNK( (I-1)*NORB +  ITRAN(I,K)) = CTRAN(I,K)
         ENDDO
      ENDDO
C
      IF(NORXB.LT.NORB)THEN
         IOLD=MXCBF
         DO I=NORXB+1,NORB
            J = JTRAN(IOLD)
            DO K=1,J
               CJUNK( (I-1)*NORB + ITRAN(IOLD,K)) = CTRAN(IOLD,K)
            ENDDO
            IOLD = IOLD-1
         ENDDO
      ENDIF
C
      CALL PUTREC(20,'JOBARC','FULLSOAO',NORB*NORB*IINTFP,CJUNK)
cSSS      WRITE(6,*) '  @READIN-I, SYMTRANS matrix '
cSSS      CALL OUTPUT(CJUNK,1,NORB,1,NORB,NORB,NORB,1)
C
C This is the full cartesian basis AO -> SO transformation.
C Include in this info the number of function in each irrep too,
C because there is nowhere else to get it.
C NOTE: Write this record ONLY if it is different from the
C computational basis.
C
      If ( LHarm .OR. NOrb .ne. Norxb) then
C
         DO 3022 I = 1, NOrb
            J = JTranC(I)
 3022    Continue
      EndIf
C
      CALL ZERO(CJUNK,NORB*NORB)
      DO I=1,NORB
         J=JTranC(I)
         DO K=1,J
            CJUNK( (I-1)*NORB + ITranC(I,K) ) = CTranC(I,K)
         ENDDO
      ENDDO
      
cSSS      WRITE(6,*) '  @READIN-I, CSYMTRAN matrix '
cSSS      CALL OUTPUT(CJUNK,1,NORB,1,NORB,NORB,NORB,1)
C
C Write some of the Cartesian stuff to JOBARC.
C
      CALL PUTREC(20,'JOBARC','NAOBASFN',      1,NORB)
C
C     Population of symmetry blocks in full SO space.
C
      CALL PUTREC(20,'JOBARC','FAOBASIR', MAXLOP,NBAS)
C
      CALL PUTREC(20,'JOBARC','CSYMTRAN',NORB*NORB*IINTFP,CJUNK)

      OPEN(UNIT=1, FILE='LABBASIS',STATUS='NEW')
      DO I = 1, NORXB
         WRITE(1,'(3X, A4,3X,A4)') MAMN(I), JPRX(I)
      ENDDO
      CLOSE(UNIT=1, STATUS='KEEP')

C
      RETURN
      END
