      SUBROUTINE TRNLST(IRREP, LISTG, DISSIZE, NUMDIS, SCR, MAXSIZE)
C
C This routine transpose a ACES II list and write it back to disk. 
C This version handle both in and out-of-core transpositions and 
C work with both square and rectangular matrices (NMO < NAO). The 
C original implementation (JG and JFS) was limited to square matrices,
C and generalized version was implemented by Ajith Perera 11/2001. 
C 
C Implementation Notes: In the case of NMO < NAO (DROPMO, 
C SPHERICAL=ON) the DISSIZE  > NUMDIS. What is done here 
C is to add the appropriate numeber of empty distributions to
C make DISSIZE=NUMDIS. With that both in-core and out-of-core
C algorithms can proceed painlessly.  
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER DISSIZE, DISLEFT, DISREAD
      DIMENSION SCR(MAXSIZE)
C
C Determine the maximum number of distributions which can be held 
C in available memory.
C
      MAXDIS = MAXSIZE/DISSIZE
C
C Check to see if we can fit the whole list into available memory.
C
      IF (MAXDIS .GE. DISSIZE) THEN
C
          CALL ZERO(SCR, DISSIZE*DISSIZE)
C
C Do it in core by using mtran2 (in place in transpose routine)
C Notice the dimension in mtran2 is dissize (not numdis). This is
C to make a square matrix with the padded zeros. 
C
         CALL GETLST(SCR, 1, NUMDIS, 1, IRREP, LISTG)
         CALL MTRAN2(SCR, DISSIZE, DISSIZE) 
         CALL PUTLST(SCR, 1, DISSIZE, 1, IRREP, LISTG)
C
      ELSE
C 
C General out-of-core algorithm. Pad the LISTG with empty
C distributions for cases where DISSIZE not equal to NUMDIS.
C
          IF (DISSIZE .NE. NUMDIS) THEN
              NPADD = DISSIZE - NUMDIS
              IF (NPADD .GT. MAXDIS) CALL INSMEM('@TRNLST',
     &            DISSIZE*NPADD, MAXSIZE)
              CALL ZERO(SCR, DISSIZE*NPADD)
              CALL PUTLST(SCR, NUMDIS+1, NPADD, 1, IRREP,
     &                    LISTG) 
          ENDIF
C
C Leave a length of one distribution for buffer space. 
C
         MAXDIS = MAXDIS - 1
C
         IF (MAXDIS .LE. 0 ) THEN
            WRITE(6,1000) 2*DISSIZE
 1000       FORMAT(T3,'@TRNLST-F, Less than ',I10,' words available',
     &        ' for transposition')
            CALL ERREX
         ENDIF
C
C Compute offset for the buffer area
C
         IBUF = 1 + MAXDIS*DISSIZE
C
C Set the number of distributions need to be treated.
C
         DISREAD = DISSIZE
C
C Set the number of distributions left.
C
         DISLEFT = DISSIZE
C
C Set OFFSET for list LISTG 
C 
         IOFFSET = 1
         IPASS   = 0
         IREAD   = 0
C
 10      CONTINUE

         IPASS = IPASS + 1
         IREAD = IREAD + DISSIZE - IOFFSET + 1
C
C Determine the number of distributions processed during this pass.
C
         DISREAD = MIN(DISLEFT, MAXDIS)
         DISLEFT = DISLEFT - DISREAD
C 
C Read the appropriate number distributions for this path.
C
         CALL GETLST(SCR, IOFFSET, DISREAD, 1, IRREP, LISTG)
C
C Set IFIRST (first distribution in core)
C
         IFIRST = IOFFSET
C
C Update the offset
C
         IOFFSET = IOFFSET + DISREAD
C
C Set ILAST (last distribution in core)
C
         ILAST = IOFFSET - 1
C
C We have now three different categories of distributions 
C
C 1. The distributions in memory 
C 2. The distributions on disk which are already the transpose of the 
C    original matrix. 
C 3. The distributions on disk which are still not the transpose of the
C    original matrix 
C
C First, do an in-place transposition of the elements that are already 
C in memory. 
C
         CALL MTRAN3(SCR(IFIRST), DISSIZE, DISREAD)
C
C Loop over the remaining distributions 
C
         DO 100 IDIS = IOFFSET, DISSIZE 

            CALL GETLST(SCR(IBUF), IDIS, 1, 1, IRREP, LISTG)
C
C  Move A(I,J) [i<j] in main memory to A(i,j) in buffer. 
C
            DO 150 IPOS = IFIRST, ILAST
               X = SCR(IBUF+IPOS-1)
               SCR(IBUF + IPOS - 1) = SCR(IDIS + (IPOS-IFIRST)*DISSIZE)
               SCR(IDIS + (IPOS-IFIRST)*DISSIZE) = X
 150        CONTINUE
C
C Write the distribution in the buffer to the disk. 
C
            CALL PUTLST(SCR(IBUF), IDIS, 1, 1, IRREP, LISTG)
C
 100     CONTINUE
C
C Also, write the transposed main area of the memory to the disk.
C
      CALL PUTLST(SCR, IFIRST, DISREAD, 1, IRREP, LISTG)
C
C Continue untill we exhaust all the distributions 
C
         IF (DISLEFT .NE. 0) GO TO 10

cSSS      WRITE(6,2000)IPASS,IREAD
cSSS 2000  FORMAT(T3,'@TRNLST-I, Transpostion required ',I5,' passes and ',
cSSS     &       I10,' reads.')
cSSS
      ENDIF
C
      RETURN
      END
