
c This routine returns the number of elements and the offset of the first
c element that a numbered process would control in a parallel environment.
c The main purpose of this routine is to localize the load-balancing algorithm.

c INPUT
c int iNode  : the node ID (from 0) whose statistics are to be measured
c int nNodes : the total number of nodes in the workspace
c int iRoot  : the node ID of the virtual root (who gets extra elements)
c int nTotal : the total number of elements to divvy up

c OUTPUT
c int iOff : the offset of the first element in the batch
c int nEls : the number of elements in the batch

      subroutine paces_batch_stat(iNode,nNodes,iRoot,nTotal,iOff,nEls)
      implicit none

c ARGUMENTS
      integer iNode, nNodes, iRoot, nTotal, iOff, nEls

c INTERNAL VARIABLES
      integer iFloor, iLeft
      integer iNode2
      integer iTmp

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert 0 <= (iNode|iRoot) < nNodes
      if ((nNodes.le.iNode).or.(iNode.lt.0).or.
     &    (nNodes.le.iRoot).or.(iRoot.lt.0)    ) then
         print *, '@PACES_BATCH_STAT: Assertion failed.'
         print *, '   nNodes = ',nNodes
         print *, '   iNode  = ',iNode
         print *, '   iRoot  = ',iRoot
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif

c ----------------------------------------------------------------------

      if (nTotal.lt.0) then
         print *, '@PACES_BATCH_STAT: The total number of elements ',
     &            'must be >= 0.'
         print *, '                   nTotal = ',nTotal
         call aces_exit(1)
      end if

      iFloor = nTotal/nNodes
      iLeft  = mod(nTotal,nNodes)

c   o the sexy version
c      iNode2 = mod(iNode+nNodes-iRoot,nNodes)
c      iOff   = iFloor*iNode2 + min(iNode2,iLeft)
c      nEls   = iFloor + min(1,(iLeft/(1+iNode2)))

c   o the naive version is a lot faster w/o the mod and integer divide
      iNode2 = iNode-iRoot
      if (iNode.lt.iRoot) iNode2 = iNode2+nNodes
      nEls = iFloor
      iOff = iFloor*iNode2 + iLeft
      if (iNode2.lt.iLeft) then
         nEls = iFloor+1
         iOff = nEls*iNode2
      end if

      return
c     end subroutine paces_batch_stat
      end

