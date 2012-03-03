void QuickSort(ElementType A[], int N)
{
   int pivotIndex;

   if (N > 1)
   {
      pivotIndex = Partition(A, N);

      QuickSort(A, pivotIndex);

      QuickSort(&A[pivotIndex+1], (N-pivotIndex-1));
   }
}

int Partition(ElementType A[], int N)
{
   int i, j, incr = 0, decr = 1, swap;
   ElementType Tmp;

   i = 0;
   j = N-1;

   while (i != j)
   {
      if (A[i] > A[j])
      {
         Tmp = A[i];
         A[i] = A[j];
         A[j] = Tmp;

         swap = incr;
         incr = decr;
         decr = swap;
      }

      i += incr;
      j -= decr;
   }

   return j;
}
