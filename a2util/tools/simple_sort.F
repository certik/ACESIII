C
      Subroutine simple_sort(a, b, n)
C      
C This is a very simple routine that can sort eigenvalues and
C eigenvectors (or any other such combination) in ascending order. 
C The sorted eigenvalues and vectors are returned in the respective
C incomming arrays a and b. 11/09, Ajith Perera. 
C

      Implicit None
      Integer i, j, k, n
      Double Precision a(n), b(n,n), dtmp
C
      do i = 1, n-1
            do j = i+1, n
               if ( a(i) .gt. a(j) ) then
                  dtmp = a(i)
                  a(i) = a(j)
                  a(j) = dtmp
                  do k = 1, n
                     dtmp   = b(k,i)
                     b(k,i) = b(k,j)
                     b(k,j) = dtmp
                  end do
               end if 
            end do 
         end do 
C
      Return
      End

