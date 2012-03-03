
      include 'maxdim.h'

      integer mx_prefetch_context
      parameter (mx_prefetch_context = 1000)

      integer context
      integer nptr_context, ptr_context
      logical prefetch_flag
      common /prefetch/context(mx_array_index,mx_prefetch_context),
     *                 nptr_context, ptr_context(mx_prefetch_context),
     *                 prefetch_flag
