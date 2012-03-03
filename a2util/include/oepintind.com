c
c This file defines two-electron integral index statement functions of
c various kinds.  It must be included right after the variables
c declarations and before the first executable statement, since it
c contains the variable declarations, statement function definitions and
c variable definitions (that are executable statements)
c
c For the operators to work, two macros N_BAS and N_OCC have to be
c defined to the names of the variables that give the number of basis
c functions and the number of occupied orbitals correspondingly. So a
c routine that includes these definitions must have these variables
c declared and defined. This piece of code also affects the caller name
c space since it needs the statement function dummy arguments to be
c declared. The implementation could be much clear if Fortran allowed
c the external functions to be specified as inline. Most of the modern
c compiler will do that automatically during the optimization step, but
c it is not guaranteed. 
c
c See also the documentation file OEPINTIND.TEX for the description of
c the offset operators for the two-electron integral indices.
c 
c Igor Schweigert, Jan 2004
c $Id: oepintind.com,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
#ifdef DEFINE_OPERATORS
c
c     Include machine-dependent definitions for VMOL-related operators.
c     
#include <machsp.com>
c
c     Declare the operators and the dummy variables. This will affect
c     the namespace of the caller, so be careful.
c     
      integer
     &     i, j, k, l, n_a, n_p, n_h,
     &     i_aa, i_ma, i_ah, i_pa, i_hh, i_pp, i_ph,
     &     n_aa, n_ma, n_ah, n_pa, n_hh, n_pp, n_ph,
     &     i_aaah, i_paah, i_ppah, i_phah, 
     &     n_aaah, n_paah, n_ppah, n_phah, 
     &     i_ppph, i_pphh, i_phph, i_phhh,
     &     n_ppph, n_pphh, n_phph, n_phhh,
     &     i1_aa, i1_pp, i1_hh, i1_ma, i1_ah, i1_pa, i1_ph,
     &     i2_aa, i2_pp, i2_hh, i2_ma, i2_ah, i2_pa, i2_ph,
     &     i1_aaah, i1_paah, i1_ppah, i1_phah, 
     &     i2_aaah, i2_paah, i2_ppah, i2_phah,
     &     i3_aaah, i3_paah, i3_ppah, i3_phah, 
     &     i4_aaah, i4_paah, i4_ppah, i4_phah,
     &     i1_ppph, i1_pphh, i1_phph, i1_phhh,
     &     i2_ppph, i2_pphh, i2_phph, i2_phhh,
     &     i3_ppph, i3_pphh, i3_phph, i3_phhh,
     &     i4_ppph, i4_pphh, i4_phph, i4_phhh,
     &     iupki, iupkj, iupkk, iupkl,
     &     i_aa_cr, i_pp_cr, i_hh_cr,
     &     i_aaah_cr, i_ppah_cr,
     &     i_ppph_cr, i_pphh_cr, i_phph_cr, i_phhh_cr

c
c     Define the operators that pack two orbital indices into the
c     combined index.
c     
      i_aa (i, j) = i + (j*(j-1))/2
      i_hh (i, j) = i + (j*(j-1))/2
      i_pp (i, j) = i + (j*(j-1))/2
c
      i_ma (i, j) = i + n_a * (j-1)
      i_ah (i, j) = i + n_a * (j-1)
      i_pa (i, j) = i + n_p * (j-1)
      i_ph (i, j) = i + n_p * (j-1)
c
c     Define the operators that pack four orbital indices into the
c     combined index.
c 
      i_aaah (i, j, k, l) = i_aa (i, j) + n_aa * (i_ah (k, l) - 1)
      i_paah (i, j, k, l) = i_pa (i, j) + n_pa * (i_ah (k, l) - 1)
      i_ppah (i, j, k, l) = i_pp (i, j) + n_pp * (i_ah (k, l) - 1)
      i_phah (i, j, k, l) = i_ph (i, j) + n_ph * (i_ah (k, l) - 1)
      i_ppph (i, j, k, l) = i_pp (i, j) + n_pp * (i_ph (k, l) - 1)
      i_pphh (i, j, k, l) = i_pp (i, j) + n_pp * (i_hh (k, l) - 1)
      i_phph (i, j, k, l) = i_aa (i_ph (i, j), i_ph (k, l))
      i_phhh (i, j, k, l) = i_ph (i, j) + n_ph * (i_hh (k, l) - 1)
c
c     Define the operators that pack two orbital indices into the
c     combined index. These operators ensure the proper range of
c     indices, and hence a bit slower. Some operators do not impose the
c     inequality condition on their arguments, so they do not have the
c     "_cr" versions.
c
      i_aa_cr (i, j) = min(i,j) + (max(i,j)*(max(i,j)-1))/2
      i_hh_cr (i, j) = min(i,j) + (max(i,j)*(max(i,j)-1))/2
      i_pp_cr (i, j) = min(i,j) + (max(i,j)*(max(i,j)-1))/2
      i_aaah_cr (i, j, k, l) = i_aa_cr (i, j) + n_aa * (i_ah (k, l) - 1)
      i_ppah_cr (i, j, k, l) = i_pp_cr (i, j) + n_pp * (i_ah (k, l) - 1)
      i_ppph_cr (i, j, k, l) = i_pp_cr (i, j) + n_pp * (i_ph (k, l) - 1)
      i_pphh_cr (i, j, k, l) = i_pp_cr (i, j) + n_pp * (i_hh_cr(k,l)- 1)
      i_phph_cr (i, j, k, l) = i_aa_cr (i_ph (i, j), i_ph (k, l))
      i_phhh_cr (i, j, k, l) = i_ph (i, j) + n_ph * (i_hh_cr (k, l) - 1)
c     
c     Define the operators that unpack a combined index into the two
c     orbital indices.
c     
      i2_aa (i) = 1 + (-1 + int (dsqrt (8.d0*i+0.999d0)))/2
      i2_pp (i) = 1 + (-1 + int (dsqrt (8.d0*i+0.999d0)))/2
      i2_hh (i) = 1 + (-1 + int (dsqrt (8.d0*i+0.999d0)))/2
      i2_ma (i) = (i-1) / n_a + 1
      i2_ah (i) = (i-1) / n_a + 1
      i2_pa (i) = (i-1) / n_p + 1
      i2_ph (i) = (i-1) / n_p + 1
      i1_aa (i, j) = i - (j*(j-1))/2
      i1_pp (i, j) = i - (j*(j-1))/2
      i1_hh (i, j) = i - (j*(j-1))/2
      i1_ma (i, j) = i - n_a * (j-1)
      i1_ah (i, j) = i - n_a * (j-1)
      i1_pa (i, j) = i - n_p * (j-1)
      i1_ph (i, j) = i - n_p * (j-1)
c
c     Define the operators that unpack a combined index into the four
c     orbital indices.
c     
      i4_aaah (i) = i2_ah ((i-1)/n_aa + 1)
      i4_paah (i) = i2_ah ((i-1)/n_pa + 1)
      i4_ppah (i) = i2_ah ((i-1)/n_pp + 1)
      i4_phah (i) = i2_ah ((i-1)/n_ph + 1)
      i4_ppph (i) = i2_ph ((i-1)/n_pp + 1)
      i4_pphh (i) = i2_hh ((i-1)/n_pp + 1)
      i4_phph (i) = i2_ph (i2_aa (i) )
      i4_phhh (i) = i2_hh ((i-1)/n_ph + 1)
c 
      i3_aaah (i, j) = i1_ah ((i-1)/n_aa + 1, j)
      i3_paah (i, j) = i1_ah ((i-1)/n_pa + 1, j)
      i3_ppah (i, j) = i1_ah ((i-1)/n_pp + 1, j)
      i3_phah (i, j) = i1_ah ((i-1)/n_ph + 1, j)
      i3_ppph (i, j) = i1_ph ((i-1)/n_pp + 1, j)
      i3_pphh (i, j) = i1_hh ((i-1)/n_pp + 1, j)
      i3_phph (i, j) = i1_ph (i2_aa (i),  j)
      i3_phhh (i, j) = i1_hh ((i-1)/n_ph + 1, j)
c 
      i2_aaah (i, j, k) = i2_aa (i - n_aa * (i_ah (j, k) - 1))
      i2_paah (i, j, k) = i2_pa (i - n_pa * (i_ah (j, k) - 1))
      i2_ppah (i, j, k) = i2_pp (i - n_pp * (i_ah (j, k) - 1))
      i2_phah (i, j, k) = i2_ph (i - n_ph * (i_ah (j, k) - 1))
      i2_ppph (i, j, k) = i2_pp (i - n_pp * (i_ph (j, k) - 1))
      i2_pphh (i, j, k) = i2_pp (i - n_pp * (i_hh (j, k) - 1))
      i2_phph (i, j, k) = i2_ph (i1_aa (i, i_ph (j,k)))
      i2_phhh (i, j, k) = i2_ph (i - n_ph * (i_hh (j, k) - 1))
c 
      i1_aaah (i, j, k, l) = i1_aa (i - n_aa * (i_ah (k, l) - 1), j)
      i1_paah (i, j, k, l) = i1_pa (i - n_pa * (i_ah (k, l) - 1), j)
      i1_ppah (i, j, k, l) = i1_pp (i - n_pp * (i_ah (k, l) - 1), j)
      i1_phah (i, j, k, l) = i1_ph (i - n_ph * (i_ah (k, l) - 1), j)
      i1_ppph (i, j, k, l) = i1_pp (i - n_pp * (i_ph (k, l) - 1), j)
      i1_pphh (i, j, k, l) = i1_pp (i - n_pp * (i_hh (k, l) - 1), j)
      i1_phph (i, j, k, l) = i1_ph (i1_aa (i, i_ph (k, l)), j)
      i1_phhh (i, j, k, l) = i1_ph (i - n_ph * (i_hh (k, l) - 1), j)
c
c     Define the operators that unpack the VMOL index into the four 2e
c     integral indices.
c
      iupki (i) = iand (i,ialone)
      iupkj (i) = iand (ishft(i,-ibitwd),ialone)
      iupkk (i) = iand (ishft(i,-2*ibitwd),ialone)
      iupkl (i) = iand (ishft(i,-3*ibitwd),ialone)
#undef DEFINE_OPERATORS
#endif /* DEFINE OPERATORS */
#ifdef DEFINE_VARIABLES
c
c     Define the sizes of arrays. Note these definitions rely on two
c     macros that has to be defined, N_BAS and N_OCC
c     
      n_a = N_BAS
      n_h = N_OCC
      n_p = n_a - n_h
c     
      n_aa = i_aa (n_a, n_a)
      n_pp = i_pp (n_p, n_p)
      n_hh = i_hh (n_h, n_h)
      n_ma = i_ma (n_a, n_a)
      n_ah = i_ah (n_a, n_h)
      n_pa = i_pa (n_p, n_a)
      n_ph = i_ph (n_p, n_h)
c 
      n_aaah = i_aaah (n_a, n_a, n_a, n_h)
      n_paah = i_paah (n_p, n_a, n_a, n_h)
      n_ppah = i_ppah (n_p, n_p, n_a, n_h)
      n_phah = i_phah (n_p, n_h, n_a, n_h)
      n_ppph = i_ppph (n_p, n_p, n_p, n_h)
      n_pphh = i_pphh (n_p, n_p, n_h, n_h)
      n_phph = i_phph (n_p, n_h, n_p, n_h)
      n_phhh = i_phhh (n_p, n_h, n_h, n_h)
c
#undef DEFINE_VARIABLES
#endif /* DEFINE_VARIABLES */

