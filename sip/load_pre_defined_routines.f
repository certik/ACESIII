C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine load_pre_defined_routines()
      implicit none

      integer dummy, load_user_sub
      external sip_barrier
      external energy_denominator
      external energy_product
      external energy_reg_denominator
      external energy_reg_product
      external energy_denominator_reg_deriv
      external init_c
      external init_eps
      external print_scalar
      external trace_on, trace_off
      external load_balance_on, load_balance_off
      external dump_block
      external blocks_to_list
      external list_to_blocks
      external der_int_setup
      external compute_derivative_integrals
      external fmult 
      external server_barrier
      external scont1, hcont1, dcont2
      external eigen_calc, eig_inv, eig_sr, eig_sr_inv
      external copy_fock
      external diis_setup 
      external compute_diis
      external write_blocks_to_list
      external read_list_to_blocks
      external remove_diagonal
      external return_diagonal
      external fock_denominator
      external set_flags
      external set_flags2
      external der2_comp
      external fock_der
      external overlap_der
      external scontxy
      external hcontxy
      external compute_sderivative_integrals
      external read_hess
      external return_h1
      external smon_on
      external smon_off
      external scf_rhf
      external array_copy
      external c1_print, c1b_print, c2aa_print, c2ab_print, c2bb_print

      external comp_ovl3c,dump_amp,open_amp,copy_ff,copy_ba,copy_ab,
     *         remove_dx, remove_sx, remove_sd, removevv_dd,
     *         removeoo_dd, remove_dd, remove_ss, remove_ds,
     *         remove_xd, remove_xs, remove_xxdd, remove_xdxd,
     *         remove_xxsd, remove_xsxd, remove_xdxs, remove_xxss,
     *         remove_xsxs, remove_dddd, remove_ssss,
     *         remove_single_double, remove_double, remove_single, 
     *         set_index, udenominator, read_grad, 
     *         energy_adenominator, energy_bdenominator,
     *         energy_abdenominator
      external check_dconf
      external dropmo_expand_basis
      external square_root, norm_fac, return_sval, return_diagonal4,
     *         symm_force_a, symm_force_i,
     *         place_sval,
     *         place_one,
     *         place_one2,
     *         place_one3,
     *         place_one4,
     *         place_one5,
     *         place_one6,
     *         place_oneb,
     *         place_one2b,
     *         place_one3b,
     *         place_one4b,
     *         place_one5b,
     *         place_one6b,
     *         smooth,smooth4,
     *         eigen_nonsymm_calc
      external 
     *         eomroot_print,  
     *         eomroot_print_new, 
     *         place_one_dip,
     *         place_one_dip_2,
     *         place_one_dip_3,
     *         place_one_dip_4,
     *         place_one_dip_5,
     *         place_one_dip_6,
     *         place_one_dip_7,
     *         place_one_dip_8,
     *         place_one_dea,
     *         place_one_dea_2,
     *         place_one_dea_3,
     *         place_one_dea_4,
     *         place_one_dea_5,
     *         place_one_dea_6,
     *         place_one_dea_7,
     *         place_one_dea_8
 
      external apply_den2, apply_den4, apply_den4_nodiag
      external energy_tdenominator
      external compute_aaaa_batch, compute_aaab_batch, 
     *         compute_aabb_batch, compute_aabc_batch,
     *         compute_abab_batch, compute_abac_batch,
     *         compute_abcd_batch, compute_no4c_batch,
     *         remove_atom_rud1, remove_atom_rud2
      external der4_comp, set_flags4
      external timestamp
      external pardo_sects, get_ijk, set_ijk_aaa, stripi,
     *         sum_64ss, set_ijk_aab
      external get_my_rank, broadcast_array
      external checkpoint, get_restart_status, commit_checkpoint
      external crash
C
C
C             ...Watson: XYZ moment integrals
C
      external print_eom_dens_info
      external return_x,  return_y,  return_z
      external return_xx, return_yy, return_zz
      external return_xy, return_xz, return_yz
      external return_derv_xyz, return_ovl,
     *         dipole_moment, nuc_dipole_moment, nuc_dipole_derivative,
     *         second_moment
      external return_1st_mom, return_2nd_mom 
      external energy_ty_denominator, reorder_energy

c 
c VFL SCF instructions 

      external compute_batch1
      external compute_batch2
      external compute_batch3
      external compute_batch4
      external compute_batch5
      external compute_batch6
      external compute_batch7
      external compute_batch8

      external compute_ubatch1
      external compute_ubatch2
      external compute_ubatch3
      external compute_ubatch4
      external compute_ubatch5
      external compute_ubatch6
      external compute_ubatch7
      external compute_ubatch8

      external form_iad
      external form_ibd
      external set_itol
      external fassign 
      external prefetch_on, prefetch_off
c
c -------------------------------------------------------------------- 
c VFL Instructions needed to perform 'super' T calculations 
c -------------------------------------------------------------------- 
c
      external set_t3blocks_a
      external set_t3blocks_i
      external compt3_a
      external compt3_i
c
c -------------------------------------------------------------------- 
c VFL Instruction needed to perform atomic density guess  
c -------------------------------------------------------------------- 
c
      external scf_atom 
c
      dummy = load_user_sub('sip_barrier' // char(0),
     *                       sip_barrier)
      dummy = load_user_sub('energy_denominator' // char(0),
     *                       energy_denominator)
      dummy = load_user_sub('energy_product' // char(0),
     *                       energy_product)
      dummy = load_user_sub('energy_reg_denominator' // char(0),
     *                       energy_reg_denominator)
      dummy = load_user_sub('energy_reg_product' // char(0),
     *                       energy_reg_product)
      dummy = load_user_sub('energy_denominator_reg_deriv' // char(0),
     *                       energy_denominator_reg_deriv)
      dummy = load_user_sub('print_scalar' // char(0),
     *                       print_scalar)
      dummy = load_user_sub('trace_on' // char(0),
     *                       trace_on)
      dummy = load_user_sub('trace_off' // char(0),
     *                       trace_off)
      dummy = load_user_sub('dump_block' // char(0),
     *                       dump_block)
      dummy = load_user_sub('load_balance_on' // char(0),
     *                       load_balance_on)
      dummy = load_user_sub('load_balance_off' // char(0),
     *                       load_balance_off)
      dummy = load_user_sub('blocks_to_list' // char(0),
     *                       blocks_to_list)
      dummy = load_user_sub('list_to_blocks' // char(0),
     *                       list_to_blocks)
      dummy = load_user_sub('fmult' // char(0),
     *                       fmult)
      dummy = load_user_sub('der_int_setup' // char(0),
     *                       der_int_setup)
      dummy = load_user_sub('compute_derivative_integrals'//char(0),
     *                       compute_derivative_integrals)
      dummy = load_user_sub('server_barrier'//char(0),
     *                       server_barrier)
      dummy = load_user_sub('scont1'//char(0), scont1)
      dummy = load_user_sub('hcont1'//char(0), hcont1)
      dummy = load_user_sub('dcont2'//char(0), dcont2)
      dummy = load_user_sub('eig'//char(0), eigen_calc)
      dummy = load_user_sub('eig_inv'//char(0), eig_inv)
      dummy = load_user_sub('eig_sr'//char(0), eig_sr)
      dummy = load_user_sub('eig_sr_inv'//char(0), eig_sr_inv)
      dummy = load_user_sub('copy_fock'//char(0), copy_fock)
      dummy = load_user_sub('diis_setup'//char(0), diis_setup)
      dummy = load_user_sub('compute_diis'//char(0), compute_diis)
      dummy = load_user_sub('write_blocks_to_list'//char(0), 
     *             write_blocks_to_list)
      dummy = load_user_sub('read_list_to_blocks'//char(0), 
     *             read_list_to_blocks)
      dummy = load_user_sub('remove_diagonal'//char(0), remove_diagonal)
      dummy = load_user_sub('return_diagonal'//char(0), return_diagonal)
      dummy = load_user_sub('fock_denominator'//char(0),
     *            fock_denominator)
      dummy = load_user_sub('set_flags'//char(0), set_flags)
      dummy = load_user_sub('set_flags2'//char(0), set_flags2)
      dummy = load_user_sub('der2_comp'//char(0), der2_comp)
      dummy = load_user_sub('fock_der'//char(0), fock_der)
      dummy = load_user_sub('overlap_der'//char(0), overlap_der)
      dummy = load_user_sub('scontxy'//char(0), scontxy)
      dummy = load_user_sub('hcontxy'//char(0), hcontxy)
      dummy = load_user_sub('compute_sderivative_integrals'//char(0),
     *                       compute_sderivative_integrals)
      dummy = load_user_sub('read_hess'//char(0), read_hess)
      dummy = load_user_sub('return_h1'//char(0), return_h1)

      dummy = load_user_sub('smon_on'//char(0), 
     *                       smon_on)
      dummy = load_user_sub('smon_off'//char(0), 
     *                       smon_off)
      dummy = load_user_sub('scf_rhf'//char(0), scf_rhf)
      dummy = load_user_sub('array_copy'//char(0), array_copy)

      dummy = load_user_sub('energy_adenominator'//char(0),
     *                       energy_adenominator) 
      dummy = load_user_sub('energy_bdenominator'//char(0),
     *                       energy_bdenominator) 
      dummy = load_user_sub('energy_abdenominator'//char(0),
     *                       energy_abdenominator) 
      dummy = load_user_sub('read_grad'//char(0), read_grad)
      dummy = load_user_sub('udenominator'//char(0), udenominator)
      dummy = load_user_sub('set_index'//char(0), set_index)
      dummy = load_user_sub('remove_single'//char(0), remove_single)
      dummy = load_user_sub('remove_double'//char(0), remove_double)
      dummy = load_user_sub('remove_single_double'//char(0), 
     *                       remove_single_double)
      dummy = load_user_sub('remove_ssss'//char(0), remove_ssss)
      dummy = load_user_sub('remove_dddd'//char(0), remove_dddd)
      dummy = load_user_sub('remove_xsxs'//char(0), remove_xsxs)
      dummy = load_user_sub('remove_xxss'//char(0), remove_xxss)
      dummy = load_user_sub('remove_xdxs'//char(0), remove_xdxs)
      dummy = load_user_sub('remove_xsxd'//char(0), remove_xsxd)
      dummy = load_user_sub('remove_xxsd'//char(0), remove_xxsd)
      dummy = load_user_sub('remove_xdxd'//char(0), remove_xdxd)
      dummy = load_user_sub('remove_xxdd'//char(0), remove_xxdd)
      dummy = load_user_sub('remove_xs'//char(0), remove_xs)
      dummy = load_user_sub('remove_xd'//char(0), remove_xd)
      dummy = load_user_sub('remove_ds'//char(0), remove_ds)
      dummy = load_user_sub('remove_ss'//char(0), remove_ss)
      dummy = load_user_sub('remove_dd'//char(0), remove_dd)
      dummy = load_user_sub('removeoo_dd'//char(0), removeoo_dd)
      dummy = load_user_sub('removevv_dd'//char(0), removevv_dd)
      dummy = load_user_sub('remove_sd'//char(0), remove_sd)
      dummy = load_user_sub('remove_sx'//char(0), remove_sx)
      dummy = load_user_sub('remove_dx'//char(0), remove_dx)
      dummy = load_user_sub('copy_ab'//char(0), copy_ab)
      dummy = load_user_sub('copy_ba'//char(0), copy_ba)
      dummy = load_user_sub('copy_ff'//char(0), copy_ff)
      dummy = load_user_sub('open_amp'//char(0), open_amp)
      dummy = load_user_sub('dump_amp'//char(0), dump_amp)
      dummy = load_user_sub('comp_ovl3c'//char(0), comp_ovl3c)
      dummy = load_user_sub('check_dconf'//char(0), check_dconf)
      dummy = load_user_sub('dropmo_expand_basis'//char(0), 
     *                       dropmo_expand_basis)
      dummy = load_user_sub('square_root'//char(0), square_root)
      dummy = load_user_sub('norm_fac'//char(0), norm_fac)
      dummy = load_user_sub('return_sval'//char(0), return_sval)
      dummy = load_user_sub('return_diagonal4'//char(0),
     *                      return_diagonal4)
      dummy = load_user_sub('symm_force_a'//char(0),
     *                      symm_force_a)
      dummy = load_user_sub('symm_force_i'//char(0),
     *                      symm_force_i)
      dummy = load_user_sub('place_sval'//char(0), place_sval)
      dummy = load_user_sub('place_one'//char(0), place_one)
      dummy = load_user_sub('place_one2'//char(0), place_one2)
      dummy = load_user_sub('place_one3'//char(0), place_one3)
      dummy = load_user_sub('place_one4'//char(0), place_one4)
      dummy = load_user_sub('place_one5'//char(0), place_one5)
      dummy = load_user_sub('place_one6'//char(0), place_one6)
      dummy = load_user_sub('place_oneb'//char(0), place_oneb)
      dummy = load_user_sub('place_one2b'//char(0), place_one2b)
      dummy = load_user_sub('place_one3b'//char(0), place_one3b)
      dummy = load_user_sub('place_one4b'//char(0), place_one4b)
      dummy = load_user_sub('place_one5b'//char(0), place_one5b)
      dummy = load_user_sub('place_one6b'//char(0), place_one6b)
      dummy = load_user_sub('place_one_dip'//char(0),
     *                                      place_one_dip)
      dummy = load_user_sub('place_one_dip_2'//char(0),
     *                                      place_one_dip_2)
      dummy = load_user_sub('place_one_dip_3'//char(0),
     *                                      place_one_dip_3)
      dummy = load_user_sub('place_one_dip_4'//char(0),
     *                                      place_one_dip_4)
      dummy = load_user_sub('place_one_dip_5'//char(0),
     *                                      place_one_dip_5)
      dummy = load_user_sub('place_one_dip_6'//char(0),
     *                                      place_one_dip_6)
      dummy = load_user_sub('place_one_dip_7'//char(0),
     *                                      place_one_dip_7)
      dummy = load_user_sub('place_one_dip_8'//char(0),
     *                                      place_one_dip_8)
      dummy = load_user_sub('place_one_dea'//char(0),
     *                                      place_one_dea)
      dummy = load_user_sub('place_one_dea_2'//char(0),
     *                                      place_one_dea_2)
      dummy = load_user_sub('place_one_dea_3'//char(0),
     *                                      place_one_dea_3)
      dummy = load_user_sub('place_one_dea_4'//char(0),
     *                                      place_one_dea_4)
      dummy = load_user_sub('place_one_dea_5'//char(0),
     *                                      place_one_dea_5)
      dummy = load_user_sub('place_one_dea_6'//char(0),
     *                                      place_one_dea_6)
      dummy = load_user_sub('place_one_dea_7'//char(0),
     *                                      place_one_dea_7)
      dummy = load_user_sub('place_one_dea_8'//char(0),
     *                                      place_one_dea_8)
      dummy = load_user_sub('eomroot_print'//char(0), eomroot_print)
      dummy = load_user_sub('eomroot_print_new'//char(0), 
     *                       eomroot_print_new)
      dummy = load_user_sub('smooth'//char(0), smooth)
      dummy = load_user_sub('smooth4'//char(0), smooth4)
      dummy = load_user_sub('eig_nonsymm'//char(0),
     *                       eigen_nonsymm_calc)
      dummy = load_user_sub('apply_den2'//char(0), apply_den2)
      dummy = load_user_sub('apply_den4'//char(0), apply_den4)
      dummy = load_user_sub('apply_den4_nodiag'//char(0),
     *                       apply_den4_nodiag)
      dummy = load_user_sub('energy_tdenominator'//char(0),
     *                       energy_tdenominator)
      
      dummy = load_user_sub('compute_aaaa_batch'//char(0),
     *                       compute_aaaa_batch)
      dummy = load_user_sub('compute_aaab_batch'//char(0),
     *                       compute_aaab_batch)
      dummy = load_user_sub('compute_aabb_batch'//char(0),
     *                       compute_aabb_batch)
      dummy = load_user_sub('compute_aabc_batch'//char(0),
     *                       compute_aabc_batch)
      dummy = load_user_sub('compute_abab_batch'//char(0),
     *                       compute_abab_batch)
      dummy = load_user_sub('compute_abac_batch'//char(0),
     *                       compute_abac_batch)
      dummy = load_user_sub('compute_abcd_batch'//char(0),
     *                       compute_abcd_batch)
      dummy = load_user_sub('compute_no4c_batch'//char(0),
     *                       compute_no4c_batch)
      dummy = load_user_sub('remove_atom_rud1'//char(0),
     *                       remove_atom_rud1)
      dummy = load_user_sub('remove_atom_rud2'//char(0),
     *                       remove_atom_rud2)
      dummy = load_user_sub('der4_comp'//char(0), der4_comp)
      dummy = load_user_sub('set_flags4'//char(0), set_flags4)
      dummy = load_user_sub('timestamp'//char(0), timestamp)
      dummy = load_user_sub('c1_print'//char(0), c1_print)
      dummy = load_user_sub('c1b_print'//char(0), c1b_print)
      dummy = load_user_sub('c2aa_print'//char(0), c2aa_print)
      dummy = load_user_sub('c2ab_print'//char(0), c2ab_print)
      dummy = load_user_sub('c2bb_print'//char(0), c2bb_print)
      dummy = load_user_sub('pardo_sects'//char(0), pardo_sects)
      dummy = load_user_sub('get_ijk'//char(0), get_ijk)
      dummy = load_user_sub('set_ijk_aaa'//char(0), set_ijk_aaa)
      dummy = load_user_sub('stripi'//char(0), stripi)
      dummy = load_user_sub('sum_64ss'//char(0), sum_64ss)
      dummy = load_user_sub('set_ijk_aab'//char(0), set_ijk_aab)
      dummy = load_user_sub('get_my_rank'//char(0), get_my_rank)
      dummy = load_user_sub('broadcast_array'//char(0), broadcast_array)
      dummy = load_user_sub('checkpoint'//char(0), checkpoint)
      dummy = load_user_sub('get_restart_status'//char(0), 
     *                       get_restart_status)
      dummy = load_user_sub('commit_checkpoint'//char(0), 
     *                       commit_checkpoint)
      dummy = load_user_sub('crash'//char(0), crash)
C
C
C             ...Watson: XYZ moment integrals
C
      dummy = load_user_sub('return_ovl'//char(0), return_ovl)
      dummy = load_user_sub('return_derv_xyz'//char(0), return_derv_xyz)
      dummy = load_user_sub('dipole_moment'//char(0), dipole_moment)
      dummy = load_user_sub('second_moment'//char(0), second_moment)
      dummy = load_user_sub('energy_ty_denominator'//char(0),
     *                       energy_ty_denominator)
      dummy = load_user_sub('reorder_energy'//char(0),reorder_energy)
      dummy = load_user_sub('return_1st_mom'//char(0), return_1st_mom)
      dummy = load_user_sub('return_2nd_mom'//char(0), return_2nd_mom)

      dummy = load_user_sub('nuc_dipole_moment'//char(0),
     *                       nuc_dipole_moment)
      dummy = load_user_sub('nuc_dipole_derivative'//char(0),
     *                       nuc_dipole_derivative)

      dummy = load_user_sub('return_x'//char(0), return_x)
      dummy = load_user_sub('return_y'//char(0), return_y)
      dummy = load_user_sub('return_z'//char(0), return_z)

      dummy = load_user_sub('return_xx'//char(0), return_xx)
      dummy = load_user_sub('return_yy'//char(0), return_yy)
      dummy = load_user_sub('return_zz'//char(0), return_zz)

      dummy = load_user_sub('return_xy'//char(0), return_xy)
      dummy = load_user_sub('return_xz'//char(0), return_xz)
      dummy = load_user_sub('return_yz'//char(0), return_yz)

      dummy = load_user_sub('print_eom_dens_info'//char(0),
     +                       print_eom_dens_info)
c
c -------------------------------------------------------------------- 
c VFL Instructions needed to perform 'fast' scf calculations 
c -------------------------------------------------------------------- 
c
      dummy = load_user_sub('compute_batch1'//char(0), compute_batch1)
      dummy = load_user_sub('compute_batch2'//char(0), compute_batch2)
      dummy = load_user_sub('compute_batch3'//char(0), compute_batch3)
      dummy = load_user_sub('compute_batch4'//char(0), compute_batch4)
      dummy = load_user_sub('compute_batch5'//char(0), compute_batch5)
      dummy = load_user_sub('compute_batch6'//char(0), compute_batch6)
      dummy = load_user_sub('compute_batch7'//char(0), compute_batch7)
      dummy = load_user_sub('compute_batch8'//char(0), compute_batch8)

      dummy = load_user_sub('compute_ubatch1'//char(0), compute_ubatch1)
      dummy = load_user_sub('compute_ubatch2'//char(0), compute_ubatch2)
      dummy = load_user_sub('compute_ubatch3'//char(0), compute_ubatch3)
      dummy = load_user_sub('compute_ubatch4'//char(0), compute_ubatch4)
      dummy = load_user_sub('compute_ubatch5'//char(0), compute_ubatch5)
      dummy = load_user_sub('compute_ubatch6'//char(0), compute_ubatch6)
      dummy = load_user_sub('compute_ubatch7'//char(0), compute_ubatch7)
      dummy = load_user_sub('compute_ubatch8'//char(0), compute_ubatch8)

      dummy = load_user_sub('form_iad'//char(0), form_iad)
      dummy = load_user_sub('form_ibd'//char(0), form_ibd)

      dummy = load_user_sub('set_itol'//char(0), set_itol)
      dummy = load_user_sub('fassign'//char(0), fassign)
c
c -------------------------------------------------------------------- 
c VFL Instructions needed to perform 'super' T calculations 
c -------------------------------------------------------------------- 
c
      dummy = load_user_sub('set_t3blocks_a'//char(0), set_t3blocks_a)
      dummy = load_user_sub('set_t3blocks_i'//char(0), set_t3blocks_i)
      dummy = load_user_sub('compt3_a'//char(0), compt3_a)
      dummy = load_user_sub('compt3_i'//char(0), compt3_i)
      dummy = load_user_sub('prefetch_on'//char(0), prefetch_on)
      dummy = load_user_sub('prefetch_off'//char(0), prefetch_off)
c
c -------------------------------------------------------------------- 
c VFL Instruction needed to perform atomic density guess  
c -------------------------------------------------------------------- 
c
      dummy = load_user_sub('scf_atom'//char(0), scf_atom) 
c
c -------------------------------------------------------------------- 


      return
      end
