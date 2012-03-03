
#ifndef _LISTNUM_H_
#define _LISTNUM_H_ /* ACES LIST HANDLES */

#define h_GS_IjKa_Ij_Ka		10
#define h_GS_IjAk_Ij_Ak		 9 /* only for UHF */

#define h_GS_iAjB_Bi_Aj		25

#define h_GS_aIbJ_bI_aJ		26 /* only for UHF */

#define h_GS_ABIJ_AI_BJ		 5
#define h_GS_abij_ai_bj		 6 /* only for UHF */
#define h_GS_AbIj_Ab_Ij		16
#define h_GS_AbIj_AI_bj		18
#define h_GS_AbIj_Aj_bI		21
#define h_GS_AbIj_bI_Aj		22 /* only for UHF */
#define h_GS_AbIj_bj_AI		17 /* only for UHF */

#define h_GS_AbCi_Ab_Ci		30
#define h_GS_AbIc_Ab_Ic		29 /* only for UHF */

#define h_GS_IjKl_Ij_Kl		13
#define h_HH_IjKl_Ij_Kl		13
#define h_HH_IJKL_IltJ_KltL	11 /* only for UHF */
#define h_HH_ijkl_iltj_kltl	12 /* only for UHF */

#define h_HH_IJKA_IltJ_KA	 7
#define h_HH_ijka_iltj_ka	 8 /* only for UHF */

#define h_HH_IAJB_BI_AJ		23
#define h_HH_iajb_bi_aj		24 /* only for UHF */

#define h_HH_ABIJ_AltB_IltJ	14
#define h_HH_ABIJ_AI_BJ		19
#define h_HH_abij_altb_iltj	15 /* only for UHF */
#define h_HH_abij_ai_bj		20 /* only for UHF */

#define h_HH_ABCI_AltB_CI	27 /* only for UHF */
#define h_HH_abci_altb_ci	28 /* only for UHF */

#define h_T1_IA_AI		90 /* "irrep" 1 */
#define h_T1_ia_ai		90 /* "irrep" 2 */ /* only for UHF */

#define h_T2_IJAB_AJ_BI		34
#define h_T2_IJAB_AltB_IltJ	44
#define h_T2_IjAb_AI_bj		37
#define h_T2_IjAb_Aj_bI		39
#define h_T2_IjAb_Ab_Ij		46
#define h_T2_IjAb_bI_Aj		38 /* only for UHF */
#define h_T2_IjAb_bj_AI		36 /* only for UHF */
#define h_T2_ijab_aj_bi		35 /* only for UHF */
#define h_T2_ijab_altb_iltj	45 /* only for UHF */

#ifdef _DO_STUFF_THAT_IS_UNCHECKED
#define h_D_IJAB_AltB_IltJ	48
#define h_D_IjAb_Ab_Ij		50
#define h_D_ijab_altb_iltj	49

#define h_d_IJAB_AltB_IltJ	64
#define h_d_IA_AI		64 /* only for "irrep" 9 */
#define h_d_ijab_altb_iltj	65
#define h_d_ia_ai		65 /* only for "irrep" 9 */
#define h_d_IjAb_Ab_Ij		66
#endif /* _DO_STUFF_THAT_IS_UNCHECKED */

#ifdef _DO_STUFF_THAT_IS_UNCHECKED
#define h_f_MI_MI		91 /* "irrep" 3 */
#define h_f_mi_mi		91 /* "irrep" 4 */

#define h_f_EA_EA		92 /* "irrep" 3 */
#define h_f_ea_ea		92 /* "irrep" 4 */

#define h_f_IA_AI		93 /* "irrep" 3 */
#define h_f_ia_ai		93 /* "irrep" 4 */
#endif /* _DO_STUFF_THAT_IS_UNCHECKED */

#endif /* _LISTNUM_H_ */

