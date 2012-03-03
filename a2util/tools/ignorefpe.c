
/*
 * This routine turns off IEEE floating-point exception trapping
 * if there are no compiler flags to do it automatically (e.g., g77).
 */

#ifdef _TRAPFPE

#include <ieeefp.h>

void ignorefpe ()
{ fp_except_t i = fpsetmask(0); return; }

void ignorefpe_()
{ fp_except_t i = fpsetmask(0); return; }

void IGNOREFPE ()
{ fp_except_t i = fpsetmask(0); return; }

#endif /* _TRAPFPE */

