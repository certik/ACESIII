
/*
 * This routine turns on IEEE floating-point exception trapping
 * if there are no compiler flags to do it automatically (e.g., g77).
 */

#ifdef _TRAPFPE

#include <ieeefp.h>
#include <signal.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

#include "f77_name.h"

void stupid_backtrace()
{
    char bt[256];
    char command[256];
    char prgname[256];
    char elink[256];
    pid_t pid = getpid();
    int lsize;
    FILE * gdb;

    sprintf(elink,"/proc/%i/file",pid);
    lsize = readlink(elink,prgname,99);
    if (lsize < 0)
    {
        fflush(NULL);
        fprintf(stderr,"\n@STUPID_BACKTRACE: readlink failed\n");
        exit(-1);
    }
    prgname[lsize] = 0;

    sprintf(command,
            "printf bt | gdb -n -q %s %i 2>/dev/null",prgname,pid);

    gdb = popen(command,"r");
    do
    {
        fgets(bt,200,gdb);
    } while ( ! (strncmp(bt,"#",1) == 0) );
    fprintf(stderr,"@STUPID_BACKTRACE: stack dump\n");
    do
    {
        fprintf(stderr,bt);
        fgets(bt,200,gdb);
    } while ( (strncmp(bt,"#",1) == 0) || (strncmp(bt,"  ",2) == 0) );
    pclose(gdb);
    fflush(stdout);
    fflush(stderr);
}

/******************************************************************************/

void
fpe_handler(int sig, siginfo_t * sip, void * uap)
{
    fp_except_t i; /* for resetting the mask */
    char *label;
    switch (sip->si_code)
    {
        case FPE_INTOVF: label = "integer overflow";                 break;
        case FPE_INTDIV: label = "integer divide by zero";           break;
        case FPE_FLTDIV: label = "floating point divide by zero";    break;
        case FPE_FLTOVF: label = "floating point overflow";          break;
        case FPE_FLTUND: label = "floating point underflow";         break;
        case FPE_FLTRES: label = "floating point inexact result";    break;
        case FPE_FLTINV: label = "invalid floating point operation"; break;
        case FPE_FLTSUB: label = "subscript out of range";           break;
        default:         label = "???";                              break;
    }
      /*
       * For your cutting and pasting pleasure!
         (sip->si_code == FPE_INTOVF) ||
         (sip->si_code == FPE_INTDIV) ||
         (sip->si_code == FPE_FLTINV) ||
         (sip->si_code == FPE_FLTSUB) ||
       */
    if ( (sip->si_code == FPE_FLTDIV) ||
         (sip->si_code == FPE_FLTOVF) ||
         (sip->si_code == FPE_FLTRES)    )
    {
        fflush(NULL);
        fprintf(stderr,
                "A floating-point exception has been caught.\n"
                "      type:    %s (0x%x)\n"
                "      address: %p\n"
                "      signal:  %d\n\n",
                label, sip->si_code, sip->si_addr, sig);
        /* stupid_backtrace(); */ /* This is broken right now. */
        exit(-1);
    }
    else
    {
        fpresetsticky( FP_X_INV | FP_X_DNML | FP_X_DZ | FP_X_OFL |
                       FP_X_UFL | FP_X_IMP  | FP_X_STK );
        /*
        i = fpsetmask( FP_X_INV | FP_X_OFL | FP_X_DZ );
        */
#define FPE_MASK FP_X_INV | FP_X_OFL | FP_X_DZ
        i = fpsetmask(FPE_MASK);
        return;
    }
}

/******************************************************************************/

void
F77_NAME(trapfpe,TRAPFPE)
()
{
    /*
     * set the FPE mask
     * FP_X_INV  : invalid operation
     * FP_X_DNML : denormal
     * FP_X_DZ   : zero divide
     * FP_X_OFL  : overflow
     * FP_X_UFL  : underflow <- BEWARE, this will kill almost any ACES prog!
     * FP_X_IMP  : (im)precision <- BEWARE, this will kill almost ANY  prog!
     * FP_X_STK  : stack fault
     */
    /*
    fp_except_t i = fpsetmask( FP_X_INV | FP_X_OFL | FP_X_DZ );
    */
    fp_except_t i = fpsetmask(FPE_MASK); /* This is reset in fpe_handler(). */

    /* change the FPE handler */
    struct sigaction sa_fpe;
    sa_fpe.sa_flags     = SA_SIGINFO; /* use sigaction instead of handler */
    sa_fpe.sa_sigaction = fpe_handler;
    sigemptyset(&sa_fpe.sa_mask);
    if (sigaction(SIGFPE,&sa_fpe,NULL) < 0)
        printf("@TRAPFPE: Floating-point exceptions cannot be caught.\n");
    return;
}

#endif /* _TRAPFPE */

