
#include <unistd.h> /* for chdir and access */
#include <stdio.h>  /* for snprintf */
#include <sys/utsname.h> /* for utsname */

#include "f77_name.h"
#include "f_types.h"

void
F77_NAME(cd_noderank,CD_NODERANK)
(f_int * iRank, f_int * iErr)
{
    char szSymLink[128];
    struct utsname UName;
    if (0>uname(&UName))
    {
        printf("@cd_noderank: unable to determine node name\n");
        *iErr = 1;
        return;
    }
    snprintf(&szSymLink[0],128,"%s.%i",UName.nodename,*iRank);
 /* switch to shared.<grank> if nodename.<grank> does not exist */
    if (access(&szSymLink[0],F_OK))
        snprintf(&szSymLink[0],128,"shared.%i",*iRank);
#ifdef _DEBUG
    printf("@cd_noderank: attempting to cd to %s\n",&szSymLink[0]);
#endif
    *iErr = chdir(&szSymLink[0]);
    return;
}

