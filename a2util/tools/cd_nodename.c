
#include <stdio.h>
#include <unistd.h> /* for chdir */
#include <sys/utsname.h> /* for utsname */

#include "f77_name.h"
#include "f_types.h"

f_int
F77_NAME(cd_nodename,CD_NODENAME)
()
{
    struct utsname UName;
    if (0>uname(&UName))
    {
        printf("@cd_nodename: unable to determine node name\n");
        return (f_int)1;
    }
    return (f_int)chdir(&UName.nodename[0]);
}

