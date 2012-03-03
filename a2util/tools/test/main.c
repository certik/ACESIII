
#include <string.h>

void main()
{
    double la[15]; long i=255; char a=0xff;
memset(&la,i,sizeof(la));for(i=0;i<15;i++)printf("%e ",la[i]);printf("\n");
memset(&la,a,sizeof(la));for(i=0;i<15;i++)printf("%e ",la[i]);printf("\n");
    return;
}

