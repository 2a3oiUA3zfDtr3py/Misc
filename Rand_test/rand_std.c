#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i,*j,k,l;
    j = (int *)malloc(1000000*sizeof(int));
    if(j == NULL)
    {
        return 0;
    }
    srand(0);
    k = 0;
    for(l=0;l<10;l++)
    {
        for(i=0;i<1000000;i++)
        {
            j[i] = rand();
            if(j[i] < (RAND_MAX >> 1))
                k++;
        }
    }
    for(i=0;i<30;i++)
    {
        *(unsigned int *)(j+i) = (*(unsigned int *)(j+i)) << 1;
        printf("%d\n",j[i]);
    }
    printf("%d\n",k);    
    return 0;
}
