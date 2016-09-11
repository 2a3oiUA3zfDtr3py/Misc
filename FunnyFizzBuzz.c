#include <stdio.h>
#include <stdlib.h>
int main(int argc,char **argv) {
    char a[] = "ACCFC@DHIAKENIKHGKJHLNGLIHOHDKFDHJCHEDC@OCB@DFMNGBJELNNLIBKJOKIJAJCB";
    int i,*j;
    i^=i;
    for(++i;i<sizeof(a) -!!i;i++) a[i] ^= i&((*a>>(!!i+!!i))-!!i);
    i^=i;
    for(++i;i<sizeof(a) -!!i;i++) a[i] -= *a;
    *a ^= *a;
    for(i^=i;(i+i)<sizeof(a);i++) a[i] = (a[i+i] << (!!i+!!i+!!i+!!i)) + a[i+i+!!i];
    ++*a;
    j = malloc(sizeof(int) * (*a + *(a + *a) + *(a + *a + *a)));
    
    i = (int) (*a) ^ (*a);
    j[i] = i++;
    j[i] = i >> i;
    j[i << i] = j[i];
    j[*a << *a | *a] = i;
    j[*(a + (*a))] = i ^ i;
    j[*(a + (*a)) + *a] = i << i;
    for(;i <= (int) *(a + (*a) + (*a)); i++) {
        j[i + *(a + (*a)) + *a] = i!=i;
        j[i + *(a + (*a)) + *a] |= j[i] & (*a << *a);
        j[i + *(a + (*a)) - *a] |= j[i] & (*a);
        j[i] &= *a << *a | *a;
        printf(a + a[j[i] + *(a + (*a)) - *a],i);
    }
    free(j);
    return i - i;
}
