#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
uint64_t *rand_pool_a;
uint32_t *rand_pool_b;
uint8_t *s_box;
#define lag_a 127
#define lag_b 30
#define lag_c 16

void apply_sbox(uint32_t *a)
{
    uint8_t *temp;
    temp = (uint8_t *)a;
    temp[0] = s_box[temp[0]];
    temp[1] = s_box[temp[1]];
    temp[2] = s_box[temp[2]];
    temp[3] = s_box[temp[3]];
    return;
}

void apply_i_round(uint32_t *a)
{
    uint32_t b;
    apply_sbox(a);
    b = 1 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0x4743afb2;
    apply_sbox(a);
    b = 3 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0xde068085;
    apply_sbox(a);
    b = 4 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0x3cb8c837;
    apply_sbox(a);
    b = 2 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0x02f0047a;
    apply_sbox(a);
    b = 1 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0xec800aa9;
    apply_sbox(a);
    b = 1 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0x3e5c36dc;
    apply_sbox(a);
    b = 3 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0xb0d276d8;
    apply_sbox(a);
    b = 4 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0x3c58be19;
    apply_sbox(a);
    b = 2 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0x926cac14;
    apply_sbox(a);
    b = 1 + (*a & 0x03);
    *a = (*a << b) | (*a >> (32 - b));
    *a ^= 0xb46b7789;
    return;
}

void My_del_rand()
{
    if(rand_pool_a != NULL)
    {
        free(rand_pool_a);
    }
    rand_pool_a = NULL;
    rand_pool_b = NULL;
    s_box = NULL;
}

int My_init_rand(int seed_input)
{
    int i,j,k,ret_code;
    uint32_t seed,seed_temp;
    uint8_t sbox_temp;
    ret_code = 0;
    
    i = sizeof(uint64_t);
    i += (lag_a + 1) * sizeof(uint32_t);
    i += 0x100 * sizeof(uint8_t);
    rand_pool_a = (uint64_t *)malloc(i);
    if(rand_pool_a == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    rand_pool_b = (uint32_t *)(rand_pool_a + 1);
    s_box = (uint8_t *)(rand_pool_b + lag_a + 1);
    
    for(i=0;i<0x100;i++)
    {
        sbox_temp = i;
        s_box[i] = 0;
        for(k=0;k<5;k++)
        {
            s_box[i] ^= sbox_temp;
            sbox_temp = (sbox_temp << 1) | (sbox_temp >> 7);
        }
        s_box[i] ^= 0x63;
    }
    
    seed = *(unsigned int *)(&seed_input);
    apply_i_round(&seed);
    rand_pool_a[0] = seed;
    apply_i_round(&seed);
    rand_pool_a[0] |= ((uint64_t)seed) << 32;

    for(i=0;i<(lag_a + 1);i++)
    {
        apply_i_round(&seed);
        rand_pool_b[i] = seed;
    }
    
    rand_pool_a[0] %= 0xffffda35;
    
    error01:
    return ret_code;
}

void My_rand(uint32_t *a,int size)
{
    static int i=0;
    uint32_t b,c;
    int j,k,l,n;
    if(a == NULL)
    {
        return;
    }
    
    for(n=0;n<size;n++)
    {
        if(i>=lag_a) {
            i = 0;
            rand_pool_a[0] *= rand_pool_a[0];
            rand_pool_a[0] %= 0xffffda35;
            rand_pool_a[0] |= ((uint64_t)rand_pool_b[lag_c]) << 32;
            rand_pool_a[0] %= 0xffffda35;
            rand_pool_b[lag_a] = (rand_pool_b[lag_a] << 1) | (rand_pool_a[0] & 0x01);
            rand_pool_a[0] *= rand_pool_a[0];
            rand_pool_a[0] %= 0xffffda35;
            rand_pool_b[lag_a] = (rand_pool_b[lag_a] << 1) | (rand_pool_a[0] & 0x01);
            rand_pool_a[0] *= rand_pool_a[0];
            rand_pool_a[0] %= 0xffffda35;
            rand_pool_b[lag_a] = (rand_pool_b[lag_a] << 1) | (rand_pool_a[0] & 0x01);
            rand_pool_a[0] *= rand_pool_a[0];
            rand_pool_a[0] %= 0xffffda35;
            rand_pool_b[lag_a] = (rand_pool_b[lag_a] << 1) | (rand_pool_a[0] & 0x01);
        }
        j = i + lag_b;
        if(j>=lag_a) {
            j = 0;
        }
        k = i + 1;
        if(k>=lag_a)
        {
            k = 0;
        }
        apply_sbox(rand_pool_b + i);
        rand_pool_b[i] = (rand_pool_b[i] << 1) | (rand_pool_b[i] >> 31);
        rand_pool_b[i] ^= rand_pool_b[lag_a];
        apply_sbox(rand_pool_b + i);
        rand_pool_b[i] = (rand_pool_b[i] << 3) | (rand_pool_b[i] >> 29);
        rand_pool_b[i] ^= rand_pool_b[j];
        apply_sbox(rand_pool_b + i);
        rand_pool_b[i] = (rand_pool_b[i] << 4) | (rand_pool_b[i] >> 28);
        rand_pool_b[i] ^= rand_pool_b[k];
        
        a[n] = rand_pool_b[i];
        i++;
    }
    return;
}

int main(int argc,char **argv)
{
    int i,j,k,l;
    int32_t *temp;
    My_init_rand(argc);
    
    temp = (int32_t *)malloc(1000000*sizeof(int32_t));
    if(temp == NULL)
    {
        return 0;
    }
    
    k=0;
    l=0;
    for(i=0;i<10;i++)
    {
        My_rand((uint32_t *)temp,1000000);
        for(j=0;j<1000000;j+=2)
        {
            l++;
            if(temp[j]%2 != 0 && temp[j+1] %2 == 0)
                k++;
        }
    }
    My_rand((uint32_t *)temp,30);
    for(i=0;i<30;i++)
    {
        printf("%d\n",temp[i]);
    }
    printf("%d\n%d\n",k,l);
    
    free(temp);
    My_del_rand();
    return 0;
}
