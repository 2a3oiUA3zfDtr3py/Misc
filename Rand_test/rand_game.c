#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
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
    i += (lag_a + 2) * sizeof(uint32_t);
    i += 0x100 * sizeof(uint8_t);
    rand_pool_a = (uint64_t *)malloc(i);
    if(rand_pool_a == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    rand_pool_b = (uint32_t *)(rand_pool_a + 1);
    s_box = (uint8_t *)(rand_pool_b + lag_a + 2);
    
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

    for(i=0;i<(lag_a + 2);i++)
    {
        apply_i_round(&seed);
        rand_pool_b[i] = seed;
    }
    
    rand_pool_a[0] %= 0xffffda35;
    error01:
    return ret_code;
}

void updata_event_key(int event_key)
{
    uint32_t *a;
    a = (uint32_t *)(&event_key);
    apply_sbox(rand_pool_b + lag_a + 1);
    rand_pool_b[lag_a+1] = (rand_pool_b[lag_a+1] << 4) | (rand_pool_b[lag_a+1] >> 28);
    rand_pool_b[lag_a+1] ^= *a;
    apply_sbox(rand_pool_b + lag_a + 1);
    rand_pool_b[lag_a+1] = (rand_pool_b[lag_a+1] << 4) | (rand_pool_b[lag_a+1] >> 28);
    apply_sbox(rand_pool_b + lag_a + 1);
    rand_pool_b[lag_a+1] = (rand_pool_b[lag_a+1] << 4) | (rand_pool_b[lag_a+1] >> 28);
    return;
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
        rand_pool_b[i] ^= rand_pool_b[lag_a+1];
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
    int i,j,k,ret_code;
    long l;
    uint32_t temp;
    char *str;
    const char *chosen[3] = {"グー","チョキ","パー"};
    ret_code = 0;
    
    if(My_init_rand(argc))
    {
        ret_code = -1;
        goto error01;
    }
    str = (char *)malloc(1024*sizeof(char));
    if(str == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    
    printf("0:グー\n1:チョキ\n2:パー\n");
    for(i=0;i<10;i++)
    {
        if(read(0,str,1023) <= 1)
        {
            break;
        }
        l = strtol(str,NULL,10);
        updata_event_key(l);
        if(l >= 3)
        {
            continue;
        }
        My_rand(&temp,1);
        temp %= 3;
        printf("Player=%s COM=%s\n",chosen[l],chosen[temp]);
        temp = (temp + 30 - l)%3;
        switch(temp)
        {
            case 0:
                printf("引き分け\n");
                break;
            case 1:
                printf("勝ち\n");
                break;
            case 2:
                printf("負け\n");
                break;
        }
    }
    
    free(str);
    error02:
    My_del_rand();
    error01:
    return ret_code;
}
