#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifdef SIZE64BIT
#define RealNumber double
#else
#define RealNumber float
#endif
#define UrandomBuff 100000
int urandom_fd;
uint64_t *lp_urandom;

int get_random(uint64_t *data,int size)
{
    static int i=UrandomBuff;
    int j,k,ret_code;
    ret_code = 0;
    if(lp_urandom == NULL || urandom_fd < 0)
    {
        ret_code = -1;
        goto urandom_error;
    }
    
    k=0;
    while(k < size)
    {
        if(i == UrandomBuff)
        {
            j = read(urandom_fd,lp_urandom,UrandomBuff*sizeof(uint64_t));
            if(j != UrandomBuff*sizeof(uint64_t))
            {
                ret_code = -1;
                goto urandom_error;
            }
            i=0;
        }
        data[k++] = lp_urandom[i++];
    }

    urandom_error:
    return ret_code;
}

int get_sd(RealNumber *data,int size)
{
    static int i=2;
    static uint64_t n[2];
    static RealNumber f[2];
    RealNumber temp[2];
    int k,ret_code;
    ret_code = 0;
    
    k=0;
    while(k < size)
    {
        if(i == 2)
        {
            if(get_random(n,2))
            {
                ret_code = -1;
                goto error;
            }
            if(n[0] == 0)
            {
                n[0]++;
            }
            else if(n[0] == UINT64_MAX)
            {
                n[0]--;
            }
            if(n[1] == 0)
            {
                n[1]++;
            }
            else if(n[1] == UINT64_MAX)
            {
                n[1]--;
            }
            temp[0] = ((RealNumber)n[0])/((RealNumber)UINT64_MAX);
            temp[1] = ((RealNumber)n[1])/((RealNumber)UINT64_MAX);
            temp[0] = sqrt(-2.0 * log(temp[0]));
            temp[1] *= 2.0 * M_PI;
            if(!isfinite(temp[0]))
            {
                ret_code = -1;
                goto error;
            }
            if(!isfinite(temp[1]))
            {
                ret_code = -1;
                goto error;
            }
            f[0] = temp[0] * cos(temp[1]);
            f[1] = temp[0] * sin(temp[1]);
            i=0;
        }
        data[k++] = 1.0/(3.0 * f[i++]);
    }

    error:
    return ret_code;
}

int main(int argc,char **argv)
{
    int i,j,fd,ret_code;
    RealNumber *data;
    ret_code = 0;
    
    if(argc != 3)
    {
        ret_code = -1;
        goto error01;
    }
    i = strtol(argv[1],NULL,10);
    if(i <= 0)
    {
        ret_code = -1;
        goto error01;
    }
    
    urandom_fd = open("/dev/urandom",O_RDONLY);
    if(urandom_fd < 0)
    {
        ret_code = -1;
        goto error01;
    }
    lp_urandom = (uint64_t *)malloc(UrandomBuff*sizeof(uint64_t));
    if(lp_urandom == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    data = (RealNumber *)malloc(i*sizeof(RealNumber));
    if(data == NULL)
    {
        ret_code = -1;
        goto error03;
    }
    fd = open(argv[2],O_WRONLY|O_CREAT|O_TRUNC,0666);
    if(fd < 3)
    {
        ret_code = -1;
        goto error04;
    }
    get_sd(data,i);
    j=write(fd,data,i*sizeof(RealNumber));
    if(i*sizeof(RealNumber) != j)
    {
        ret_code = -1;
    }
    close(fd);
    error04:
    free(data);
    error03:
    free(lp_urandom);
    error02:
    close(urandom_fd);
    error01:
    return ret_code;
}
