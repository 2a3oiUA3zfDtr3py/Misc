#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#define RealNumber double
#define UrandomBuff 100000
#define Dimension 8
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

int get_sud(RealNumber *data,int size)
{
    int64_t n;
    int i,ret_code;
    ret_code = 0;
    
    for(i=0;i<size;i++)
    {
        if(get_random((uint64_t *)&n,1))
        {
            ret_code = -1;
            goto error;
        }
        data[i] = ((RealNumber)n)/((RealNumber)INT64_MAX);
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
    i *= Dimension;
    
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
    get_sud(data,i);
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
