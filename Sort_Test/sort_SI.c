#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#ifdef SIZE64BIT
#define Integer int64_t
#define UInteger uint64_t
#else
#define Integer int32_t
#define UInteger uint32_t
#endif
#ifndef Times
#define Times 1
#endif

#ifdef Merge
int sort(Integer *input_data,int size)
{
    int i,j,k,l,m,n,m_max,n_max,ret_code;
    Integer *data[2];
    ret_code = 0;
    data[0] = input_data;
    data[1] = (Integer *)malloc(size*sizeof(Integer));
    if(data[1] == NULL)
    {
        ret_code = -1;
        goto error01;
    }

    l = 0;
    for(i=1;i<size;i*=2)
    {
        for(j=0;j<size;j += i*2)
        {
            m = j;
            n = j + i;
            m_max = ((m + i) < size)? (m + i) : size;
            n_max = ((n + i) < size)? (n + i) : size;
            k = j;
            while(m < m_max && n < n_max)
            {
                if(data[l][m] < data[l][n])
                {
                    data[l^0x01][k++] = data[l][m++];
                }
                else
                {
                    data[l^0x01][k++] = data[l][n++];
                }
            }
            while(m < m_max)
            {
                data[l^0x01][k++] = data[l][m++];
            }
            while(n < n_max)
            {
                data[l^0x01][k++] = data[l][n++];
            }
        }
        l ^= 0x01;
    }

    if(l == 0x01)
    {
        for(i=0;i<size;i++)
        {
            data[0][i] = data[1][i];
        }
    }

    free(data[1]);
    error01:
    return ret_code;
}
#elif defined Radix
#define Alpha 0x100
#define Bravo 8
#define Charlie 1
int sort(Integer *input_data,int size)
{
    int i,j,k,ret_code;
    UInteger *data[2];
    int *v;
    ret_code = 0;
    data[0] = (UInteger *)input_data;
    data[1] = (UInteger *)malloc(size*sizeof(UInteger));
    if(data[1] == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    v = (int *)malloc((Alpha + 1)*sizeof(int));
    if(v == NULL)
    {
        ret_code = -1;
        goto error02;
    }

    k = 0;
    for(i=0;i<(sizeof(UInteger)/Charlie)-1;i++)
    {
        for(j=0;j<=Alpha;j++)
        {
            v[j] = 0;
        }
        for(j=0;j<size;j++)
        {
            v[((data[k][j] >> (i*Bravo)) & (Alpha - 1)) + 1]++;
        }
        for(j=1;j<=Alpha;j++)
        {
            v[j] += v[j-1];
        }
        for(j=0;j<size;j++)
        {
            data[k^0x01][v[((data[k][j] >> (i*Bravo)) & (Alpha - 1))]++] =  data[k][j];
        }
        k ^= 0x01;
    }
    for(j=0;j<=Alpha;j++)
    {
        v[j] = 0;
    }
    for(j=0;j<size;j++)
    {
        v[(((data[k][j] >> (i*Bravo)) & (Alpha - 1)) ^ (Alpha>>1)) + 1]++;
    }
    for(j=1;j<=Alpha;j++)
    {
        v[j] += v[j-1];
    }
    for(j=0;j<size;j++)
    {
        data[k^0x01][v[(((data[k][j] >> (i*Bravo)) & (Alpha - 1)) ^ (Alpha>>1))]++] =  data[k][j];
    }
    k ^= 0x01;
    
    if(k == 0x01)
    {
        for(i=0;i<size;i++)
        {
            data[0][i] = data[1][i];
        }
    }

    free(v);
    error02:
    free(data[1]);
    error01:
    return ret_code;
}
#elif defined Hybrid
#define Alpha 0x100
#define Bravo 8
#define Charlie 1
#define Delta 1048576
#define Echo 30
int sort(Integer *input_data,int size)
{
    int a,b,c,d,e,d_max,e_max,i,j,k,l,m,ret_code;
    UInteger *data[2],*td[2];
    int *v,*w;
    ret_code = 0;
    
    data[0] = (UInteger *)malloc(size*(sizeof(int) + sizeof(UInteger)));
    if(data[0] == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    if(sizeof(int) < sizeof(UInteger))
    {
        data[1] = data[0];
        w = (int *)(data[1] + size);
    }
    else
    {
        w = (int *)data[0];
        data[1] = (UInteger *)(w + size);
    }
    data[0] = (UInteger *)input_data;
    v = (int *)malloc((Alpha + 1) * sizeof(int));
    if(v == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    for(j=0;j<size;j++)
    {
        data[0][j] = data[0][j] ^ (1 << (BitSize - 1));
    }
    for(i=1;i<size;i++)
    {
        w[i] = 0;
    }
    w[0] = size;
    l = 0;
    k = sizeof(UInteger);
    while(k > 0)
    {
        k -= Charlie;
        if(k < 0)
        {
            k = 0;
        }
        i = 0;
        while(i < size)
        {
            m = w[i];
            if(m < 0)
            {
                i = -w[i];
                continue;
            }
            if((m - i) * sizeof(UInteger) < Delta)
            {
                m -= i;
                j = 0;
                td[0] = (data[l] + i);
                td[1] = (data[l^0x01] + i);
                e = (k + Charlie < sizeof(UInteger))? k + Charlie : sizeof(UInteger);
                if(Echo * e > m)
                {
                    for(a=1;a<m;a*=2)
                    {
                        for(b=0;b<m;b+=2*a)
                        {
                            c = b;
                            d = b;
                            e = b + a;
                            d_max = (d + a > m)? m : d + a;
                            e_max = (e + a > m)? m : e + a;
                            while(d < d_max && e < e_max)
                            {
                                if(td[j][d] < td[j][e])
                                {
                                    td[j^0x01][c++] = td[j][d++];
                                }
                                else
                                {
                                    td[j^0x01][c++] = td[j][e++];
                                }
                            }
                            while(d < d_max)
                            {
                                td[j^0x01][c++] = td[j][d++];
                            }
                            while(e < e_max)
                            {
                                td[j^0x01][c++] = td[j][e++];
                            }
                        }
                        j ^= 0x01;
                    }
                }
                else
                {
                    for(a=0;a<e;a++)
                    {
                        for(b=0;b<=Alpha;b++)
                        {
                            v[b] = 0;
                        }
                        for(b=0;b<m;b++)
                        {
                            v[((td[j][b] >> (a*Bravo)) & (Alpha - 1)) + 1]++;
                        }
                        for(b=1;b<=Alpha;b++)
                        {
                            v[b] += v[b-1];
                        }
                        for(b=0;b<m;b++)
                        {
                            td[j^0x01][v[(td[j][b] >> (a*Bravo)) & (Alpha - 1)]++] = td[j][b];
                        }
                        j ^= 0x01;
                    }
                }
                for(a=0;a<m;a++)
                {
                    td[j^0x01][a] = td[j][a];
                }
                w[i] = -(m + i);
                i += m;
                continue;
            }
            for(j=0;j<=Alpha;j++)
            {
                v[j] = 0;
            }
            for(j=i;j<m;j++)
            {
                v[((data[l][j] >> (k*Bravo)) & (Alpha - 1)) + 1]++;
            }
            v[0] = i;
            for(j=1;j<=Alpha;j++)
            {
                v[j] += v[j-1];
            }
            for(j=0;j<Alpha && v[j] < m;j++)
            {
                w[v[j]] = v[j+1];
            }
            for(j=i;j<m;j++)
            {
                data[l^0x01][v[(data[l][j] >> (k*Bravo)) & (Alpha - 1)]++] =  data[l][j];
            }
            i = m;
        }
        l ^= 0x01;
    }
    for(j=0;j<size;j++)
    {
        data[l][j] = data[l][j] ^ (1 << (BitSize - 1));
    }
    
    if(l == 0x01)
    {
        for(i=0;i<size;i++)
        {
            data[0][i] = data[1][i];
        }
    }
    
    free(v);
    error02:
    if(sizeof(int) < sizeof(UInteger))
    {
        free(data[1]);
    }
    else
    {
        free(w);
    }
    error01:
    return ret_code;
}
#else
int comp(const void *a,const void *b)
{
    if(*(Integer *)a > *(Integer *)b)
    {
        return 1;
    }
    else if(*(Integer *)a < *(Integer *)b)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}
int sort(Integer *input_data,int size)
{
    int ret_code;
    ret_code = 0;
    qsort(input_data,size,sizeof(Integer),comp);
    return ret_code;
}
#endif

int main()
{
    int i,j,fd,ret_code;
    off_t file_size,vel;
    Integer *data;
    struct stat file_stat;
    ret_code = 0;

    data = NULL;
    fd = open("temp.bin",O_RDONLY);
    if(fd < 3)
    {
        ret_code = -1;
        goto error01;
    }
    if(fstat(fd,&file_stat))
    {
        ret_code = -1;
        goto error02;
    }
    if(file_stat.st_size <= 0)
    {
        ret_code = -1;
        goto error02;
    }
    data = (Integer *)malloc(file_stat.st_size);
    if(data == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    file_size = read(fd,data,file_stat.st_size);
    if(file_size != file_stat.st_size || file_size % (Times*sizeof(Integer)) != 0)
    {
        ret_code = -1;
        goto error02;
    }
    error02:
    close(fd);
    if(ret_code)
    {
        goto error03;
    }

    i = file_size/(Times*sizeof(Integer));
    for(j=0;j<Times;j++)
    {
        sort(data + i*j,i);
    }
    printf("%d\n",i);

    
    fd = open("test.bin",O_WRONLY|O_CREAT|O_TRUNC,0666);
    if(fd < 3)
    {
        ret_code = -1;
        goto error03;
    }
    vel = write(fd,data,file_size);
    if(vel != file_size)
    {
        ret_code = -1;
    }
    close(fd);
    error03:
    if(data != NULL)
    {
        free(data);
    }
    error01:
    return ret_code;
}
