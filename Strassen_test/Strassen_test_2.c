#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#define RealNumber float
#define UrandomBuff 100000
//matA_row is a data chanck
#define matA_col ((1024))
#define matA_row ((1024))
#define matX_col ((1024))
#define matX_row ((1024))
#define strassen_factor ((1))
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

int mul_mat(RealNumber **data_a,RealNumber **data_b,RealNumber **data_c,int x,int y,int z)
{
    int i,j,k;
    
    for(i=0;i<x;i++)
    {
        for(k=0;k<z;k++)
        {
            data_c[i][k] = 0.0;
        }
        for(j=0;j<y;j++)
        {
            for(k=0;k<z;k++)
            {
                data_c[i][k] += data_a[i][j] * data_b[j][k];
            }
        }
    }

    return 0;
}

int mul_mat2(RealNumber **data_a,RealNumber **data_x,RealNumber **data_b)
{
    int i,j,k,ret_code;
    ret_code = 0;

    for(i=0;i<matA_row;i++)
    {
        for(j=0;j<matX_col;j++)
        {
            data_b[i][j] = 0.0;
        }
        for(j=0;j<matA_col;j++)
        {
            for(k=0;k<matX_col;k++)
            {
                data_b[i][k] += data_a[i][j] * data_x[j][k];
            }
        }
    }
    
    return ret_code;
}

void strassen_mov_s2c(RealNumber **data_source,RealNumber **data_cache,int cache_row,int cache_col,int row,int col)
{
    int i,j;
    for(i=0;i<cache_row;i++)
    {
        for(j=0;j<cache_col;j++)
        {
            data_cache[i][j] = data_source[i + row][j + col];
        }
    }
    return;
}

void strassen_mov_c2s(RealNumber **data_cache,RealNumber **data_source,int cache_row,int cache_col,int row,int col)
{
    int i,j;
    for(i=0;i<cache_row;i++)
    {
        for(j=0;j<cache_col;j++)
        {
            data_source[i + row][j + col] = data_cache[i][j];
        }
    }
    return;
}

void strassen_mov_c2c(RealNumber **data_cache1,RealNumber **data_cache2,int cache_row,int cache_col)
{
    int i,j;
    for(i=0;i<cache_row;i++)
    {
        for(j=0;j<cache_col;j++)
        {
            data_cache2[i][j] = data_cache1[i][j];
        }
    }
    return;
}

void strassen_add_s2c(RealNumber **data_source,RealNumber **data_cache,int cache_row,int cache_col,int row,int col)
{
    int i,j;
    for(i=0;i<cache_row;i++)
    {
        for(j=0;j<cache_col;j++)
        {
            data_cache[i][j] += data_source[i + row][j + col];
        }
    }
    return;
}

void strassen_sub_s2c(RealNumber **data_source,RealNumber **data_cache,int cache_row,int cache_col,int row,int col)
{
    int i,j;
    for(i=0;i<cache_row;i++)
    {
        for(j=0;j<cache_col;j++)
        {
            data_cache[i][j] -= data_source[i + row][j + col];
        }
    }
    return;
}

void strassen_add_c2c(RealNumber **data_cache1,RealNumber **data_cache2,int cache_row,int cache_col)
{
    int i,j;
    for(i=0;i<cache_row;i++)
    {
        for(j=0;j<cache_col;j++)
        {
            data_cache2[i][j] += data_cache1[i][j];
        }
    }
    return;
}

void strassen_sub_c2c(RealNumber **data_cache1,RealNumber **data_cache2,int cache_row,int cache_col)
{
    int i,j;
    for(i=0;i<cache_row;i++)
    {
        for(j=0;j<cache_col;j++)
        {
            data_cache2[i][j] -= data_cache1[i][j];
        }
    }
    return;
}

void strassen_mul_mat(RealNumber **data_a,RealNumber **data_b,RealNumber **data_c,int x,int y,int z)
{
    int i,j,k;
    for(i=0;i<x;i++)
    {
        for(j=0;j<z;j++)
        {
            data_c[i][j] = 0.0;
        }
        for(j=0;j<y;j++)
        {
            for(k=0;k<z;k++)
            {
                data_c[i][k] += data_a[i][j] * data_b[j][k];
            }
        }
    }
    return;
}

int strassen(RealNumber **data_a,RealNumber **data_b,RealNumber **data_c,int x,int y,int z,RealNumber ***temp)
{
    int i,j,l,m,n,ret_code;
    int stack[strassen_factor];
    ret_code = 0;
    
    i = 1 << strassen_factor;
    if(x%i || y%i || z%i)
    {
        ret_code = -1;
        goto error;
    }

    temp[8] = data_a;
    temp[9] = data_b;
    temp[0] = data_c;
    stack[0] = 0;
    i=1;
    j=0;
    while(1)
    {
        l = x >> i;
        m = y >> i;
        n = z >> i;
        switch(j)
        {
            case 0:
                strassen_mov_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,0,0);
                strassen_add_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,l,m);
                strassen_mov_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,0,0);
                strassen_add_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,m,n);
                break;
            case 1:
                strassen_mov_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,l,0);
                strassen_add_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,l,m);
                strassen_mov_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,0,0);
                break;
            case 2:
                strassen_mov_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,0,0);
                strassen_mov_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,0,n);
                strassen_sub_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,m,n);
                break;
            case 3:
                strassen_mov_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,l,m);
                strassen_mov_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,m,0);
                strassen_sub_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,0,0);
                break;
            case 4:
                strassen_mov_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,0,0);
                strassen_add_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,0,m);
                strassen_mov_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,m,n);
                break;
            case 5:
                strassen_mov_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,l,0);
                strassen_sub_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,0,0);
                strassen_mov_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,0,0);
                strassen_add_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,0,n);
                break;
            case 6:
                strassen_mov_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,0,m);
                strassen_sub_s2c(temp[10*(i-1)+8],temp[10*i+8],l,m,l,m);
                strassen_mov_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,m,0);
                strassen_add_s2c(temp[10*(i-1)+9],temp[10*i+9],m,n,m,n);
                break;
            case 7:
                strassen_mov_c2c(temp[10*i+0],temp[10*i+7],l,n);
                strassen_add_c2c(temp[10*i+3],temp[10*i+7],l,n);
                strassen_sub_c2c(temp[10*i+4],temp[10*i+7],l,n);
                strassen_add_c2c(temp[10*i+6],temp[10*i+7],l,n);
                strassen_mov_c2s(temp[10*i+7],temp[10*(i-1)+stack[i-1]],l,n,0,0);
                strassen_mov_c2c(temp[10*i+2],temp[10*i+7],l,n);
                strassen_add_c2c(temp[10*i+4],temp[10*i+7],l,n);
                strassen_mov_c2s(temp[10*i+7],temp[10*(i-1)+stack[i-1]],l,n,0,n);
                strassen_mov_c2c(temp[10*i+1],temp[10*i+7],l,n);
                strassen_add_c2c(temp[10*i+3],temp[10*i+7],l,n);
                strassen_mov_c2s(temp[10*i+7],temp[10*(i-1)+stack[i-1]],l,n,l,0);
                strassen_mov_c2c(temp[10*i+0],temp[10*i+7],l,n);
                strassen_add_c2c(temp[10*i+2],temp[10*i+7],l,n);
                strassen_sub_c2c(temp[10*i+1],temp[10*i+7],l,n);
                strassen_add_c2c(temp[10*i+5],temp[10*i+7],l,n);
                strassen_mov_c2s(temp[10*i+7],temp[10*(i-1)+stack[i-1]],l,n,l,n);
                break;
        }

        if(j == 7 && i == 1)
        {
            break;
        }
        if(j == 7)
        {
            j = stack[--i];
            j++;
        }
        else if(i == strassen_factor && j < 7)
        {
            strassen_mul_mat(temp[10*i+8],temp[10*i+9],temp[10*i+j],l,m,n);
            j++;
        }
        else
        {
            stack[i++] = j;
            j=0;
            continue;
        }
    }

    error:
    return ret_code;
}

int main(int argc,char **argv)
{
    int i,j,k,l,m,ret_code;
    RealNumber **data_x,**data_a,**data_b,**data_c;
    RealNumber ***temp;
    clock_t time;
    ret_code = 0;

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
    
    data_x = (RealNumber **)malloc(matX_row*sizeof(RealNumber *));
    if(data_x == NULL)
    {
        ret_code = -1;
        goto error03;
    }
    data_x[0] = (RealNumber *)malloc(matX_row*matX_col*sizeof(RealNumber));
    if(data_x[0] == NULL)
    {
        ret_code = -1;
        goto error04;
    }
    for(i=1;i<matX_row;i++)
    {
        data_x[i] = data_x[i-1] + matX_col;
    }
    data_a = (RealNumber **)malloc(matA_row*sizeof(RealNumber *));
    if(data_a == NULL)
    {
        ret_code = -1;
        goto error05;
    }
    data_a[0] = (RealNumber *)malloc(matA_row*matA_col*sizeof(RealNumber));
    if(data_a[0] == NULL)
    {
        ret_code = -1;
        goto error06;
    }
    for(i=1;i<matA_row;i++)
    {
        data_a[i] = data_a[i-1] + matA_col;
    }
    data_b = (RealNumber **)malloc(matA_row*sizeof(RealNumber *));
    if(data_b == NULL)
    {
        ret_code = -1;
        goto error07;
    }
    data_b[0] = (RealNumber *)malloc(matA_row*matX_col*sizeof(RealNumber));
    if(data_b[0] == NULL)
    {
        ret_code = -1;
        goto error08;
    }
    for(i=1;i<matA_row;i++)
    {
        data_b[i] = data_b[i-1] + matX_col;
    }
    data_c = (RealNumber **)malloc(matA_row*sizeof(RealNumber *));
    if(data_c == NULL)
    {
        ret_code = -1;
        goto error09;
    }
    data_c[0] = (RealNumber *)malloc(matA_row*matX_col*sizeof(RealNumber));
    if(data_c[0] == NULL)
    {
        ret_code = -1;
        goto error10;
    }
    for(i=1;i<matA_row;i++)
    {
        data_c[i] = data_c[i-1] + matX_col;
    }

    temp = (RealNumber ***)malloc(10*(strassen_factor+1)*sizeof(RealNumber **));
    if(temp == NULL)
    {
        ret_code = -1;
        goto error11;
    }
    i = 0;
    for(j=0;j<strassen_factor;j++)
    {
        l = 2 << j;
        i += 9*(matA_row/l) + (matX_row/l);
    }
    temp[10] = (RealNumber **)malloc(i*sizeof(RealNumber *));
    if(temp[10] == NULL)
    {
        ret_code = -1;
        goto error12;
    }
    i = 0;
    for(j=0;j<strassen_factor;j++)
    {
        l = 2 << j;
        i += 8*(matA_row/l)*(matX_col/l);
        i += (matA_row/l)*(matA_col/l);
        i += (matX_row/l)*(matX_col/l);
    }
    temp[10][0] = (RealNumber *)malloc(i*sizeof(RealNumber));
    if(temp[10][0] == NULL)
    {
        ret_code = -1;
        goto error13;
    }
    i = 10;
    for(j=0;j<strassen_factor;j++)
    {
        l = 2 << j;
        for(k=0;k<9;k++)
        {
            temp[i+1] = temp[i] + (matA_row/l);
            i++;
        }
        if(j == strassen_factor -1)
        {
            break;
        }
        temp[i+1] = temp[i] + (matX_row/l);
        i++;
    }
    i = 0;
    for(j=0;j<strassen_factor;j++)
    {
        l = 2 << j;
        for(k=0;k<8*(matA_row/l);k++)
        {
            temp[10][i+1] = temp[10][i] + (matX_col/l);
            i++;
        }
        for(k=0;k<(matA_row/l);k++)
        {
            temp[10][i+1] = temp[10][i] + (matA_col/l);
            i++;
        }
        for(k=0;k<(matX_row/l)-1;k++)
        {
            temp[10][i+1] = temp[10][i] + (matX_col/l);
            i++;
        }
        
        if(j == strassen_factor-1)
        {
            break;
        }
        temp[10][i+1] = temp[10][i] + (matX_col/l);
        i++;
    }
    
    get_sud(data_a[0],matA_row*matA_col);
    get_sud(data_x[0],matX_row*matX_col);

    mul_mat(data_a,data_x,data_b,matA_row,matA_col,matX_col);
    strassen(data_a,data_x,data_c,matA_row,matA_col,matX_col,temp);
    time = clock();
    for(i=0;i<30;i++)
    {
        //mul_mat(data_a,data_x,data_b,matA_row,matA_col,matX_col);
        //mul_mat2(data_a,data_x,data_b);
        strassen(data_a,data_x,data_c,matA_row,matA_col,matX_col,temp);
    }
    time = clock() - time;


    double a,b;
    a = 0.0;
    for(i=0;i<matA_row*matX_col;i++)
    {
        b = data_b[0][i] - data_c[0][i];
        a += b*b;
    }
    
    printf("%f\n%f\n",a,((double)time)/((double)CLOCKS_PER_SEC));

    free(temp[10][0]);
    error13:
    free(temp[10]);
    error12:
    free(temp);
    error11:
    free(data_c[0]);
    error10:
    free(data_c);
    error09:
    free(data_b[0]);
    error08:
    free(data_b);
    error07:
    free(data_a[0]);
    error06:
    free(data_a);
    error05:
    free(data_x[0]);
    error04:
    free(data_x);
    error03:
    free(lp_urandom);
    error02:
    close(urandom_fd);
    error01:
    return ret_code;
}
