#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <MagickWand/MagickWand.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define num_parent 64
#define num_group_size 16
#define num_data_set 2
#define data_pool_size 2
#ifdef DOUBLE
#define RealNumber double
#define RealNumber_MAX DBL_MAX
#define RealNumber_MIN DBL_MIN
#else
#define RealNumber float
#define RealNumber_MAX FLT_MAX
#define RealNumber_MIN FLT_MIN
#endif
#define ThreadQuantity 2
#define UrandomBuff 100000
#define aux_vector_size 1
#define layer 12
#define input_datasize 3
#define output_datasize 3
#define hidden_channels 32
#define pix_ref 9
struct st_data
{
    RealNumber *input_image_data;
    RealNumber *output_image_data;
    int height;
    int width;
};

int matrix_size[layer][3] =
{
    {input_datasize,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,output_datasize,1}
};

int urandom_fd;
uint64_t *lp_urandom;
int element_size;
struct st_data *data;
int max_img_size;
int max_channel_size;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t rand_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

int vector_x_matrix(RealNumber *a,int a_size,RealNumber **b,int b_row,int b_col,RealNumber *c,int c_size)
{
    int i,j;
    for(j=0;j<b_col;j++)
    {
        c[j] = 0.0;
    }
    for(i=0;i<b_row;i++)
    {
        for(j=0;j<b_col;j++)
        {
            c[j] += b[i][j] * a[i];
        }
    }
    return 0;
}

int calc_regularization(RealNumber *gene,RealNumber *diff)
{
    int i,return_code;
    double t1,t2;
    RealNumber t3;
    return_code = 0;
    t1 = 0;/*
    for(i=0;i<element_size;i++)
    {
        t2 = (double)(gene[i]);
        if(t2 < 0.0)
        {
            t2 = -t2;
        }
        t2 = sqrt(t2);
        t1 += t2;
    }
    t1 *= t1;
    t1 *= 0.0000001;
    */

    t3 = (RealNumber)t1;
    if(isfinite(t3))
    {
        *diff = t3;
    }
    else
    {
        *diff = -1.0;
    }
    return return_code;
}

int calc(RealNumber *gene,int data_num,RealNumber *diff,RealNumber ***a,RealNumber **b,RealNumber *c)
{
    int return_code;
    int i,j,k,l,m,n,width,height;
    int p[pix_ref][2];
    double t1;
    RealNumber t2;
    return_code = 0;
    width = data[data_num].width;
    height = data[data_num].height;
    
    if(*diff < -0.5 || !isfinite(*diff))
    {
        *diff = -1.0;
        goto minus_diff;
    }

    a[0][0] = gene;
    for(j=1;j<pix_ref*matrix_size[0][0]+aux_vector_size;j++)
    {
        a[0][j] = a[0][j-1] + matrix_size[0][1];
    }
    for(i=1;i<layer;i++)
    {
        a[i][0] = a[i-1][pix_ref*matrix_size[i-1][0] + aux_vector_size -1] + matrix_size[i-1][1];
        for(j=1;j<pix_ref*matrix_size[i][0] + aux_vector_size;j++)
        {
            a[i][j] = a[i][j-1] + matrix_size[i][1];
        }
    }

    for(i=0;i<height*width*matrix_size[0][0];i++)
    {
        b[0][i] = data[data_num].input_image_data[i];
    }
    l=0;
    for(j=-1;j<=1;j++)
    {
        for(k=-1;k<=1;k++)
        {
            p[l][0] = j*matrix_size[0][2];
            p[l][1] = k*matrix_size[0][2];
            l++;
        }
    }
    for(j=0;j<height;j++)
    {
        for(k=0;k<width;k++)
        {
            for(n=0;n<pix_ref;n++)
            {
                if(j+p[n][0] >= 0 && j+p[n][0] < height && k+p[n][1] >= 0 && k+p[n][1] < width)
                {
                    for(m=0;m<matrix_size[0][0];m++)
                    {
                        c[matrix_size[0][0]*n + m] = b[0][matrix_size[0][0]*((j+p[n][0])*width+(k+p[n][1]))+m];
                    }
                }
                else
                {
                    for(m=0;m<matrix_size[0][0];m++)
                    {
                            c[matrix_size[0][0]*n + m] = -1.0;
                    }
                }
            }
            vector_x_matrix(c,pix_ref*matrix_size[0][0],a[0],pix_ref*matrix_size[0][0],matrix_size[0][1],
                            c+pix_ref*hidden_channels,matrix_size[0][1]);
            for(n=0;n<matrix_size[0][1];n++)
            {
                c[pix_ref*hidden_channels+n] += a[0][pix_ref*matrix_size[0][0]][n];
                b[1][matrix_size[0][1]*(j*width+k)+n] = c[pix_ref*hidden_channels+n];
            }
        }
    }

    for(i=1;i<layer;i++)
    {
        l=0;
        for(j=-1;j<=1;j++)
        {
            for(k=-1;k<=1;k++)
            {
                p[l][0] = j*matrix_size[i][2];
                p[l][1] = k*matrix_size[i][2];
                l++;
            }
        }
        for(j=0;j<height*width*matrix_size[i][0];j++)
        {
            if(b[i][j] < 0.0)
            {
                b[i][j] = 0.0;
            }
        }
        for(j=0;j<height;j++)
        {
            for(k=0;k<width;k++)
            {
                for(n=0;n<pix_ref;n++)
                {
                    if(j+p[n][0] < 0)
                    {
                        if(k+p[n][1] < 0)
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][m];
                            }
                        }
                        else if(k+p[n][1] < width)
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*(k+p[n][1])+m];
                            }
                        }
                        else
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*(width-1)+m];
                            }
                        }
                    }
                    else if(j+p[n][0] < height)
                    {
                        if(k+p[n][1] < 0)
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((j+p[n][0])*width)+m];
                            }
                        }
                        else if(k+p[n][1] < width)
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((j+p[n][0])*width+(k+p[n][1]))+m];
                            }
                        }
                        else
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((j+p[n][0])*width+(width-1))+m];
                            }
                        }
                    }
                    else
                    {
                        if(k+p[n][1] < 0)
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((height-1)*width)+m];
                            }
                        }
                        else if(k+p[n][1] < width)
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((height-1)*width+(k+p[n][1]))+m];
                            }
                        }
                        else
                        {
                            for(m=0;m<matrix_size[i][0];m++)
                            {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((height-1)*width+(width-1))+m];
                            }
                        }
                    }
                }
                vector_x_matrix(c,pix_ref*matrix_size[i][0],a[i],pix_ref*matrix_size[i][0],matrix_size[i][1],
                                c+pix_ref*hidden_channels,matrix_size[i][1]);
                for(n=0;n<matrix_size[i][1];n++)
                {
                    c[pix_ref*hidden_channels+n] += a[i][pix_ref*matrix_size[i][0]][n];
                    b[i+1][matrix_size[i][1]*(j*width+k)+n] = c[pix_ref*hidden_channels+n];
                }
            }
        }
    }
    
    t1 = 0.0;
    for(j=0;j<height*width*matrix_size[layer-1][1];j++)
    {
        b[layer][j] -= data[data_num].output_image_data[j];
        t1 += (double)(b[layer][j]*b[layer][j]);
    }
    t2 = (RealNumber)t1;
    if(isfinite(t2))
    {
        *diff += t2;
    }
    else
    {
        *diff = -1.0;
    }

    minus_diff:
    if(return_code)
    {
        *diff = -1.0;
    }
    return return_code;
}

int R_sort_p(RealNumber *score,int **index,int size)
{
    int i,j,k,l,m,n,m_max,n_max,ret_val;
    ret_val = 0;
    if(score == NULL || index == NULL || index[0] == NULL || index[1] == NULL)
    {
        ret_val = -1;
        goto error01;
    }

    for(i=0;i<size;i++)
    {
        index[0][i] = i;
    }
    l=0;
    for(i=1;i<size;i*=2)
    {
        for(j=0;j<size;j+=i*2)
        {
            k=j;
            m=j;
            n=j+i;
            m_max = (m+i<size)? m+i:size;
            n_max = (n+i<size)? n+i:size;
            while(m<m_max && n<n_max)
            {
                if(score[index[l][m]] >= RealNumber_MIN && score[index[l][n]] >= RealNumber_MIN)
                {
                    if(score[index[l][m]] < score[index[l][n]])
                    {
                        index[l^0x01][k++] = index[l][m++];
                    }
                    else
                    {
                        index[l^0x01][k++] = index[l][n++];
                    }
                }
                else if(score[index[l][m]] < -RealNumber_MIN && score[index[l][n]] < -RealNumber_MIN)
                {
                    index[l^0x01][k++] = index[l][m++];
                    index[l^0x01][k++] = index[l][n++];
                }
                else if(score[index[l][n]] < -RealNumber_MIN)
                {
                    index[l^0x01][k++] = index[l][m++];
                }
                else
                {
                    index[l^0x01][k++] = index[l][n++];
                }
            }
            while(m<m_max)
            {
                index[l^0x01][k++] = index[l][m++];
            }
            while(n<n_max)
            {
                index[l^0x01][k++] = index[l][n++];
            }
        }
        l ^= 0x01;
    }
    if(l == 1)
    {
        for(i=0;i<size;i++)
        {
            index[0][i] = index[1][i];
        }
    }
    error01:
    return ret_val;
}

int R_sort(RealNumber *score,int *idx,int size)
{
    int ret_val;
    int *index[2];
    ret_val = 0;
    if(idx == NULL || score == NULL)
    {
        ret_val = -1;
        goto error01;
    }
    index[0] = idx;
    index[1] = (int *)malloc(size*sizeof(int));
    if(index[1] == NULL)
    {
        ret_val = -1;
        goto error01;
    }
    if(R_sort_p(score,index,size))
    {
        ret_val = -1;
    }
    free(index[1]);
    error01:
    return ret_val;
}

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

int I_sort_p(uint64_t *r,int *idx_a,int *idx_b,int size)
{
    int i,j,k,l,m,n,m_max,n_max,ret_code;
    int *index[2];
    ret_code = 0;
    if(r == NULL || idx_a == NULL || idx_b == NULL || size <= 0)
    {
        ret_code = -1;
        goto error;
    }
    index[0] = idx_a;
    index[1] = idx_b;
    for(i=0;i<size;i++)
    {
        index[0][i] = i;
    }
    l=0;
    for(i=1;i<size;i*=2)
    {
        for(j=0;j<size;j+=i*2)
        {
            k=j;
            m=j;
            n=j+i;
            m_max = (m+i<size)? m+i:size;
            n_max = (n+i<size)? n+i:size;
            while(m<m_max && n<n_max)
            {
                if(r[index[l][m]] < r[index[l][n]])
                {
                    index[l^0x01][k++] = index[l][m++];
                }
                else
                {
                    index[l^0x01][k++] = index[l][n++];
                }
            }
            while(m<m_max)
            {
                index[l^0x01][k++] = index[l][m++];
            }
            while(n<n_max)
            {
                index[l^0x01][k++] = index[l][n++];
            }
        }
        l ^= 0x01;
    }
    if(l == 1)
    {
        for(i=0;i<size;i++)
        {
            index[0][i] = index[1][i];
        }
    }
    error:
    return ret_code;
}

int get_sud(double *data,int size)
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
        data[i] = ((double)n)/((double)INT64_MAX);
    }

    error:
    return ret_code;
}

int get_sd(double *data,int size)
{
    static int i=2;
    static uint64_t n[2];
    static double f[2];
    double temp[2];
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
            temp[0] = ((double)n[0])/((double)UINT64_MAX);
            temp[1] = ((double)n[1])/((double)UINT64_MAX);
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
        data[k++] = f[i++];
    }

    error:
    return ret_code;
}

int random_list(int *idx,int size)
{
    int i,ret_code;
    int *index;
    uint64_t *r;
    ret_code = 0;
    if(idx == NULL || size <= 0)
    {
        ret_code = -1;
        goto input_error;
    }
    index = (int *)malloc(size*sizeof(int));
    if(index == NULL)
    {
        ret_code = -1;
        goto malloc_error1;
    }
    r = (uint64_t *)malloc(size*sizeof(uint64_t));
    if(r == NULL)
    {
        ret_code = -1;
        goto malloc_error2;
    }
    i = get_random(r,size);
    if(i)
    {
        ret_code = -1;
        goto urandom_error;
    }
    
    I_sort_p(r,idx,index,size);
    
    urandom_error:
    free(r);
    malloc_error2:
    free(index);
    malloc_error1:
    input_error:
    return ret_code;
}

int save_gene(char *filename,RealNumber *gene)
{
    int fd,length,ret_code;
    ret_code = 0;
    fd = open(filename,O_WRONLY|O_CREAT|O_TRUNC,0666);
    if(fd < 0)
    {
        ret_code = -1;
        goto file_open_error;
    }
    length = write(fd,gene,element_size * sizeof(RealNumber));
    if(length != element_size * sizeof(RealNumber))
    {
        ret_code = -1;
    }
    
    close(fd);
    file_open_error:
    return ret_code;
}

int load_gene(char *filename,RealNumber *gene)
{
    int i,fd,ret_code;
    struct stat file_stat;
    ret_code = 0;
    fd = open(filename,O_RDONLY);
    if(fd < 0)
    {
        ret_code = -1;
        goto error01;
    }
    if(fstat(fd,&file_stat))
    {
        ret_code = -1;
        goto error02;
    }
    if(file_stat.st_size != element_size*sizeof(RealNumber))
    {
        ret_code = -1;
        goto error02;
    }
    i=read(fd,gene,file_stat.st_size);
    if(i != file_stat.st_size)
    {
        ret_code = -1;
        goto error02;
    }
    error02:
    close(fd);
    error01:
    return ret_code;
}

int get_filedata(struct st_data *i_data,int file_nunber)
{
    int i,ret_code = 0;
    unsigned long internal_width[2],internal_height[2];
    RealNumber *internal_data[2];
    MagickWand *img_handle;
    MagickBooleanType status;
    char file_name[256];
    char *file_base[2] = { "./data/%08da.png" , "./data/%08db.png" };

    if(file_nunber < 0 || file_nunber >= data_pool_size || i_data == NULL)
    {
        goto input_error;
        ret_code = -1;
    }

    if(i_data->input_image_data != NULL)
    {
        free(i_data->input_image_data);
        i_data->input_image_data = NULL;
    }
    if(i_data->output_image_data != NULL)
    {
        free(i_data->output_image_data);
        i_data->output_image_data = NULL;
    }

    MagickWandGenesis();
    img_handle = NewMagickWand();
    if(img_handle == NULL)
    {
        ret_code = -1;
        goto error1;
    }
    for(i=0;i<2;i++)
    {
        internal_data[i] = NULL;
    }
    for(i=0;i<2;i++)
    {
        sprintf(file_name,file_base[i],file_nunber);
        if(MagickReadImage(img_handle,file_name) == MagickFalse)
        {
            ret_code = -1;
            goto error2;
        }
        internal_width[i] = MagickGetImageWidth(img_handle);
        internal_height[i] = MagickGetImageHeight(img_handle);
        if(internal_width[i] == 0 || internal_height[i] == 0)
        {
            ret_code = -1;
            goto error2;
        }
        internal_data[i] = (RealNumber *)malloc(3*internal_width[i]*internal_height[i]*sizeof(RealNumber));
        if(internal_data[i] == NULL)
        {
            ret_code = -1;
            goto error2;
        }
#ifdef DOUBLE
        status = MagickExportImagePixels(img_handle,0,0,internal_width[i],internal_height[i],"RGB",DoublePixel,internal_data[i]);
#else
        status = MagickExportImagePixels(img_handle,0,0,internal_width[i],internal_height[i],"RGB",FloatPixel,internal_data[i]);
#endif
        if(status == MagickFalse)
        {
            ret_code = -1;
            goto error2;
        }
        ClearMagickWand(img_handle);
    }
    if(internal_width[0] != internal_width[1] || internal_height[0] != internal_height[1])
    {
        ret_code = -1;
        goto error2;
    }
    i_data->width = internal_width[0];
    i_data->height = internal_height[0];
    i_data->input_image_data = internal_data[0];
    internal_data[0] = NULL;
    i_data->output_image_data = internal_data[1];
    internal_data[1] = NULL;

    error2:
    if(ret_code)
    {
        for(i=0;i<2;i++)
        {
            if(internal_data[i] != NULL)
            {
                free(internal_data[i]);
                internal_data[i] = NULL;
            }
        }
    }
    DestroyMagickWand(img_handle);
    error1:
    MagickWandTerminus();
    input_error:
    return ret_code;
}

void* calc_p(void *arg)
{
    void **args = (void *)arg;
    RealNumber **gene  = (RealNumber **)args[0];
    RealNumber *diff  = (RealNumber  *)args[1];
    int length = *(int *)args[2];
    int *threadnum = (int *)args[3];
    int i,j;
    RealNumber ***a,**b,*c;
    
    a = (RealNumber ***)malloc(layer*sizeof(RealNumber **));
    if(a == NULL)
    {
        goto malloc_a_error1;
    }
    j=0;
    for(i=0;i<layer;i++)
    {
        j += pix_ref*matrix_size[i][0] + aux_vector_size;
    }
    a[0] = (RealNumber **)malloc(j*sizeof(RealNumber *));
    if(a[0] == NULL)
    {
        goto malloc_a_error2;
    }
    for(i=1;i<layer;i++)
    {
        a[i] = a[i-1] + pix_ref*matrix_size[i-1][0] + aux_vector_size;
    }
    b = (RealNumber **)malloc((layer+1)*sizeof(RealNumber *));
    if(b == NULL)
    {
        goto malloc_b_error1;
    }
    b[0] = (RealNumber *)malloc(2*max_channel_size*max_img_size*sizeof(RealNumber));
    if(b == NULL)
    {
        goto malloc_b_error2;
    }
    b[1] = b[0] + max_channel_size*max_img_size;
    for(i=2;i<=layer;i++)
    {
        b[i] = b[i&0x01];
    }
    c = (RealNumber *)malloc(2*pix_ref*max_channel_size*sizeof(RealNumber));
    if(c == NULL)
    {
        goto malloc_c_error;
    }

    for(j=0;j<length;j++)
    {
        calc_regularization(gene[j],diff + j);
    }
    for(i=0;i<num_data_set;i++)
    {
        for(j=0;j<length;j++)
        {
            calc(gene[j],i,diff + j,a,b,c);
        }
    }

    free(c);
    malloc_c_error:
    free(b[0]);
    malloc_b_error2:
    free(b);
    malloc_b_error1:
    free(a[0]);
    malloc_a_error2:
    free(a);
    malloc_a_error1:

    pthread_mutex_lock(&mutex);
    *threadnum = -1;
    pthread_cond_signal(&cond);
    pthread_mutex_unlock(&mutex);
    return NULL;
}

int superset_gene_calc(RealNumber **gene,RealNumber *diff) {
    int i,j,k,l,m,n,o,ret_code;
    int **threadnum,*threads_gene_number,**threads_gene_length;
    void ***lp;
    pthread_t *threads;
    RealNumber **internal_diff;
    ret_code = 0;
    lp = (void ***)malloc(ThreadQuantity*sizeof(void **));
    if(lp == NULL)
    {
        ret_code = -1;
        goto malloc_error1;
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        lp[i] = (void **)malloc(9*sizeof(void *));
        if(lp[i] == NULL)
        {
            ret_code = -1;
            goto malloc_error2;
        }
    }
    threads = (pthread_t *)malloc(ThreadQuantity*sizeof(pthread_t));
    if(threads == NULL)
    {
        ret_code = -1;
        goto malloc_error2;
    }
    internal_diff = (RealNumber **)malloc(ThreadQuantity*sizeof(RealNumber *));
    if(internal_diff == NULL)
    {
        ret_code = -1;
        goto malloc_error3;
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        internal_diff[i] = (RealNumber *)malloc(num_parent*sizeof(RealNumber));
        if(internal_diff[i] == NULL)
        {
            ret_code = -1;
            goto malloc_error4;
        }
    }
    threadnum = (int **)malloc(ThreadQuantity*sizeof(int *));
    if(threadnum == NULL)
    {
        ret_code = -1;
        goto malloc_error4;
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        threadnum[i] = (int *)malloc(sizeof(int));
        if(threadnum[i] == NULL)
        {
            ret_code = -1;
            goto malloc_error5;
        }
    }
    threads_gene_number = (int *)malloc(ThreadQuantity*sizeof(int));
    if(threads_gene_number == NULL)
    {
        ret_code = -1;
        goto malloc_error5;
    }
    threads_gene_length = (int **)malloc(ThreadQuantity*sizeof(int *));
    if(threads_gene_length == NULL)
    {
        ret_code = -1;
        goto malloc_error6;
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        threads_gene_length[i] = (int *)malloc(sizeof(int));
        if(threads_gene_length[i] == NULL)
        {
            ret_code = -1;
            goto malloc_error7;
        }
    }
    
    for(i=0;i<ThreadQuantity;i++)
    {
        lp[i][1] = (void *)internal_diff[i];
        lp[i][2] = (void *)threads_gene_length[i];
        lp[i][3] = (void *)threadnum[i];
        threadnum[i][0] = -2;
    }
    
    k = 0;
    while(1)
    {
        i = (num_parent - k)/(ThreadQuantity + 1);
        if(i == 0)
        {
            break;
        }
        for(l=0;l<ThreadQuantity;l++)
        {
            pthread_mutex_lock(&mutex);
            while(1)
            {
                m=0;
                for(n=0;n<ThreadQuantity;n++)
                {
                    m=threadnum[n][0];
                    if(m)
                    {
                        break;
                    }
                }
                if(m)
                {
                    break;
                }
                else
                {
                    pthread_cond_wait(&cond, &mutex);
                }
            }
            pthread_mutex_unlock(&mutex);
            if(m==-1)
            {
                pthread_join(threads[n],NULL);
                for(o=0;o<threads_gene_length[n][0];o++)
                {
                    diff[o+threads_gene_number[n]] = internal_diff[n][o];
                }
            }
            lp[n][0] = gene + k;
            threads_gene_number[n] = k;
            threads_gene_length[n][0] = i;
            threadnum[n][0] = 0;
            for(o=0;o<threads_gene_length[n][0];o++)
            {
                internal_diff[n][o] = diff[o+threads_gene_number[n]];
            }
            pthread_create(threads+n,NULL,(void *)calc_p,(void *)lp[n]);
            k += i;
        }
    }
    for(i=k;i<num_parent;i++)
    {
        pthread_mutex_lock(&mutex);
        while(1)
        {
            m=0;
            for(n=0;n<ThreadQuantity;n++)
            {
                m=threadnum[n][0];
                if(m)
                {
                    break;
                }
            }
            if(m)
            {
                break;
            }
            else
            {
                pthread_cond_wait(&cond, &mutex);
            }
        }
        pthread_mutex_unlock(&mutex);
        if(m == -1)
        {
            pthread_join(threads[n],NULL);
            pthread_mutex_lock(&mutex);
            threadnum[n][0] = -2;
            pthread_mutex_unlock(&mutex);
            for(o=0;o<threads_gene_length[n][0];o++)
            {
                diff[o+threads_gene_number[n]] = internal_diff[n][o];
            }
        }
        lp[n][0] = gene + i;
        threads_gene_number[n] = i;
        threads_gene_length[n][0] = 1;
        threadnum[n][0] = 0;
        for(o=0;o<threads_gene_length[n][0];o++)
        {
            internal_diff[n][o] = diff[o+threads_gene_number[n]];
        }
        pthread_create(threads+n,NULL,(void *)calc_p,(void *)lp[n]);
    }
    for(n=0;n<ThreadQuantity;n++)
    {
        pthread_mutex_lock(&mutex);
        m=threadnum[n][0];
        pthread_mutex_unlock(&mutex);
        if(m != -2)
        {
            pthread_join(threads[n],NULL);
            for(o=0;o<threads_gene_length[n][0];o++)
            {
                diff[o+threads_gene_number[n]] = internal_diff[n][o];
            }
        }
    }
    while(1)
    {
        pthread_mutex_lock(&mutex);
        j=-1;
        while(1)
        {
            m=0;
            for(n=0;n<ThreadQuantity;n++)
            {
                m=threadnum[n][0];
                if(m == 0)
                {
                    j=0;
                }
                else if(m == -1)
                {
                    break;
                }
            }
            if(m == -1  || (j!=0  && n == ThreadQuantity))
            {
                break;
            }
            else
            {
                pthread_cond_wait(&cond, &mutex);
            }
        }
        pthread_mutex_unlock(&mutex);
        if(m==-1)
        {
            pthread_join(threads[n],NULL);
            threadnum[n][0] = -2;
            for(o=0;o<threads_gene_length[n][0];o++)
            {
                diff[o+threads_gene_number[n]] = internal_diff[n][o];
            }
        }
        if(j && n == ThreadQuantity)
        {
            break;
        }
    }

    malloc_error7:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(threads_gene_length[i] != NULL)
        {
            free(threads_gene_length[i]);
        }
        else
        {
            break;
        }
    }
    free(threads_gene_length);
    malloc_error6:
    free(threads_gene_number);
    malloc_error5:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(threadnum[i] != NULL)
        {
            free(threadnum[i]);
        }
        else
        {
            break;
        }
    }
    free(threadnum);
    malloc_error4:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(internal_diff[i] != NULL)
        {
            free(internal_diff[i]);
        }
        else
        {
            break;
        }
    }
    free(internal_diff);
    malloc_error3:
    free(threads);
    malloc_error2:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(lp[i] != NULL)
        {
            free(lp[i]);
        }
        else
        {
            break;
        }
    }
    free(lp);
    malloc_error1:
    return ret_code;
}

void* calc_subset(void *arg)
{
    void **args = (void *)arg;
    RealNumber **gene  = (RealNumber **)args[0];
    RealNumber *diff = (RealNumber *)args[1];
    int *gene_index[3];
    int *threadnum = (int *)args[3];
    int i,j,l,q;
    uint64_t *urnd;
    RealNumber ***a,**b,*c,*calc_gene,temp_diff;
    double *temp_gene[2];
    gene_index[0]  = (int *)args[2];

    a = (RealNumber ***)malloc(layer*sizeof(RealNumber **));
    if(a == NULL)
    {
        goto malloc_a_error1;
    }
    j = 0;
    for(i=0;i<layer;i++)
    {
        j += pix_ref*matrix_size[i][0] + aux_vector_size;
    }
    a[0] = (RealNumber **)malloc(j*sizeof(RealNumber *));
    if(a[0] == NULL)
    {
        goto malloc_a_error2;
    }
    for(i=1;i<layer;i++)
    {
        a[i] = a[i-1] + pix_ref*matrix_size[i-1][0] + aux_vector_size;
    }
    b = (RealNumber **)malloc((layer+1)*sizeof(RealNumber *));
    if(b == NULL)
    {
        goto malloc_b_error1;
    }
    b[0] = (RealNumber *)malloc(2*max_img_size*max_channel_size*sizeof(RealNumber));
    if(b[0] == NULL)
    {
        goto malloc_b_error2;
    }
    b[1] = b[0] + max_img_size*max_channel_size;
    for(i=2;i<layer+1;i++)
    {
        b[i] = b[i&0x01];
    }
    c = (RealNumber *)malloc(2*pix_ref*max_channel_size*sizeof(RealNumber));
    if(c == NULL)
    {
        goto malloc_c_error;
    }
    urnd = (uint64_t *)malloc(num_group_size*(sizeof(uint64_t) + 2*sizeof(int)));
    if(urnd == NULL)
    {
        goto malloc_urnd_error;
    }
    gene_index[1] = (int *)(urnd + num_group_size);
    gene_index[2] = gene_index[1] + num_group_size;
    calc_gene = (RealNumber *)malloc(element_size*sizeof(RealNumber));
    if(calc_gene == NULL)
    {
        goto malloc_calc_gene_error;
    }
    temp_gene[0] = (double *)malloc(2*(2+element_size)*sizeof(double));
    if(temp_gene[0] == NULL)
    {
        goto malloc_temp_gene_error;
    }
    temp_gene[1] = temp_gene[0] + 2+element_size;

    q = (num_group_size / 2) - 1;
    if(q < 2)
    {
        q = 2;
    }

    for(l=0;l<num_group_size;l++)
    {
        R_sort_p(diff,gene_index,num_group_size);
        gene_index[1][0] = 0;
        gene_index[2][0] = num_group_size;
        for(i=1;i<num_group_size;i++)
        {
            gene_index[1][i] = i;
            gene_index[2][i] = gene_index[2][i-1] + num_group_size -i;
        }
        pthread_mutex_lock(&rand_mutex);
        get_random(urnd,q);
        pthread_mutex_unlock(&rand_mutex);
        for(i=0;i<q;i++)
        {
            urnd[i] %= gene_index[2][num_group_size-1-i];
            for(j=0;j<num_group_size-1-i;j++)
            {
                if(urnd[i] < gene_index[2][j])
                {
                    break;
                }
            }
            urnd[i] = gene_index[1][j];
            for(;j<num_group_size-1-i;j++)
            {
                gene_index[1][j] = gene_index[1][j+1];
                gene_index[2][j] = gene_index[2][j+1] -urnd[i];
            }
        }
        I_sort_p(urnd,gene_index[2],gene_index[1],q);
        
        for(i=0;i<q;i++)
        {
            gene_index[1][i] = urnd[gene_index[2][i]];
        }
        for(i=0;i<element_size+2;i++)
        {
            temp_gene[0][i] = 0.0;
        }
        for(i=0;i<q-1;i++)
        {
            temp_gene[0][element_size+1] = 2.0 *((double)(q-i)/(double)q) - 1.0;
            temp_gene[0][element_size] += temp_gene[0][element_size+1];
            for(j=0;j<element_size;j++)
            {
                temp_gene[0][j] += temp_gene[0][element_size+1] * ((double)(gene[gene_index[1][i]][j]));
            }
        }
        for(i=0;i<element_size;i++)
        {
            temp_gene[0][i] /= temp_gene[0][element_size];
            calc_gene[i] = (RealNumber)(2.0*temp_gene[0][i] - ((double)(gene[gene_index[1][q-1]][i])));
        }
        calc_regularization(calc_gene,&temp_diff);
        for(j=0;j<num_data_set;j++)
        {
            calc(calc_gene,j,&temp_diff,a,b,c);
        }
        j=0;
        if(diff[gene_index[1][q-1]] > -0.0 && temp_diff > -0.0 && diff[gene_index[1][q-1]] > temp_diff)
        {
            j=1;
        }
        if(diff[gene_index[1][q-1]] < -0.0 && temp_diff > -0.0)
        {
            j=1;
        }
        if(j == 1)
        {
            for(i=0;i<element_size;i++)
            {
                gene[gene_index[1][q-1]][i] = calc_gene[i];
            }
            diff[gene_index[1][q-1]] = temp_diff;
            continue;
        }
        for(i=0;i<element_size;i++)
        {
            calc_gene[i] = (RealNumber)(0.5*(temp_gene[0][i] + ((double)(gene[gene_index[1][q-1]][i])) ));
        }
        calc_regularization(calc_gene,&temp_diff);
        for(j=0;j<num_data_set;j++)
        {
            calc(calc_gene,j,&temp_diff,a,b,c);
        }
        j=0;
        if(diff[gene_index[1][q-1]] > -0.0 && temp_diff > -0.0 && diff[gene_index[1][q-1]] > temp_diff)
        {
            j=1;
        }
        if(diff[gene_index[1][q-1]] < -0.0 && temp_diff > -0.0)
        {
            j=1;
        }
        if(j == 1)
        {
            for(i=0;i<element_size;i++)
            {
                gene[gene_index[1][q-1]][i] = calc_gene[i];
            }
            diff[gene_index[1][q-1]] = temp_diff;
            continue;
        }
        pthread_mutex_lock(&rand_mutex);
        get_sud(temp_gene[1],2+element_size);
        pthread_mutex_unlock(&rand_mutex);
        for(i=0;i<element_size;i++)
        {
            temp_gene[1][i] *= 0.1;
            temp_gene[1][i] += temp_gene[0][i];
            gene[gene_index[1][q-1]][i] = (RealNumber)temp_gene[1][i];
        }
        calc_regularization(gene[gene_index[1][q-1]],&temp_diff);
        for(j=0;j<num_data_set;j++)
        {
            calc(gene[gene_index[1][q-1]],j,&temp_diff,a,b,c);
        }
        diff[gene_index[1][q-1]] = temp_diff;
    }
    R_sort_p(diff,gene_index,num_group_size);

    free(temp_gene[0]);
    malloc_temp_gene_error:
    free(calc_gene);
    malloc_calc_gene_error:
    free(urnd);
    malloc_urnd_error:
    free(c);
    malloc_c_error:
    free(b[0]);
    malloc_b_error2:
    free(b);
    malloc_b_error1:
    free(a[0]);
    malloc_a_error2:
    free(a);
    malloc_a_error1:

    pthread_mutex_lock(&mutex);
    *threadnum = -1;
    pthread_cond_signal(&cond);
    pthread_mutex_unlock(&mutex);
    return NULL;
}

int subset_gene_calc(RealNumber **superset_gene,RealNumber *superset_diff,int *superset_idx)
{
    int i,j,k,m,n,ret_code,grupe_num;
    RealNumber ***subset_gene,**internal_diff;
    int *threadnum,*thread_grupe_num,**subset_idx;
    void ***lp;
    pthread_t *threads;
    ret_code = 0;

    internal_diff = (RealNumber **)malloc(ThreadQuantity*sizeof(RealNumber *));
    if(internal_diff == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        internal_diff[i] = (RealNumber *)malloc(num_parent*sizeof(RealNumber));
        if(internal_diff[i] == NULL)
        {
            ret_code = -1;
            goto error01;
        }
    }
    subset_gene = (RealNumber ***)malloc(ThreadQuantity*sizeof(RealNumber **));
    if(subset_gene == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        subset_gene[i] = (RealNumber **)malloc(num_group_size*sizeof(RealNumber *));
        if(subset_gene[i] == NULL)
        {
            ret_code = -1;
            goto error03;
        }
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        subset_gene[i][0] = (RealNumber *)malloc(element_size*num_group_size*sizeof(RealNumber));
        if(subset_gene[i][0] == NULL)
        {
            ret_code = -1;
            goto error04;
        }
        for(j=1;j<num_group_size;j++)
        {
            subset_gene[i][j] = subset_gene[i][j-1] + element_size;
        }
    }
    subset_idx = (int **)malloc(ThreadQuantity*sizeof(int *));
    if(subset_idx == NULL)
    {
        ret_code = -1;
        goto error04;
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        subset_idx[i] = (int *)malloc(num_group_size*sizeof(int));
        if(subset_idx[i] == NULL)
        {
            ret_code = -1;
            goto error05;
        }
        for(j=0;j<num_group_size;j++)
        {
            subset_idx[i][j] = j;
        }
    }
    thread_grupe_num = (int *)malloc(ThreadQuantity*sizeof(int));
    if(thread_grupe_num == NULL)
    {
        ret_code = -1;
        goto error05;
    }
    threadnum = (int *)malloc(ThreadQuantity*sizeof(int));
    if(threadnum == NULL)
    {
        ret_code = -1;
        goto error06;
    }
    threads = (pthread_t *)malloc(ThreadQuantity*sizeof(pthread_t));
    if(threads == NULL)
    {
        ret_code = -1;
        goto error07;
    }
    lp = (void ***)malloc(ThreadQuantity*sizeof(void **));
    if(lp == NULL)
    {
        ret_code = -1;
        goto error08;
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        lp[i] = (void **)malloc(9*sizeof(void *));
        if(lp[i] == NULL)
        {
            ret_code = -1;
            goto error09;
        }
    }
    for(i=0;i<ThreadQuantity;i++)
    {
        lp[i][0] = (void *)subset_gene[i];
        lp[i][1] = (void *)internal_diff[i];
        lp[i][2] = (void *)subset_idx[i];
        lp[i][3] = (void *)(threadnum + i);
        threadnum[i] = -2;
    }

    grupe_num = num_parent / num_group_size;
    for(i=0;i<grupe_num;i++)
    {
        pthread_mutex_lock(&mutex);
        while(1)
        {
            m=0;
            for(n=0;n<ThreadQuantity;n++)
            {
                m=threadnum[n];
                if(m != 0)
                {
                    break;
                }
            }
            if(m != 0)
            {
                break;
            }
            else
            {
                pthread_cond_wait(&cond, &mutex);
            }
        }
        pthread_mutex_unlock(&mutex);
        if(m==-1)
        {
            pthread_join(threads[n],NULL);
            for(j=0;j<num_group_size;j++)
            {
                for(k=0;k<element_size;k++)
                {
                    superset_gene[superset_idx[grupe_num*j + thread_grupe_num[n]]][k] = subset_gene[n][subset_idx[n][j]][k];
                }
                superset_diff[superset_idx[grupe_num*j + thread_grupe_num[n]]] = internal_diff[n][subset_idx[n][j]];
            }
        }
        for(j=0;j<num_group_size;j++)
        {
            for(k=0;k<element_size;k++)
            {
                subset_gene[n][j][k] = superset_gene[superset_idx[grupe_num*j + i]][k];
            }
            internal_diff[n][j] = superset_diff[superset_idx[grupe_num*j + i]];
        }
        thread_grupe_num[n] = i;
        pthread_mutex_lock(&mutex);
        threadnum[n] = 0;
        pthread_mutex_unlock(&mutex);
        pthread_create(threads+n,NULL,(void *)calc_subset,(void *)lp[n]);
    }
    while(1)
    {
        pthread_mutex_lock(&mutex);
        while(1)
        {
            m=0;
            j=-1;
            for(n=0;n<ThreadQuantity;n++)
            {
                m=threadnum[n];
                if(m == 0)
                {
                    j=0;
                }
                else if(m == -1)
                {
                    break;
                }
            }
            if(m == -1  || (j!=0  && n == ThreadQuantity))
            {
                break;
            }
            else
            {
                pthread_cond_wait(&cond, &mutex);
            }
        }
        pthread_mutex_unlock(&mutex);
        if(m == -1)
        {
            pthread_join(threads[n],NULL);
            for(j=0;j<num_group_size;j++)
            {
                for(k=0;k<element_size;k++)
                {
                    superset_gene[superset_idx[grupe_num*j + thread_grupe_num[n]]][k] = subset_gene[n][subset_idx[n][j]][k];
                }
                superset_diff[superset_idx[grupe_num*j + thread_grupe_num[n]]] = internal_diff[n][subset_idx[n][j]];
            }
            pthread_mutex_lock(&mutex);
            threadnum[n] = -2;
            pthread_mutex_unlock(&mutex);
        }
        if(j!=0  && n == ThreadQuantity)
        {
            break;
        }
    }

    error09:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(lp[i] != NULL)
        {
            free(lp[i]);
        }
        else
        {
            break;
        }
    }
    free(lp);
    error08:
    free(threads);
    error07:
    free(threadnum);
    error06:
    free(thread_grupe_num);
    error05:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(subset_idx[i] != NULL)
        {
            free(subset_idx[i]);
        }
        else
        {
            break;
        }
    }
    free(subset_idx);
    error04:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(subset_gene[i][0] != NULL)
        {
            free(subset_gene[i][0]);
        }
        else
        {
            break;
        }
    }
    error03:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(subset_gene[i] != NULL)
        {
            free(subset_gene[i]);
        }
        else
        {
            break;
        }
    }
    free(subset_gene);
    error02:
    for(i=0;i<ThreadQuantity;i++)
    {
        if(internal_diff[i] != NULL)
        {
            free(internal_diff[i]);
        }
        else
        {
            break;
        }
    }
    free(internal_diff);
    error01:
    return ret_code;
}

int main()
{
    int i,j,ret_code,*superset_idx,*file_index;
    RealNumber **superset_gene,*superset_diff;
    double data_diff[3];
    char *file_name;
    ret_code = 0;
    element_size = 0;
    max_channel_size = matrix_size[0][0];
    for(i=0;i<layer;i++)
    {
        element_size += (pix_ref*matrix_size[i][0] + aux_vector_size) * matrix_size[i][1];
        if(max_channel_size < matrix_size[i][1])
        {
            max_channel_size = matrix_size[i][1];
        }
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

    data = (struct st_data *)malloc(num_data_set*sizeof(struct st_data));
    if(data == NULL)
    {
        ret_code = -1;
        goto error03;
    }
    for(i=0;i<num_data_set;i++)
    {
        data[i].input_image_data = NULL;
        data[i].output_image_data = NULL;
        data[i].width = 0;
        data[i].height = 0;
    }
    superset_gene = (RealNumber **)malloc(num_parent*sizeof(RealNumber *));
    if(superset_gene == NULL)
    {
        ret_code = -1;
        goto error04;
    }
    superset_gene[0] = (RealNumber *)malloc(element_size*num_parent*sizeof(RealNumber));
    if(superset_gene[0] == NULL)
    {
        ret_code = -1;
        goto error05;
    }
    for(i=1;i<num_parent;i++)
    {
        superset_gene[i] = superset_gene[i-1] + element_size;
    }
    superset_diff = (RealNumber *)malloc(num_parent*sizeof(RealNumber));
    if(superset_diff == NULL)
    {
        ret_code = -1;
        goto error06;
    }
    superset_idx = (int *)malloc((num_parent + data_pool_size)*sizeof(int));
    if(superset_idx == NULL)
    {
        ret_code = -1;
        goto error07;
    }
    file_index = superset_idx + num_parent;
    file_name = (char *)malloc(1024*sizeof(char));
    if(file_name == NULL)
    {
        ret_code = -1;
        goto error08;
    }
   
    for(i=0;i<num_parent;i++)
    {
        sprintf(file_name,"./gene/gene%08d.bin",i);
        if(load_gene(file_name,superset_gene[i]))
        {
            ret_code = -1;
            goto error09;
        }
    }

    random_list(file_index,data_pool_size);
    data_diff[0] = 0.0;
    max_img_size = 0;
    for(i=0;i<num_data_set;i++)
    {
        if(get_filedata(data + i,file_index[i]))
        {
            ret_code = -1;
            goto error10;
        }
        if(max_img_size < data[i].height * data[i].width)
        {
            max_img_size = data[i].height * data[i].width;
        }
        data_diff[1] = 0.0;
        for(j=0;j < data[i].height * data[i].width * matrix_size[0][0];j++)
        {
            data_diff[2] = ((double)(data[i].input_image_data[j])) - ((double)(data[i].output_image_data[j]));
            data_diff[1] += data_diff[2] * data_diff[2];
        }
        data_diff[0] += data_diff[1];
    }
    for(i=0;i<num_parent;i++)
    {
        superset_diff[i] = 0.0;
    }
    superset_gene_calc(superset_gene,superset_diff);
    R_sort(superset_diff,superset_idx,num_parent);

    fcntl(0,F_SETFL,O_NONBLOCK);
    while(1)
    {
        i=0;
        while(i<num_parent)
        {
            printf("%f\n",superset_diff[superset_idx[num_parent-i-1]]/((RealNumber)num_data_set));
            i++;
        }
        printf("%d\n%f\n",superset_idx[0],data_diff[0]/((RealNumber)num_data_set));
        subset_gene_calc(superset_gene,superset_diff,superset_idx);
        R_sort(superset_diff,superset_idx,num_parent);
        for(i=0;i<num_parent;i++)
        {
            sprintf(file_name,"./gene/gene%08d.bin",i);
            if(save_gene(file_name,superset_gene[superset_idx[i]]))
            {
                ret_code = -1;
                goto error10;
            }
        }
        if(read(0,file_name,1023)>=0)
        {
            break;
        }
    }

    error10:
    for(i=0;i<num_data_set;i++)
    {
        if(data[i].input_image_data != NULL)
        {
            free(data[i].input_image_data);
            data[i].input_image_data = NULL;
        }
        if(data[i].output_image_data != NULL)
        {
            free(data[i].output_image_data);
            data[i].output_image_data = NULL;
        }
    }
    error09:
    free(file_name);
    error08:
    free(superset_idx);
    error07:
    free(superset_diff);
    error06:
    free(superset_gene[0]);
    error05:
    free(superset_gene);
    error04:
    free(data);
    error03:
    free(lp_urandom);
    error02:
    close(urandom_fd);
    error01:
    return ret_code;
}
