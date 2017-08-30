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
#define num_parent 256
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
#define layer 24
#define input_datasize 3
#define output_datasize 3
#define hidden_channels 64
#define pix_ref 9
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
    {hidden_channels,hidden_channels,1},
    {hidden_channels,hidden_channels,1},
    {hidden_channels,output_datasize,1}
};

int urandom_fd;
uint64_t *lp_urandom;

int get_random(uint64_t *data,int size) {
    static int i=UrandomBuff;
    int j,k,ret_code;
    ret_code = 0;
    if(lp_urandom == NULL || urandom_fd < 0) {
        ret_code = -1;
        goto urandom_error;
    }
    
    k=0;
    while(k < size) {
        if(i == UrandomBuff) {
            j = read(urandom_fd,lp_urandom,UrandomBuff*sizeof(uint64_t));
            if(j != UrandomBuff*sizeof(uint64_t)) {
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

int get_sud(RealNumber *data,int size) {
    static int i=2;
    static int64_t n[2];
    static RealNumber f[2];
    int k,ret_code;
    ret_code = 0;
    
    k=0;
    while(k < size) {
        if(i == 2) {
            if(get_random((uint64_t *)n,2)) {
                ret_code = -1;
                goto error;
            }
            f[0] = ((RealNumber)n[0])/((RealNumber)INT64_MAX);
            f[1] = ((RealNumber)n[1])/((RealNumber)INT64_MAX);
            i=0;
        }
        data[k++] = f[i++];
    }

    error:
    return ret_code;
}

int creat_gene(char *filename) {
    int i,j,k,l,fd,m,element_size;
    RealNumber *b,temp;
    element_size = 0;
    for(i=0;i<layer;i++) {
        element_size += (pix_ref * matrix_size[i][0] + aux_vector_size)*matrix_size[i][1];
    }
    
    b = (RealNumber *)malloc(element_size*sizeof(RealNumber));
    get_sud(b,element_size);

    i=0;
    j=0;
    for(k=0;k<layer;k++) {
        j += pix_ref * matrix_size[k][0]*matrix_size[k][1];
        for(;i<j;i++) {
            b[i] *= 1.4/((RealNumber)pix_ref);
        }
        j += aux_vector_size*matrix_size[k][1];
        for(;i<j;i++) {
            b[i] *= 0.2;
        }
    }
    k=0;
    for(i=-1;i<=1;i++) {
        for(j=-1;j<=1;j++) {
            if(i==0 && j == 0) {
                l=k;
            }
            k++;
        }
    }
    
    i=0;
    j=0;
    for(k=0;k<layer;k++) {
        for(i=0;i<matrix_size[k][0] && i<matrix_size[k][1];i++) {
            b[i + i*matrix_size[k][1] + l*matrix_size[k][0]*matrix_size[k][1] +j] *= 0.05;
            b[i + i*matrix_size[k][1] + l*matrix_size[k][0]*matrix_size[k][1] +j] += 1.0;
        }
        j += (pix_ref * matrix_size[k][0] + aux_vector_size)*matrix_size[k][1];
    }

    fd = open(filename,O_WRONLY|O_CREAT|O_TRUNC,0666);
    if(fd < 0) {
        return -1;
    }
    i=write(fd,b,element_size*sizeof(RealNumber));
    free(b);
    close(fd);
    return 0;
}
int main() {
    int i,ret_code;
    char *filename;
    ret_code = 0;
    urandom_fd = open("/dev/urandom",O_RDONLY);
    if(urandom_fd < 0) {
        ret_code = -1;
        goto urandom_failed1;
    }
    lp_urandom = (uint64_t *)malloc(UrandomBuff*sizeof(uint64_t));
    if(lp_urandom == NULL) {
        ret_code = -1;
        goto urandom_failed2;
    }
    filename = (char *)malloc(256*sizeof(char));
    if(filename == NULL) {
        ret_code = -1;
        goto malloc_filename_error;
    }

    for(i=0;i<num_parent;i++) {
        sprintf(filename,"./gene/gene%08d.bin",i);
        creat_gene(filename);
    }

    free(filename);
    malloc_filename_error:
    free(lp_urandom);
    urandom_failed2:
    close(urandom_fd);
    urandom_failed1:
    return ret_code;
}
