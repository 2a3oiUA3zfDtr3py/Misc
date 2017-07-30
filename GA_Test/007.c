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
#include <wand/MagickWand.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define num_generate 3200
#define num_parent 80
#define gene_pool_size 400
#define num_data_set 6
#define data_pool_size 90
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
#define layer 8
#define input_datasize 3
#define output_datasize 3
#define hidden_channels 32
#define pix_ref 17
int matrix_size[layer][3] = {
    {input_datasize,hidden_channels,2},
    {hidden_channels,hidden_channels,4},
    {hidden_channels,hidden_channels,8},
    {hidden_channels,hidden_channels,16},
    {hidden_channels,hidden_channels,16},
    {hidden_channels,hidden_channels,8},
    {hidden_channels,hidden_channels,4},
    {hidden_channels,output_datasize,2}};
int urandom_fd;
uint64_t *lp_urandom;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t min_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

int vector_x_matrix(RealNumber *a,int a_size,RealNumber **b,int b_row,int b_col,RealNumber *c,int c_size) {
    int i,j;
    for(j=0;j<b_col;j++) {
        c[j] =0.0;
    }
    for(i=0;i<b_row;i++) {
        for(j=0;j<b_col;j++) {
            c[j] += b[i][j] * a[i];
        }
    }
    return 0;
}

int calc(RealNumber *gene,RealNumber **data,int width,int height,RealNumber *diff,RealNumber ***a,RealNumber **b,RealNumber *c) {
    int return_code;
    int i,j,k,l,m,n;
    int p[17][2];
    RealNumber temp;
    return_code = 0;

    a[0][0] = gene;
    for(j=1;j<pix_ref*matrix_size[0][0]+aux_vector_size;j++) {
        a[0][j] = a[0][j-1] + matrix_size[0][1];
    }
    for(i=1;i<layer;i++) {
        a[i] = a[i-1] + pix_ref*matrix_size[i-1][0] + aux_vector_size;
        a[i][0] = a[i-1][pix_ref*matrix_size[i-1][0] + aux_vector_size -1] + matrix_size[i-1][1];
        for(j=1;j<pix_ref*matrix_size[i][0] + aux_vector_size;j++) {
            a[i][j] = a[i][j-1] + matrix_size[i][1];
        }
    }
    
    k=0;
    for(i=-1;i<=1;i++) {
        for(j=-1;j<=1;j++) {
            p[k][0] = i;
            p[k][1] = j;
            k++;
        }
    }

    for(i=0;i<height*width*matrix_size[0][0];i++) {
        b[0][i] = data[0][i];
    }
    l=9;
    for(j=-1;j<=1;j++) {
        for(k=-1;k<=1;k++) {
            if(j!=0 || k!=0) {
                p[l][0] = j*matrix_size[0][2];
                p[l][1] = k*matrix_size[0][2];
                l++;
            }
        }
    }
    for(j=0;j<height;j++) {
        for(k=0;k<width;k++) {
            for(n=0;n<pix_ref;n++) {
                if(j+p[n][0] >= 0 && j+p[n][0] < height && k+p[n][1] >= 0 && k+p[n][1] < width) {
                    for(m=0;m<matrix_size[0][0];m++) {
                        c[matrix_size[0][0]*n + m] = b[0][matrix_size[0][0]*((j+p[n][0])*width+(k+p[n][1]))+m];
                    }
                }else{
                    for(m=0;m<matrix_size[0][0];m++) {
                            c[matrix_size[0][0]*n + m] = -1.0;
                    }
                }
            }
            vector_x_matrix(c,pix_ref*matrix_size[0][0],a[0],pix_ref*matrix_size[0][0],matrix_size[0][1],
                            c+pix_ref*hidden_channels,matrix_size[0][1]);
            for(n=0;n<matrix_size[0][1];n++) {
                c[pix_ref*hidden_channels+n] += a[0][pix_ref*matrix_size[0][0]][n];
                b[1][matrix_size[0][1]*(j*width+k)+n] = c[pix_ref*hidden_channels+n];
            }
        }
    }
    
    
    for(i=1;i<layer;i++) {
        l=9;
        for(j=-1;j<=1;j++) {
            for(k=-1;k<=1;k++) {
                if(j!=0 || k!=0) {
                    p[l][0] = j*matrix_size[i][2];
                    p[l][1] = k*matrix_size[i][2];
                    l++;
                }
            }
        }
        for(j=0;j<height*width*matrix_size[i][0];j++) {
            if(b[i][j] < 0.0) {
                b[i][j] = 0.0;
            }
        }
        for(j=0;j<height;j++) {
            for(k=0;k<width;k++) {
                for(n=0;n<pix_ref;n++) {
                    if(j+p[n][0] < 0) {
                        if(k+p[n][1] < 0) {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][m];
                            }
                        }else if(k+p[n][1] < width) {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*(k+p[n][1])+m];
                            }
                        }else {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*(width-1)+m];
                            }
                        }
                    }else if(j+p[n][0] < height) {
                        if(k+p[n][1] < 0) {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((j+p[n][0])*width)+m];
                            }
                        }else if(k+p[n][1] < width) {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((j+p[n][0])*width+(k+p[n][1]))+m];
                            }
                        }else {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((j+p[n][0])*width+(width-1))+m];
                            }
                        }
                    }else{
                        if(k+p[n][1] < 0) {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((height-1)*width)+m];
                            }
                        }else if(k+p[n][1] < width) {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((height-1)*width+(k+p[n][1]))+m];
                            }
                        }else {
                            for(m=0;m<matrix_size[i][0];m++) {
                                c[matrix_size[i][0]*n + m] = b[i][matrix_size[i][0]*((height-1)*width+(width-1))+m];
                            }
                        }
                    }
                }
                vector_x_matrix(c,pix_ref*matrix_size[i][0],a[i],pix_ref*matrix_size[i][0],matrix_size[i][1],
                                c+pix_ref*hidden_channels,matrix_size[i][1]);
                for(n=0;n<matrix_size[i][1];n++) {
                    c[pix_ref*hidden_channels+n] += a[i][pix_ref*matrix_size[i][0]][n];
                    b[i+1][matrix_size[i][1]*(j*width+k)+n] = c[pix_ref*hidden_channels+n];
                }
            }
        }
    }
    
    temp = 0.0;
    for(j=0;j<height*width*matrix_size[layer-1][1];j++) {
        b[i][j] -= data[1][j];
        temp += b[i][j]*b[i][j];
    }
    if(isfinite(temp)) {
        temp /= (RealNumber)(width*height*output_datasize);
    }
    if(isfinite(temp)) {
        *diff += temp;
    } else {
        *diff = -1.0;
    }
    

    if(return_code) {
        *diff = -1.0;
    }
    return return_code;
}

int R_sort(RealNumber *score,int *idx,int size) {
    int i,j,k,l,m,n,m_max,n_max,ret_val;
    int *index[2];
    ret_val = 0;
    if(idx == NULL || score == NULL) {
        ret_val = -1;
        goto failed_input;
    }
    index[0] = idx;
    index[1] = (int *)malloc(size*sizeof(int));
    if(index[1] == NULL){
        ret_val = -1;
        goto failed_index;
    }
    for(i=0;i<size;i++) {
        index[0][i] = i;
    }
    
    l=0;
    for(i=1;i<size;i*=2) {
        for(j=0;j<size;j+=i*2) {
            k=j;
            m=j;
            n=j+i;
            m_max = (m+i<size)? m+i:size;
            n_max = (n+i<size)? n+i:size;
            while(m<m_max && n<n_max) {
                if(score[index[l][m]] >= RealNumber_MIN && score[index[l][n]] >= RealNumber_MIN) {
                    if(score[index[l][m]] < score[index[l][n]]) {
                        index[l^0x01][k++] = index[l][m++];
                    }else{
                        index[l^0x01][k++] = index[l][n++];
                    }
                }else if(score[index[l][m]] < -RealNumber_MIN && score[index[l][n]] < -RealNumber_MIN) {
                    index[l^0x01][k++] = index[l][m++];
                    index[l^0x01][k++] = index[l][n++];
                }else if(score[index[l][n]] < -RealNumber_MIN) {
                    index[l^0x01][k++] = index[l][m++];
                }else{
                    index[l^0x01][k++] = index[l][n++];
                }
            }
            while(m<m_max) {
                index[l^0x01][k++] = index[l][m++];
            }
            while(n<n_max) {
                index[l^0x01][k++] = index[l][n++];
            }
        }
        l ^= 0x01;
    }
    if(l == 1) {
        for(i=0;i<size;i++) {
            idx[i] = index[1][i];
        }
    }
    free(index[1]);
    failed_index:
    failed_input:
    return ret_val;
}

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

int get_sd(RealNumber *data,int size) {
    static int i=2;
    static uint64_t n[2];
    static RealNumber f[2];
    RealNumber temp[2];
    int k,ret_code;
    ret_code = 0;
    
    k=0;
    while(k < size) {
        if(i == 2) {
            if(get_random(n,2)) {
                ret_code = -1;
                goto error;
            }
            if(n[0] == 0) {
                n[0]++;
            }else if(n[0] == UINT64_MAX) {
                n[0]--;
            }
            if(n[1] == 0) {
                n[1]++;
            }else if(n[1] == UINT64_MAX) {
                n[1]--;
            }
            temp[0] = ((RealNumber)n[0])/((RealNumber)UINT64_MAX);
            temp[1] = ((RealNumber)n[1])/((RealNumber)UINT64_MAX);
            temp[0] = sqrt(-2.0 * log(temp[0]));
            temp[1] *= 2.0 * M_PI;
            if(!isfinite(temp[0])) {
                ret_code = -1;
                goto error;
            }
            if(!isfinite(temp[1])) {
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

int random_list(int *idx,int size) {
    int i,j,k,l,m,n,m_max,n_max,ret_code;
    int *index[2];
    uint64_t *r;
    ret_code = 0;
    if(idx == NULL || size <= 0) {
        ret_code = -1;
        goto input_error;
    }
    index[0] = idx;
    index[1] = (int *)malloc(size*sizeof(int));
    if(index[1] == NULL) {
        ret_code = -1;
        goto malloc_error1;
    }
    r = (uint64_t *)malloc(size*sizeof(uint64_t));
    if(r == NULL) {
        ret_code = -1;
        goto malloc_error2;
    }
    i = get_random(r,size);
    if(i) {
        ret_code = -1;
        goto urandom_error;
    }
    for(i=0;i<size;i++) {
        index[0][i] = i;
    }
    l=0;
    for(i=1;i<size;i*=2) {
        for(j=0;j<size;j+=i*2) {
            k=j;
            m=j;
            n=j+i;
            m_max = (m+i<size)? m+i:size;
            n_max = (n+i<size)? n+i:size;
            while(m<m_max && n<n_max) {
                if(r[index[l][m]] < r[index[l][n]]) {
                    index[l^0x01][k++] = index[l][m++];
                }else{
                    index[l^0x01][k++] = index[l][n++];
                }
            }
            while(m<m_max) {
                index[l^0x01][k++] = index[l][m++];
            }
            while(n<n_max) {
                index[l^0x01][k++] = index[l][n++];
            }
        }
        l ^= 0x01;
    }
    if(l == 1) {
        for(i=0;i<size;i++) {
            idx[i] = index[1][i];
        }
    }
    urandom_error:
    free(r);
    malloc_error2:
    free(index[1]);
    malloc_error1:
    input_error:
    return ret_code;
}

int gene_generator(RealNumber **a,RealNumber **b,int *parent_gene_idx,RealNumber *parent_gene_diff) {
    int i,j,k,l,element_size,ret_code;
    uint64_t *urdm;
    int64_t *srdm;
    RealNumber **temp;
    RealNumber alpha,beta,gamma,delta,*tiger;
    ret_code = 0;
    element_size = 0;
    for(i=0;i<layer;i++) {
        element_size += (pix_ref * matrix_size[i][0] + aux_vector_size) * matrix_size[i][1];
    }
    if(a == NULL || b ==NULL) {
        ret_code = -1;
        goto input_error;
    }
    if(b[0] == NULL) {
        ret_code = -1;
        goto input_error;
    }

    urdm = (uint64_t *)malloc(2*element_size*sizeof(uint64_t));
    if(urdm == NULL) {
        ret_code = -1;
        goto malloc_urdm_error;
    }
    srdm = (int64_t *)malloc(2*element_size*sizeof(int64_t));
    if(srdm == NULL) {
        ret_code = -1;
        goto malloc_srdm_error;
    }
    temp = (RealNumber **)malloc((num_parent+4)*sizeof(RealNumber *));
    if(temp == NULL) {
        ret_code = -1;
        goto malloc_temp_error1;
    }
    temp[0] = (RealNumber *)malloc((num_parent+4)*element_size*sizeof(RealNumber));
    if(temp[0] == NULL) {
        ret_code = -1;
        goto malloc_temp_error2;
    }
    for(i=1;i<num_parent+4;i++) {
        temp[i] = temp[i-1] + element_size;
    }
    i = (num_parent > element_size)? num_parent : element_size;
    tiger = (RealNumber *)malloc(i*sizeof(RealNumber));
    if(tiger == NULL) {
        ret_code = -1;
        goto malloc_tiger_error;
    }
    
    for(j=0;j<element_size;j++) {
        temp[num_parent][j] = 0.0;
        temp[num_parent+1][j] = 0.0;
        temp[num_parent+2][j] = 0.0;
        temp[num_parent+3][j] = 0.0;
    }
    for(i=0;i<num_parent;i++) {
        for(j=0;j<element_size;j++) {
            temp[num_parent][j] += a[i][j];
        }
    }
    for(j=0;j<element_size;j++) {
        temp[num_parent][j] /= num_parent;
    }
    for(i=0;i<num_parent;i++) {
        for(j=0;j<element_size;j++) {
            temp[num_parent+1][j] += (a[i][j] - temp[num_parent][j])*(a[i][j] - temp[num_parent][j]);
        }
    }
    for(k=0;k<element_size;k++) {
        temp[num_parent+1][k] /= num_parent;
        temp[num_parent+1][k] = sqrt(temp[num_parent+1][k]);
    }
    k=num_parent/8;
    beta = 0.0;
    for(i=0;i<k;i++) {
        alpha = ((RealNumber)(k-i))/((RealNumber)k);
        alpha *= alpha;
        for(j=0;j<element_size;j++) {
            temp[num_parent+2][j] += alpha * a[parent_gene_idx[i]][j];
        }
        beta += alpha;
    }
    for(j=0;j<element_size;j++) {
        temp[num_parent+2][j] /= beta;
    }
    beta = 0.0;
    for(i=0;i<k;i++) {
        alpha = ((RealNumber)(i+1))/((RealNumber)k);
        alpha *= alpha;
        for(j=0;j<element_size;j++) {
            temp[num_parent+3][j] += alpha * a[parent_gene_idx[i]][j];
        }
        beta += alpha;
    }
    for(j=0;j<element_size;j++) {
        temp[num_parent+3][j] /= beta;
    }
    delta = 0.7;
    alpha = parent_gene_diff[parent_gene_idx[0]]/parent_gene_diff[parent_gene_idx[num_parent/2]];
    alpha = alpha / delta;
    if(!isfinite(alpha)) {
        alpha = 1.0;
    }
    delta = 0.5*(alpha - 1.0)*delta/(delta - 1.0) + 1.0;
    if(alpha < 1.0){
        alpha = 1.0;
    }
    if(delta > 1.0)
    {
        delta = 1.0;
    }
    if(delta < 0.5)
    {
        delta = 0.5;
    }
    printf("%f\n%f\n",alpha,delta);
    gamma = alpha;
    alpha *= sqrt(1.0/num_parent);

    for(i=0;i<num_parent;i++) {
        for(j=0;j<element_size;j++) {
            temp[i][j] = a[i][j] - temp[num_parent][j];
        }
    }
    l = num_generate / 3;
    for(i=0;i<l*1;i++) {
        get_sd(tiger,num_parent);
        for(j=0;j<element_size;j++) {
            b[i][j] = 0;
        }
        for(k=0;k<num_parent;k++) {
            for(j=0;j<element_size;j++) {
                b[i][j] += alpha * tiger[k] * temp[k][j];
            }
        }
        get_random(urdm,1);
        beta = ((RealNumber)urdm[0])/((RealNumber)UINT64_MAX);
        if(beta < 0.02) {
            get_random(urdm,element_size);
            get_sd(tiger,element_size);
            for(j=0;j<element_size;j++) {
                beta = ((RealNumber)urdm[j])/((RealNumber)UINT64_MAX);
                if(beta < 0.01) {
                    b[i][j] += 2.0 * gamma * tiger[j] * temp[num_parent+1][j];
                }
            }
        }else if(beta < 0.04) {
            get_random(urdm,1);
            get_sd(tiger,1);
            urdm[0] %= num_parent;
            for(j=0;j<element_size;j++) {
                b[i][j] += 2.0 * gamma * tiger[0] * temp[urdm[0]][j];
            }
        }
        get_sd(tiger,1);
        tiger[0] *= 2.0;
        for(j=0;j<element_size;j++) {
            b[i][j] += (delta+tiger[0]) * temp[num_parent+2][j] + (1.0 -delta -tiger[0]) * temp[num_parent+3][j];
        }
    }
    for(i=l*1;i<l*2;i++) {
        get_sd(tiger,num_parent);
        for(j=0;j<element_size;j++) {
            b[i][j] = 0;
        }
        for(k=0;k<num_parent;k++) {
            for(j=0;j<element_size;j++) {
                b[i][j] += alpha * tiger[k] * temp[k][j];
            }
        }
        get_sd(tiger,element_size);
        get_random(urdm,element_size);
        for(j=0;j<element_size;j++) {
            beta = ((RealNumber)urdm[j])/((RealNumber)UINT64_MAX);
            if(beta < 0.5) {
                tiger[j] *= 0.9;
                tiger[j] += 1.0;
                if(beta < 0.25) {
                    tiger[j] *= -1.0;
                }
                b[i][j] = gamma * tiger[j] * temp[num_parent+1][j];
            }
        }
        get_sd(tiger,1);
        tiger[0] *= 2.0;
        for(j=0;j<element_size;j++) {
            b[i][j] += (delta+tiger[0]) * temp[num_parent+2][j] + (1.0 -delta -tiger[0]) * temp[num_parent+3][j];
        }
    }
    for(i=l*2;i<num_generate;i++) {
        get_sd(tiger,element_size);
        get_random(urdm,element_size);
        for(j=0;j<element_size;j++) {
            beta = ((RealNumber)urdm[j])/((RealNumber)UINT64_MAX);
            tiger[j] *= 0.9;
            tiger[j] += 1.0;
            if(beta < 0.5) {
                tiger[j] *= -1.0;
            }
            b[i][j] = gamma * tiger[j] * temp[num_parent+1][j];
        }
        get_sd(tiger,1);
        tiger[0] *= 2.0;
        for(j=0;j<element_size;j++) {
            b[i][j] += (delta+tiger[0]) * temp[num_parent+2][j] + (1.0 -delta -tiger[0]) * temp[num_parent+3][j];
        }
    }
    
    free(tiger);
    malloc_tiger_error:
    free(temp[0]);
    malloc_temp_error2:
    free(temp);
    malloc_temp_error1:
    free(srdm);
    malloc_srdm_error:
    free(urdm);
    malloc_urdm_error:
    input_error:
    return ret_code;
}

int get_filedata(char *filename,RealNumber **a,int *width,int *height) {
    int ret_code;
    MagickWand *img_handle;
    MagickBooleanType status;
    ret_code = 0;
    if(filename == NULL || a == NULL || width == NULL || height == NULL) {
        goto input_error;
    }
    if(*a != NULL) {
        free(*a);
        *a = NULL;
    }
    MagickWandGenesis();
    img_handle = NewMagickWand();
    status = MagickReadImage(img_handle,filename);
    if(status == MagickFalse) {
        ret_code = -1;
        goto error;
    }
    *width = (int)MagickGetImageWidth(img_handle);
    *height = (int)MagickGetImageHeight(img_handle);
    *a = (RealNumber *)malloc(3*(*width)*(*height)*sizeof(RealNumber));
    if(*a == NULL) {
        ret_code = -1;
        goto error;
    }
#ifdef DOUBLE
    status = MagickExportImagePixels(img_handle,0,0,(unsigned long)(*width),(unsigned long)(*height),"RGB",DoublePixel,*a);
#else
    status = MagickExportImagePixels(img_handle,0,0,(unsigned long)(*width),(unsigned long)(*height),"RGB",FloatPixel,*a);
#endif
    if(status == MagickFalse) {
        ret_code = -1;
        goto error;
    }
    error:
    if(ret_code && *a != NULL) {
        free(*a);
        *a = NULL;
    }
    DestroyMagickWand(img_handle);
    MagickWandTerminus();
    input_error:
    return ret_code;
}

void* calc_p(void *arg) {
    void **args = (void *)arg;
    RealNumber **gene  = (RealNumber **)args[0];
    RealNumber **data = (RealNumber **)args[1];
    int width = *(int *)args[2];
    int height = *(int *)args[3];
    RealNumber *diff  = (RealNumber  *)args[4];
    int length = *(int *)args[5];
    int *threadnum = (int *)args[6];
    int i,j;
    RealNumber ***a,**b,*c;
    
    a = (RealNumber ***)malloc(layer*sizeof(RealNumber **));
    if(a == NULL) {
        goto malloc_a_error1;
    }
    j = 0;
    for(i=0;i<layer;i++) {
        j += pix_ref*matrix_size[i][0] + aux_vector_size;
    }
    a[0] = (RealNumber **)malloc(j*sizeof(RealNumber *));
    if(a[0] == NULL) {
        goto malloc_a_error2;
    }
    b = (RealNumber **)malloc((layer+1)*sizeof(RealNumber *));
    if(b == NULL) {
        goto malloc_b_error1;
    }
    b[0] = (RealNumber *)malloc(2*hidden_channels*width*height*sizeof(RealNumber));
    if(b[0] == NULL) {
        goto malloc_b_error2;
    }
    b[1] = b[0] + hidden_channels*width*height;
    for(i=2;i<layer+1;i++) {
        b[i] = b[i&0x01];
    }
    c = (RealNumber *)malloc(2*pix_ref*hidden_channels*sizeof(RealNumber));
    if(c == NULL) {
        goto malloc_c_error;
    }
    
    for(i=0;i<length;i++) {
        calc(gene[i],data,width,height,diff + i,a,b,c);
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

int gene_calc(RealNumber **gene,RealNumber *diff,RealNumber **data,int width,int heigth,int length) {
    int i,j,k,l,m,n,o,element_size,ret_code;
    int **threadnum,*threads_gene_number,**threads_gene_length;
    void ***lp;
    pthread_t *threads;
    RealNumber **internal_diff;
    ret_code = 0;
    element_size = 0;
    for(i=0;i<layer;i++) {
        element_size += (pix_ref * matrix_size[i][0] + aux_vector_size) * matrix_size[i][1];
    }
    lp = (void ***)malloc(ThreadQuantity*sizeof(void **));
    if(lp == NULL) {
        ret_code = -1;
        goto malloc_error1;
    }
    for(i=0;i<ThreadQuantity;i++) {
        lp[i] = (void **)malloc(9*sizeof(void *));
        if(lp[i] == NULL) {
            ret_code = -1;
            goto malloc_error2;
        }
    }
    threads = (pthread_t *)malloc(ThreadQuantity*sizeof(pthread_t));
    if(threads == NULL) {
        ret_code = -1;
        goto malloc_error2;
    }
    internal_diff = (RealNumber **)malloc(ThreadQuantity*sizeof(RealNumber *));
    if(internal_diff == NULL) {
        ret_code = -1;
        goto malloc_error3;
    }
    for(i=0;i<ThreadQuantity;i++) {
        internal_diff[i] = (RealNumber *)malloc(length*sizeof(RealNumber));
        if(internal_diff[i] == NULL) {
            ret_code = -1;
            goto malloc_error4;
        }
    }
    threadnum = (int **)malloc(ThreadQuantity*sizeof(int *));
    if(threadnum == NULL) {
        ret_code = -1;
        goto malloc_error4;
    }
    for(i=0;i<ThreadQuantity;i++) {
        threadnum[i] = (int *)malloc(sizeof(int));
        if(threadnum[i] == NULL) {
            ret_code = -1;
            goto malloc_error5;
        }
    }
    threads_gene_number = (int *)malloc(ThreadQuantity*sizeof(int));
    if(threads_gene_number == NULL) {
        ret_code = -1;
        goto malloc_error5;
    }
    threads_gene_length = (int **)malloc(ThreadQuantity*sizeof(int *));
    if(threads_gene_length == NULL) {
        ret_code = -1;
        goto malloc_error6;
    }
    for(i=0;i<ThreadQuantity;i++) {
        threads_gene_length[i] = (int *)malloc(sizeof(int));
        if(threads_gene_length[i] == NULL) {
            ret_code = -1;
            goto malloc_error7;
        }
    }
    
    for(i=0;i<ThreadQuantity;i++) {
        lp[i][1] = (void *)data;
        lp[i][2] = (void *)(&width);
        lp[i][3] = (void *)(&heigth);
        lp[i][4] = (void *)internal_diff[i];
        lp[i][5] = (void *)threads_gene_length[i];
        lp[i][6] = (void *)threadnum[i];
        threadnum[i][0] = -2;
    }
    
    k = 0;
    while(1) {
        i = (length - k)/(ThreadQuantity + 1);
        if(i == 0) {
            break;
        }
        for(l=0;l<ThreadQuantity;l++) {
            pthread_mutex_lock(&mutex);
            while(1) {
                m=0;
                for(n=0;n<ThreadQuantity;n++) {
                    m=threadnum[n][0];
                    if(m) {
                        break;
                    }
                }
                if(m) {
                    break;
                } else {
                    pthread_cond_wait(&cond, &mutex);
                }
            }
            pthread_mutex_unlock(&mutex);
            if(m==-1) {
                pthread_join(threads[n],NULL);
                for(o=0;o<threads_gene_length[n][0];o++) {
                    diff[o+threads_gene_number[n]] = internal_diff[n][o];
                }
            }
            lp[n][0] = gene + k;
            threads_gene_number[n] = k;
            threads_gene_length[n][0] = i;
            threadnum[n][0] = 0;
            for(o=0;o<threads_gene_length[n][0];o++) {
                internal_diff[n][o] = diff[o+threads_gene_number[n]];
            }
            pthread_create(threads+n,NULL,(void *)calc_p,(void *)lp[n]);
            k += i;
        }
    }
    for(i=k;i<length;i++) {
        pthread_mutex_lock(&mutex);
        while(1) {
            m=0;
            for(n=0;n<ThreadQuantity;n++) {
                m=threadnum[n][0];
                if(m) {
                    break;
                }
            }
            if(m) {
                break;
            }else {
                pthread_cond_wait(&cond, &mutex);
            }
        }
        pthread_mutex_unlock(&mutex);
        if(m == -1) {
            pthread_join(threads[n],NULL);
            pthread_mutex_lock(&mutex);
            threadnum[n][0] = -2;
            pthread_mutex_unlock(&mutex);
            for(o=0;o<threads_gene_length[n][0];o++) {
                diff[o+threads_gene_number[n]] = internal_diff[n][o];
            }
        }
        lp[n][0] = gene + i;
        threads_gene_number[n] = i;
        threads_gene_length[n][0] = 1;
        threadnum[n][0] = 0;
        for(o=0;o<threads_gene_length[n][0];o++) {
            internal_diff[n][o] = diff[o+threads_gene_number[n]];
        }
        pthread_create(threads+n,NULL,(void *)calc_p,(void *)lp[n]);
    }
    for(n=0;n<ThreadQuantity;n++) {
        pthread_mutex_lock(&mutex);
        m=threadnum[n][0];
        pthread_mutex_unlock(&mutex);
        if(m != -2) {
            pthread_join(threads[n],NULL);
            for(o=0;o<threads_gene_length[n][0];o++) {
                diff[o+threads_gene_number[n]] = internal_diff[n][o];
            }
        }
    }
    while(1) {
        pthread_mutex_lock(&mutex);
        j=-1;
        while(1) {
            m=0;
            for(n=0;n<ThreadQuantity;n++) {
                m=threadnum[n][0];
                if(m == 0) {
                    j=0;
                }else if(m == -1) {
                    break;
                }
            }
            if(m == -1  || (j!=0  && n == ThreadQuantity)) {
                break;
            }else{
                pthread_cond_wait(&cond, &mutex);
            }
        }
        pthread_mutex_unlock(&mutex);
        if(m==-1) {
            pthread_join(threads[n],NULL);
            threadnum[n][0] = -2;
            for(o=0;o<threads_gene_length[n][0];o++) {
                diff[o+threads_gene_number[n]] = internal_diff[n][o];
            }
        }
        if(j && n == ThreadQuantity) {
            break;
        }
    }

    malloc_error7:
    for(i=0;i<ThreadQuantity;i++) {
        if(threads_gene_length[i] != NULL) {
            free(threads_gene_length[i]);
        }else{
            break;
        }
    }
    free(threads_gene_length);
    malloc_error6:
    free(threads_gene_number);
    malloc_error5:
    for(i=0;i<ThreadQuantity;i++) {
        if(threadnum[i] != NULL) {
            free(threadnum[i]);
        }else{
            break;
        }
    }
    free(threadnum);
    malloc_error4:
    for(i=0;i<ThreadQuantity;i++) {;
        if(internal_diff[i] != NULL) {
            free(internal_diff[i]);
        }else{
            break;
        }
    }
    free(internal_diff);
    malloc_error3:
    free(threads);
    malloc_error2:
    for(i=0;i<ThreadQuantity;i++) {
        if(lp[i] != NULL) {
            free(lp[i]);
        }else{
            break;
        }
    }
    free(lp);
    malloc_error1:
    return ret_code;
}

int save_gene(char *filename,RealNumber *gene,int element_size) {
    int fd,length,ret_code;
    ret_code = 0;
    fd = open(filename,O_WRONLY|O_CREAT|O_TRUNC,0666);
    if(fd < 0) {
        ret_code = -1;
        goto file_open_error;
    }
    length = write(fd,gene,element_size * sizeof(RealNumber));
    if(length != element_size * sizeof(RealNumber)) {
        ret_code = -1;
    }
    
    close(fd);
    file_open_error:
    return ret_code;
}

int load_gene(char *filename,RealNumber *gene,int element_size) {
    int i,fd,ret_code;
    struct stat file_stat;
    ret_code = 0;
    fd = open(filename,O_RDONLY);
    if(fd < 0) {
        ret_code = -1;
        goto file_error;
    }
    if(fstat(fd,&file_stat)) {
        ret_code = -1;
        close(fd);
        goto file_error;
    }
    if(file_stat.st_size != element_size*sizeof(RealNumber)) {
        ret_code = -1;
        close(fd);
        goto file_error;
    }
    i=read(fd,gene,file_stat.st_size);
    close(fd);
    if(i != file_stat.st_size) {
        ret_code = -1;
        goto file_error;
    }
    file_error:
    return ret_code;
}

int main() {
    int i,j,k,l,m,element_size,ret_code;
    int **file_size,*p;
    char *file_name;
    RealNumber ***data,**parent_gene,**child_gene;
    RealNumber *parent_gene_diff,*child_gene_diff,f_diff,g_diff;
    int *gene_pool_idx,*parent_gene_idx,*child_gene_idx;

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

    ret_code = 0;
    element_size = 0;
    for(i=0;i<layer;i++) {
        element_size += (pix_ref * matrix_size[i][0] + aux_vector_size) * matrix_size[i][1];
    }
    file_size = (int **)malloc(num_data_set*sizeof(int *));
    if(file_size == NULL) {
        ret_code = -1;
        goto malloc_file_size_error1;
    }
    data = (RealNumber ***)malloc(num_data_set*sizeof(RealNumber **));
    if(data == NULL) {
        ret_code = -1;
        goto malloc_data_error1;
    }
    parent_gene = (RealNumber **)malloc(num_parent*sizeof(RealNumber *));
    if(parent_gene == NULL) {
        ret_code = -1;
        goto malloc_parent_gene_error1;
    }
    child_gene = (RealNumber **)malloc(num_generate*sizeof(RealNumber *));
    if(child_gene == NULL) {
        ret_code = -1;
        goto malloc_child_gene_error1;
    }
    parent_gene_diff = (RealNumber *)malloc(num_parent*sizeof(RealNumber));
    if(parent_gene_diff == NULL) {
        ret_code = -1;
        goto malloc_parent_gene_diff_error;
    }
    child_gene_diff = (RealNumber *)malloc(num_generate*sizeof(RealNumber));
    if(child_gene_diff == NULL) {
        ret_code = -1;
        goto malloc_child_gene_diff_error;
    }
    gene_pool_idx = (int *)malloc(gene_pool_size*sizeof(int));
    if(gene_pool_idx == NULL) {
        ret_code = -1;
        goto malloc_gene_pool_idx_error;
    }
    parent_gene_idx = (int *)malloc(num_parent*sizeof(int));
    if(parent_gene_idx == NULL) {
        ret_code = -1;
        goto malloc_parent_gene_idx_error;
    }
    child_gene_idx = (int *)malloc(num_generate*sizeof(int));
    if(child_gene_idx == NULL) {
        ret_code = -1;
        goto malloc_child_gene_idx_error;
    }
    file_size[0] = (int *)malloc(num_data_set*4*sizeof(int));
    if(file_size[0] == NULL) {
        ret_code = -1;
        goto malloc_file_size_error2;
    }
    data[0] = (RealNumber **)malloc(num_data_set*2*sizeof(RealNumber *));
    if(data[0] == NULL) {
        ret_code = -1;
        goto malloc_data_error2;
    }
    parent_gene[0] = (RealNumber *)malloc(element_size*num_parent*sizeof(RealNumber));
    if(parent_gene[0] == NULL) {
        ret_code = -1;
        goto malloc_parent_gene_error2;
    }
    child_gene[0] = (RealNumber *)malloc(element_size*num_generate*sizeof(RealNumber));
    if(child_gene[0] == NULL) {
        ret_code = -1;
        goto malloc_child_gene_error2;
    }
    file_name = (char *)malloc(256*sizeof(char));
    if(file_name == NULL) {
        ret_code = -1;
        goto malloc_file_name_error;
    }
    p = (int *)malloc(data_pool_size*sizeof(int));
    if(p == NULL) {
        ret_code = -1;
        goto malloc_p_error;
    }
    j=0;
    for(i=0;i<layer;i++) {
        j+= matrix_size[i][0];
    }
    for(i=1;i<num_data_set;i++) {
        file_size[i] = file_size[i-1] + 4;
    }
    for(i=1;i<num_data_set;i++) {
        data[i] = data[i-1] + 2;
    }
    for(i=1;i<num_parent;i++) {
        parent_gene[i] = parent_gene[i-1] + element_size;
    }
    for(i=1;i<num_generate;i++) {
        child_gene[i] = child_gene[i-1] + element_size;
    }
    
    fcntl(0,F_SETFL,O_NONBLOCK);
    while(1) {
        random_list(p,data_pool_size);
        random_list(gene_pool_idx,gene_pool_size);
        for(i=0;i<num_data_set;i++) {
            sprintf(file_name,"./data/%08da.png",p[i]);
            if(get_filedata(file_name,data[i] + 0,file_size[i] + 0,file_size[i] + 1)) {
                ret_code = -1;
                goto loop_error;
            }
            sprintf(file_name,"./data/%08db.png",p[i]);
            if(get_filedata(file_name,data[i] + 1,file_size[i] + 2,file_size[i] + 3)) {
                ret_code = -1;
                goto loop_error;
            }
            if(file_size[i][0] != file_size[i][2] || file_size[i][1] != file_size[i][3]) {
                ret_code = -1;
                goto loop_error;
            }
            if(file_size[i][0] == 0 || file_size[i][1] == 0) {
                ret_code = -1;
                goto loop_error;
            }
        }
        for(i=0;i<num_parent;i++) {
            sprintf(file_name,"./gene/gene%08d.bin",gene_pool_idx[i]);
            load_gene(file_name,parent_gene[i],element_size);
            parent_gene_diff[i] = 0.0;
        }
        for(i=0;i<num_data_set;i++) {
            gene_calc(parent_gene,parent_gene_diff,data[i],file_size[i][0],file_size[i][1],num_parent);
        }
        R_sort(parent_gene_diff,parent_gene_idx,num_parent);
        gene_generator(parent_gene,child_gene,parent_gene_idx,parent_gene_diff);
        for(i=0;i<num_generate;i++) {
            child_gene_diff[i] = 0.0;
        }
        for(i=0;i<num_data_set;i++) {
            gene_calc(child_gene,child_gene_diff,data[i],file_size[i][0],file_size[i][1],num_generate);
        }
        R_sort(child_gene_diff,child_gene_idx,num_generate);
        f_diff = 0.0;
        for(i=0;i<num_data_set;i++) {
            g_diff = 0.0;
            for(j=0;j<file_size[i][0]*file_size[i][1]*matrix_size[layer-1][1];j++) {
                data[i][0][j] -= data[i][1][j];
                g_diff += data[i][0][j]*data[i][0][j];
            }
            g_diff /= (RealNumber)(file_size[i][0]*file_size[i][1]*matrix_size[layer-1][1]);
            f_diff += g_diff;
        }
        f_diff /= (RealNumber)num_data_set;
        f_diff *= 255.0*255.0;
        printf("%f\n",f_diff);
        for(i=0;i<num_parent;i++) {
            for(j=0;j<element_size;j++) {
                parent_gene[parent_gene_idx[i]][j] = child_gene[child_gene_idx[i]][j];
            }
            sprintf(file_name,"./gene/gene%08d.bin",gene_pool_idx[parent_gene_idx[i]]);
            save_gene(file_name,parent_gene[parent_gene_idx[i]],element_size);
            f_diff = parent_gene_diff[parent_gene_idx[i]] / ((RealNumber)num_data_set);
            f_diff *= 255.0*255.0;
            g_diff = child_gene_diff[child_gene_idx[i]]  / ((RealNumber)num_data_set);
            g_diff *= 255.0*255.0;
            printf("%f -> %f\n",f_diff,g_diff);
            parent_gene_diff[parent_gene_idx[i]] = child_gene_diff[child_gene_idx[i]];
        }
        R_sort(parent_gene_diff,parent_gene_idx,num_parent);
        f_diff = parent_gene_diff[parent_gene_idx[0]] / ((RealNumber)num_data_set);
        f_diff *= 255.0*255.0;
        printf("\ngene%08d.bin\n%f\n\n",gene_pool_idx[parent_gene_idx[0]],f_diff);
        if(read(0,file_name,256)>=0) {
            break;
        }
    }
    loop_error:
    for(i=0;i<num_data_set;i++) {
        if(data[i][0] != NULL) {
            free(data[i][0]);
        }
        if(data[i][1] != NULL) {
            free(data[i][1]);
        }
    }
    
    free(p);
    malloc_p_error:
    free(file_name);
    malloc_file_name_error:
    free(child_gene[0]);
    malloc_child_gene_error2:
    free(parent_gene[0]);
    malloc_parent_gene_error2:
    free(data[0]);
    malloc_data_error2:
    free(file_size[0]);
    malloc_file_size_error2:
    free(child_gene_idx);
    malloc_child_gene_idx_error:
    free(parent_gene_idx);
    malloc_parent_gene_idx_error:
    free(gene_pool_idx);
    malloc_gene_pool_idx_error:
    free(child_gene_diff);
    malloc_child_gene_diff_error:
    free(parent_gene_diff);
    malloc_parent_gene_diff_error:
    free(child_gene);
    malloc_child_gene_error1:
    free(parent_gene);
    malloc_parent_gene_error1:
    free(data);
    malloc_data_error1:
    free(file_size);
    malloc_file_size_error1:
    free(lp_urandom);
    urandom_failed2:
    close(urandom_fd);
    urandom_failed1:
    return ret_code;
}
