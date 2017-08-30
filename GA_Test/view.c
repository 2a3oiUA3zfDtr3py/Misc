#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <MagickWand/MagickWand.h>
#ifdef DOUBLE
#define RealNumber double
#define RealNumber_MAX DBL_MAX
#define RealNumber_MIN DBL_MIN
#else
#define RealNumber float
#define RealNumber_MAX FLT_MAX
#define RealNumber_MIN FLT_MIN
#endif
#define aux_vector_size 1
#define layer 12
#define input_datasize 3
#define output_datasize 3
#define hidden_channels 32
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
    {hidden_channels,output_datasize,1}
};

int element_size;

int vector_x_matrix(RealNumber *a,int a_size,RealNumber **b,int b_row,int b_col,RealNumber *c,int c_size)
{
    int i,j;
    for(j=0;j<b_col;j++)
    {
        c[j] =0.0;
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

int calc(RealNumber *gene,RealNumber *data,int width,int height)
{
    int return_code;
    int i,j,k,l,m,n;
    int p[pix_ref][2];
    RealNumber ***a;
    RealNumber **b;
    RealNumber *c;
    return_code = 0;
    
    if(gene == NULL || data == NULL || width <= 0 || height <= 0)
    {
        return_code = -1;
        goto error01;
    }

    a = (RealNumber ***)malloc(layer*sizeof(RealNumber **));
    if(a == NULL)
    {
        return_code = -1;
        goto error01;
    }
    j=0;
    for(i=0;i<layer;i++)
    {
        j += pix_ref*matrix_size[i][0] + aux_vector_size;
    }
    a[0] = (RealNumber **)malloc(j*sizeof(RealNumber *));
    if(a[0] == NULL)
    {
        return_code = -1;
        goto error02;
    }
    for(i=1;i<layer;i++)
    {
        a[i] = a[i-1] + pix_ref*matrix_size[i-1][0] + aux_vector_size;
    }
    j = matrix_size[0][0];
    for(i=0;i<layer;i++)
    {
        j = (j<matrix_size[i][1])? matrix_size[i][1] : j;
    }
    b = (RealNumber **)malloc((layer+1)*sizeof(RealNumber *));
    if(b == NULL)
    {
        return_code = -1;
        goto error03;
    }
    b[0] = (RealNumber *)malloc(2*j*width*height*sizeof(RealNumber));
    if(b == NULL)
    {
        return_code = -1;
        goto error04;
    }
    b[1] = b[0] + width*height;
    for(i=2;i<=layer;i++)
    {
        b[i] = b[i&0x01];
    }
    c = (RealNumber *)malloc(2*pix_ref*j*sizeof(RealNumber));
    if(c == NULL)
    {
        return_code = -1;
        goto error05;
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
        b[0][i] = data[i];
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

    for(i=0;i<height*width*matrix_size[layer-1][1];i++)
    {
        data[i] = b[layer][i];
    }

    free(c);
    error05:
    free(b[0]);
    error04:
    free(b);
    error03:
    free(a[0]);
    error02:
    free(a);
    error01:
    return return_code;
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
        goto file_error;
    }
    if(fstat(fd,&file_stat))
    {
        ret_code = -1;
        close(fd);
        goto file_error;
    }
    if(file_stat.st_size != element_size*sizeof(RealNumber))
    {
        ret_code = -1;
        close(fd);
        goto file_error;
    }
    i=read(fd,gene,file_stat.st_size);
    close(fd);
    if(i != file_stat.st_size)
    {
        ret_code = -1;
        goto file_error;
    }
    file_error:
    return ret_code;
}

int main(int argc,char **argv)
{
    int i,return_code;
    unsigned long width,height;
    RealNumber *data,*gene;
    MagickWand *img_handle;
    MagickBooleanType status;
    return_code = 0;
    if(argc != 3)
    {
        return_code = -1;
        goto error01;
    }
    
    element_size = 0;
    for(i=0;i<layer;i++)
    {
        element_size += (pix_ref*matrix_size[i][0] + aux_vector_size) * matrix_size[i][1];
    }
    gene = (RealNumber *)malloc(element_size*sizeof(RealNumber));
    if(gene == NULL)
    {
        return_code = -1;
        goto error01;
    }
    if(load_gene("gene/gene00000000.bin",gene))
    {
        return_code = -1;
        goto error02;
    }
    
    MagickWandGenesis();
    img_handle = NewMagickWand();
    if(img_handle == NULL)
    {
        return_code = -1;
        goto error03;
    }
    if(MagickReadImage(img_handle,argv[1]) == MagickFalse)
    {
        return_code = -1;
        goto error04;
    }
    width = MagickGetImageWidth(img_handle);
    height = MagickGetImageHeight(img_handle);
    if(width == 0 || height == 0)
    {
        return_code = -1;
        goto error04;
    }
    data = (RealNumber *)malloc(3*width*height*sizeof(RealNumber));
    if(data == NULL)
    {
        return_code = -1;
        goto error04;
    }
#ifdef DOUBLE
    status = MagickExportImagePixels(img_handle,0,0,width,height,"RGB",DoublePixel,data);
#else
    status = MagickExportImagePixels(img_handle,0,0,width,height,"RGB",FloatPixel,data);
#endif
    if(status == MagickFalse)
    {
        return_code = -1;
        goto error05;
    }/*
    if(calc(gene,data,(int)width,(int)height))
    {
        return_code = -1;
        goto error05;
    }*/
    i=0;
    while(argv[2][i] != 0)
    {
        i++;
    }
    if( (argv[2][i-3] == 'P' || argv[2][i-3] == 'p') && 
        (argv[2][i-2] == 'N' || argv[2][i-2] == 'n') && 
        (argv[2][i-1] == 'G' || argv[2][i-1] == 'g') )
    {
        if(MagickSetFormat(img_handle,"PNG64") == MagickFalse)
        {
            return_code = -1;
            goto error05;
        }
    }
#ifdef DOUBLE
    status = MagickImportImagePixels(img_handle,0,0,width,height,"RGB",DoublePixel,data);
#else
    status = MagickImportImagePixels(img_handle,0,0,width,height,"RGB",FloatPixel,data);
#endif
    if(status == MagickFalse)
    {
        return_code = -1;
        goto error05;
    }
    
    if(MagickWriteImage(img_handle,argv[2]) == MagickFalse)
    {
        return_code = -1;
        goto error05;
    }
    
    error05:
    free(data);
    error04:
    DestroyMagickWand(img_handle);
    error03:
    MagickWandTerminus();
    error02:
    free(gene);
    error01:
    return return_code;
}

