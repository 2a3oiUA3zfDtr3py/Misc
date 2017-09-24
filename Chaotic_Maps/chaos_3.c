#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <MagickWand/MagickWand.h>
#define width 800
#define height 600

void chaotic_map(unsigned int *v,double k)
{
    int i,j;
    double a[4],b[4],c[4],d;
    
    for(i=0;i<height;i++) {
        v[i] = 0;
    }
    
    a[0] = 0.320103;
    a[1] = 0.420103;
    a[2] = 0.5;
    a[3] = 0.5;
 
    for(i=0;i<100;i++)
    {
        b[0] = 1.0 - 2.0 * a[0];
        if(b[0] < 0.0)
        {
            b[0] = -b[0];
        }
        b[1] = b[0] * b[0];
        b[2] = b[1] * b[0];
        b[3] = b[1] * b[1];
        b[0] = 1.0 - b[0];
        b[1] = 1.0 - b[1];
        b[2] = 1.0 - b[2];
        b[3] = 1.0 - b[3];
        c[0] = a[1]*a[2]*b[0] + a[1]*(1-a[2])*b[1] + (1.0-a[1])*a[3]*b[2] + (1.0-a[1])*(1.0-a[3])*b[3];
        c[1] = a[1]*(1-a[2])*b[0] + (1.0-a[1])*a[3]*b[1] + (1.0-a[1])*(1.0-a[3])*b[2] + a[1]*a[2]*b[3];
        c[2] = (1.0-a[1])*a[3]*b[0] + (1.0-a[1])*(1.0-a[3])*b[1] + a[1]*a[2]*b[2] + a[1]*(1-a[2])*b[3];
        c[3] = (1.0-a[1])*(1.0-a[3])*b[0] + a[1]*a[2]*b[1] + a[1]*(1-a[2])*b[2] + (1.0-a[1])*a[3]*b[3];
        a[0] = k*c[0];
        a[1] = c[1];
        a[2] = c[2];
        a[3] = c[3];
    }
    for(i=0;i<10000;i++) {
        b[0] = 1.0 - 2.0 * a[0];
        if(b[0] < 0.0)
        {
            b[0] = -b[0];
        }
        b[1] = b[0] * b[0];
        b[2] = b[1] * b[0];
        b[3] = b[1] * b[1];
        b[0] = 1.0 - b[0];
        b[1] = 1.0 - b[1];
        b[2] = 1.0 - b[2];
        b[3] = 1.0 - b[3];
        c[0] = a[1]*a[2]*b[0] + a[1]*(1-a[2])*b[1] + (1.0-a[1])*a[3]*b[2] + (1.0-a[1])*(1.0-a[3])*b[3];
        c[1] = a[1]*(1-a[2])*b[0] + (1.0-a[1])*a[3]*b[1] + (1.0-a[1])*(1.0-a[3])*b[2] + a[1]*a[2]*b[3];
        c[2] = (1.0-a[1])*a[3]*b[0] + (1.0-a[1])*(1.0-a[3])*b[1] + a[1]*a[2]*b[2] + a[1]*(1-a[2])*b[3];
        c[3] = (1.0-a[1])*(1.0-a[3])*b[0] + a[1]*a[2]*b[1] + a[1]*(1-a[2])*b[2] + (1.0-a[1])*a[3]*b[3];
        a[0] = k*c[0];
        a[1] = c[1];
        a[2] = c[2];
        a[3] = c[3];
        d = a[0] * height;
        j = (int)d;
        if(j < 0) {
            j = 0;
        }
        else if(j >= height) {
            j = height -1;
        }
        v[j]++;
    }
    a[0] = 0.340103;
    a[1] = 0.420103;
    a[2] = 0.5;
    a[3] = 0.5;
    for(i=0;i<100;i++)
    {
        b[0] = 1.0 - 2.0 * a[0];
        if(b[0] < 0.0)
        {
            b[0] = -b[0];
        }
        b[1] = b[0] * b[0];
        b[2] = b[1] * b[0];
        b[3] = b[1] * b[1];
        b[0] = 1.0 - b[0];
        b[1] = 1.0 - b[1];
        b[2] = 1.0 - b[2];
        b[3] = 1.0 - b[3];
        c[0] = a[1]*a[2]*b[0] + a[1]*(1-a[2])*b[1] + (1.0-a[1])*a[3]*b[2] + (1.0-a[1])*(1.0-a[3])*b[3];
        c[1] = a[1]*(1-a[2])*b[0] + (1.0-a[1])*a[3]*b[1] + (1.0-a[1])*(1.0-a[3])*b[2] + a[1]*a[2]*b[3];
        c[2] = (1.0-a[1])*a[3]*b[0] + (1.0-a[1])*(1.0-a[3])*b[1] + a[1]*a[2]*b[2] + a[1]*(1-a[2])*b[3];
        c[3] = (1.0-a[1])*(1.0-a[3])*b[0] + a[1]*a[2]*b[1] + a[1]*(1-a[2])*b[2] + (1.0-a[1])*a[3]*b[3];
        a[0] = k*c[0];
        a[1] = c[1];
        a[2] = c[2];
        a[3] = c[3];
    }
    for(i=0;i<10000;i++) {
        b[0] = 1.0 - 2.0 * a[0];
        if(b[0] < 0.0)
        {
            b[0] = -b[0];
        }
        b[1] = b[0] * b[0];
        b[2] = b[1] * b[0];
        b[3] = b[1] * b[1];
        b[0] = 1.0 - b[0];
        b[1] = 1.0 - b[1];
        b[2] = 1.0 - b[2];
        b[3] = 1.0 - b[3];
        c[0] = a[1]*a[2]*b[0] + a[1]*(1-a[2])*b[1] + (1.0-a[1])*a[3]*b[2] + (1.0-a[1])*(1.0-a[3])*b[3];
        c[1] = a[1]*(1-a[2])*b[0] + (1.0-a[1])*a[3]*b[1] + (1.0-a[1])*(1.0-a[3])*b[2] + a[1]*a[2]*b[3];
        c[2] = (1.0-a[1])*a[3]*b[0] + (1.0-a[1])*(1.0-a[3])*b[1] + a[1]*a[2]*b[2] + a[1]*(1-a[2])*b[3];
        c[3] = (1.0-a[1])*(1.0-a[3])*b[0] + a[1]*a[2]*b[1] + a[1]*(1-a[2])*b[2] + (1.0-a[1])*a[3]*b[3];
        a[0] = k*c[0];
        a[1] = c[1];
        a[2] = c[2];
        a[3] = c[3];
        d = a[0] * height;
        j = (int)d;
        if(j < 0) {
            j = 0;
        }
        else if(j >= height) {
            j = height -1;
        }
        v[j]++;
    }
    
    return;
}

void sort(unsigned int *data,unsigned *s)
{
    int i,j,k,l;
    unsigned int *index[2],*radix[4];
    index[0] = s;
    index[1] = index[0] + height;
    radix[0] = index[1] + height;
    for(i=1;i<sizeof(unsigned int);i++) {
        radix[i] = radix[i-1] + 257;
    }
    for(i=0;i<height;i++) {
        index[0][i] = i;
    }
    for(i=0;i<257*sizeof(unsigned int);i++) {
        radix[0][i] = 0;
    }
    for(i=0;i<height;i++) {
        for(j=0;j<sizeof(unsigned int);j++) {
            k = *(((unsigned char *)(data+i)) + j);
            radix[j][k+1]++;
        }
    }
    for(i=0;i<sizeof(unsigned int);i++) {
        for(j=0;j<256;j++) {
            radix[i][j+1] += radix[i][j];
        }
    }
    l = 0;
    for(i=0;i<sizeof(unsigned int);i++) {
        for(j=0;j<height;j++) {
            k = *(((unsigned char *)(data+index[l][j])) + i);
            index[l^0x01][radix[i][k]++] = index[l][j];
        }
        l ^= 0x01;
    }
    if(l == 0x01) {
        for(i=0;i<height;i++) {
            index[l^0x01][j] = index[l][j];
        }
    }
    return;
}

int main()
{
    int i,j,k,ret_code;
    unsigned int *v,*s;
    double **img,temp;
    char filename[64];
    MagickWand *img_handle;
    PixelWand *pixel;
    ret_code = 0;
    
    i = height*3 + 257*sizeof(unsigned int);
    v = (unsigned int *)malloc(i*sizeof(unsigned int));
    if(v == NULL) {
        ret_code = -1;
        goto error01;
    }
    s = v + height;
    img = (double **)malloc(height*sizeof(double *));
    if(img == NULL) {
        ret_code = -1;
        goto error02;
    }
    img[0] = (double *)malloc(3*height*width*sizeof(double));
    if(img[0] == NULL) {
        ret_code = -1;
        goto error03;
    }
    for(i=1;i<height;i++) {
        img[i] = img[i-1] + 3*width;
    }
    
    for(i=0;i<width;i++)
    {
        temp = ((double)i)/width;
        chaotic_map(v,temp);
        sort(v,s);
        k = v[s[height-1]];
        for(j=0;j<height;j++)
        {
            if(k > v[height - j - 1]) {
                temp = ((double)v[height - j - 1])/((double)k);
            }
            else if(v[height - j - 1] == 0) {
                temp = 0.0;
            }
            else {
                temp = 1.0;
            }
            img[j][3*i + 0] = 0.0;
            img[j][3*i + 1] = temp;
            img[j][3*i + 2] = temp;
        }
    }
    
    MagickWandGenesis();
    img_handle = NewMagickWand();
    if(IsMagickWand(img_handle) == MagickFalse) {
        ret_code = -1;
        goto error04;
    }
    pixel = NewPixelWand();
    if(IsPixelWand(pixel) == MagickFalse) {
        ret_code = -1;
        goto error05;
    }
    if(PixelSetColor(pixel,"#000000") == MagickFalse) {
        ret_code = -1;
        goto error06;
    }
    if(MagickNewImage(img_handle,width,height,pixel) == MagickFalse) {
        ret_code = -1;
        goto error06;
    }
    if(MagickImportImagePixels(img_handle,0,0,width,height,"RGB",DoublePixel,img[0]) == MagickFalse) {
        ret_code = -1;
        goto error06;
    }
    sprintf(filename,"chaostic_3_map.png");
    if(MagickWriteImage(img_handle,filename) == MagickFalse) {
        ret_code = -1;
        goto error06;
    }
    
    
    
    error06:
    DestroyPixelWand(pixel);
    error05:
    DestroyMagickWand(img_handle);
    error04:
    MagickWandTerminus();
    free(img[0]);
    error03:
    free(img);
    error02:
    free(v);
    error01:
    return ret_code;
}
