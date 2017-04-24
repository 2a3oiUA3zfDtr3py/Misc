#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#define RealNum double

int main(void) {
    int i,j,fd;
    volatile int k;
    RealNum *b,c[2];
    unsigned char *temp;
    struct stat file_stat;
    c[0]=1.0;
    c[1]=-1.0;
    fd = open("temp1.bin",O_RDONLY);
    fstat(fd,&file_stat);
    i = file_stat.st_size / sizeof(RealNum);
    b = (RealNum *)malloc(i*sizeof(RealNum));
    read(fd,b,i*sizeof(RealNum));
    close(fd);
    for(k=0;k<100;k++) {
        for(j=0;j<i;j++) {
#if defined type1
            *(sizeof(RealNum) -1 + (unsigned char *)(b + j)) &= 0x7f;
#elif defined type2
            if(*(sizeof(RealNum) -1 + (unsigned char *)(b + j))&0x7f) {
                b[j] = -b[j];
            }
#elif defined type3
            temp = sizeof(RealNum) -1 + (unsigned char *)(b + j);
            if((*temp)>>7) {
                *temp &= 0x7f;
            }
#elif defined type4
            temp = sizeof(RealNum) -1 + (unsigned char *)(b + j);
            if(*temp&0x7f) {
                *temp &= 0x7f;
            }
#elif defined type5
            temp = sizeof(RealNum) -1 + (unsigned char *)(b + j);
            *temp&=0x7f;
#elif defined type6
            if(*(sizeof(RealNum) -1 + (unsigned char *)(b + j))>>7) {
                b[j] = -b[j];
            }
#elif defined type7
            b[j] *= c[(*(sizeof(RealNum) -1 + (unsigned char *)(b + j))>>7)];
#elif defined type8
            if(b[j] < 0.0) {
                b[j] = -b[j];
            }
#endif
        }
    }
    fd = open("temp2.bin",O_WRONLY|O_CREAT|O_TRUNC,0666);
    write(fd,b,i*sizeof(RealNum));
    close(fd);    
    return 0;
}
