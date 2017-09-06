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
#define Dimension 2

int get_world(char *file_name,RealNumber **points,int *stars)
{
    int i,fd,ret_code;
    struct stat file_stats;
    ret_code = 0;
    
    if(file_name == NULL || stars == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    if(*points != NULL)
    {
        free(NULL);
    }
    *points = NULL;
    *stars = 0;
    
    fd = open(file_name,O_RDONLY);
    if(fd <= 2)
    {
        ret_code = -1;
        goto error01;
    }
    if(fstat(fd,&file_stats))
    {
        ret_code = -1;
        goto error02;
    }
    if(file_stats.st_size % (Dimension * sizeof(RealNumber)) != 0 || file_stats.st_size == 0)
    {
        ret_code = -1;
        goto error02;
    }
    *stars = file_stats.st_size / (Dimension * sizeof(RealNumber));
    *points = (RealNumber *)malloc(file_stats.st_size);
    if(*points == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    i = read(fd,*points,file_stats.st_size);
    if(i != file_stats.st_size)
    {
        ret_code = -1;
        goto error02;
    }
    error02:
    close(fd);
    
    if(ret_code)
    {
        if(*points != NULL)
        {
            free(*points);
        }
        *points = NULL;
    }
    error01:
    return ret_code;
}

int main()
{
    int i,j,k,l,ret_code;
    RealNumber *world_points;
    int world_stars;
    int *index;
    RealNumber distance,t1,t2,t3;
    ret_code = 0;
    
    world_points = NULL;
    if(get_world("world01.bin",&world_points,&world_stars))
    {
        ret_code = -1;
        goto error01;
    }
    index = (int *)malloc(world_stars * sizeof(int));
    if(index == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    for(i=0;i<world_stars;i++)
    {
        index[i] = i;
    }
    
    distance = 0.0;

    for(i=0;i<world_stars-1;i++)
    {
        k = i+1;
        t1 = 0.0;
        for(l=0;l<Dimension;l++)
        {
            t3 = world_points[Dimension*index[k] + l] - world_points[Dimension*index[i] + l];
            t1 += t3 * t3;
        }
        t1 = sqrt(t1);
        for(j=i+2;j<world_stars;j++)
        {
            t2 = 0.0;
            for(l=0;l<Dimension;l++)
            {
                t3 = world_points[Dimension*index[j] + l] - world_points[Dimension*index[i] + l];
                t2 += t3 * t3;
            }
            t2 = sqrt(t2);
            if(t1 > t2)
            {
                t1 = t2;
                k = j;
            }
        }
        l = index[i+1];
        index[i+1] = index[k];
        index[k] = l;
        distance += t1;
    }
    t1 = 0.0;
    for(l=0;l<Dimension;l++)
    {
        t3 = world_points[Dimension*index[0] + l] - world_points[Dimension*index[world_stars-1] + l];
        t1 += t3 * t3;
    }
    t1 = sqrt(t1);
    distance += t1;
    printf("%f\n",distance);

    for(i=0;i<world_stars;i++)
    {
        printf("%d\n",index[i]);
    }
    distance = 0.0;
    for(i=0;i<world_stars-1;i++)
    {
        t2 = 0.0;
        for(l=0;l<Dimension;l++)
        {
            t3 = world_points[Dimension*index[i+1] + l] - world_points[Dimension*index[i] + l];
            t2 += t3 * t3;
        }
        t2 = sqrt(t2);
        distance += t2;
    }
    t2 = 0.0;
    for(l=0;l<Dimension;l++)
    {
        t3 = world_points[Dimension*index[0] + l] - world_points[Dimension*index[world_stars-1] + l];
        t2 += t3 * t3;
    }
    t2 = sqrt(t2);
    distance += t2;
    printf("%f\n",distance);
    

    free(index);
    error02:
    free(world_points);
    error01:
    return ret_code;
}
