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
RealNumber **distance_table;

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
    int i,j,k,l,m,n,ret_code;
    RealNumber *world_points;
    int world_stars;
    int *index[2];
    RealNumber distance,t1,t2;
    ret_code = 0;
    
    world_points = NULL;
    if(get_world("world01.bin",&world_points,&world_stars))
    {
        ret_code = -1;
        goto error01;
    }
    index[0] = (int *)malloc(2 * world_stars * sizeof(int));
    if(index[0] == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    index[1] = index[0] + world_stars;
    distance_table = (RealNumber **)malloc(world_stars*sizeof(RealNumber *));
    if(distance_table == NULL)
    {
        ret_code = -1;
        goto error03;
    }
    distance_table[0] = (RealNumber *)malloc(world_stars*world_stars*sizeof(RealNumber));
    if(distance_table[0] == NULL)
    {
        ret_code = -1;
        goto error04;
    }
    for(i=1;i<world_stars;i++)
    {
        distance_table[i] = distance_table[i-1] + world_stars;
    }
    for(i=0;i<world_stars;i++)
    {
        for(j=0;j<world_stars;j++)
        {
            t1 = 0.0;
            for(k=0;k<Dimension;k++)
            {
                t2 = world_points[Dimension*i + k] - world_points[Dimension*j + k];
                t1 += t2 * t2;
            }
            distance_table[i][j] = sqrt(t1);
        }
    }
    
    for(i=0;i<world_stars;i++)
    {
        index[0][i] = i;
    }
    distance = 0.0;

    for(i=0;i<world_stars-1;i++)
    {
        k = i+1;
        t1 = distance_table[index[0][k]][index[0][i]];
        for(j=i+2;j<world_stars;j++)
        {
            t2 = distance_table[index[0][j]][index[0][i]];
            if(t1 > t2)
            {
                t1 = t2;
                k = j;
            }
        }
        l = index[0][i+1];
        index[0][i+1] = index[0][k];
        index[0][k] = l;
        distance += t1;
    }
    t1 = distance_table[index[0][0]][index[0][world_stars-1]];
    distance += t1;
    printf("%f\n",distance);
    for(i=0;i<world_stars;i++)
    {
        index[1][i] = index[0][i];
    }

    int q[4];
    int flag = 1;
    while(flag)
    {
        flag = 0;
        for(i=0;i<world_stars;i++)
        {
            for(j=i+2;j<world_stars;j++)
            {
                q[0] = i;
                q[1] = i+1;
                q[2] = j;
                q[3] = j+1;
                if(q[1] >= world_stars)
                {
                    q[1] = 0;
                }
                if(q[3] >= world_stars)
                {
                    q[3] = 0;
                }
                m = 0;
                t1 = distance_table[index[0][q[0]]][index[0][q[1]]];
                t1 += distance_table[index[0][q[2]]][index[0][q[3]]];
                for(k=1;k<2;k++)
                {
                    t2 = distance_table[index[0][q[0]]][index[0][q[2]]];
                    t2 += distance_table[index[0][q[1]]][index[0][q[3]]];
                    if(t1 > t2)
                    {
                        t1 = t2;
                        m = k;
                    }
                }
                if(m != 0)
                {
                    l = q[1] + q[2];
                    for(n=q[1];n<=q[2];n++)
                    {
                        index[1][l-n] = index[0][n];
                    }
                    for(n=q[1];n<=q[2];n++)
                    {
                        index[0][n] = index[1][n];
                    }
                    flag = 1;
                }
            }
        }
    }

    t1 = 0.0;
    for(i=0;i<world_stars-1;i++)
    {
        t1 += distance_table[index[0][i+1]][index[0][i]];
    }
    t1 += distance_table[index[0][0]][index[0][world_stars-1]];
    printf("%f\n",t1);

    free(distance_table[0]);
    error04:
    free(distance_table);
    error03:
    free(index[0]);
    error02:
    free(world_points);
    error01:
    return ret_code;
}
