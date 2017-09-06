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
int urandom_fd;
uint64_t *lp_urandom;
RealNumber *world_points;
int world_stars;
int num_parent;
int mutation_timer;
double *entropy_table;
RealNumber **distance_table;

int get_random(uint64_t *data,int size)
{
    static int i=UrandomBuff;
    int j,k,ret_code;
    ret_code = 0;
    if(lp_urandom == NULL || urandom_fd < 0)
    {
        ret_code = -1;
        goto error01;
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
                goto error01;
            }
            i=0;
        }
        data[k++] = lp_urandom[i++];
    }
    error01:
    return ret_code;
}

int I_sort_p(uint64_t *r,int *idx_a,int *idx_b,int *s,int size)
{
    int i,j,k,l,m,n,m_max,n_max,ret_code;
    int *index[2];
    ret_code = 0;
    if(r == NULL || idx_a == NULL || idx_b == NULL || idx_a == idx_b || size <= 0)
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
    if(s == NULL)
    {
        if(l == 1)
        {
            for(i=0;i<size;i++)
            {
                index[0][i] = index[1][i];
            }
        }
    }
    else
    {
        *s = l;
    }
    error:
    return ret_code;
}

int roulette_wheel_selector(uint64_t *urnd,int *i_temp_a,int *i_temp_b,int o_size,int size)
{
    int i,j,ret_code;
    ret_code = 0;
    if(urnd == NULL || i_temp_a == NULL || i_temp_b == NULL || i_temp_a == i_temp_b ||o_size < 0 || size < o_size || size < 0)
    {
        ret_code = -1;
        goto error01;
    }
    i_temp_a[0] = 0;
    i_temp_b[0] = size;
    for(i=1;i<size;i++)
    {
        i_temp_a[i] = i;
        i_temp_b[i] = i_temp_b[i-1] + size - i;
    }
    get_random(urnd,o_size);
    for(i=0;i<o_size;i++)
    {
        urnd[i] %= i_temp_b[size-i-1];
        for(j=0;j<size-i;j++)
        {
            if(urnd[i] < i_temp_b[j])
            {
                break;
            }
        }
        urnd[i] = i_temp_a[j];
        for(;j<size-i-1;j++)
        {
            i_temp_a[j] = i_temp_a[j+1];
            i_temp_b[j] = i_temp_b[j+1] - size + urnd[i];
        }
    }
    error01:
    return ret_code;
}

int evaluate_score(int **gene,RealNumber *score)
{
    int i,j,ret_code;
    ret_code = 0;
    if(gene == NULL || score == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    for(i=0;i<num_parent;i++)
    {
        score[i] = 0.0;
        for(j=0;j<world_stars-1;j++)
        {
            score[i] += distance_table[gene[i][j+1]][gene[i][j]];
        }
        score[i] += distance_table[gene[i][0]][gene[i][world_stars-1]];
    }
    error01:
    return ret_code;
}

int evaluate_entropy(int **gene,RealNumber *entropy)
{
    int i,j,ret_code;
    int **count;
    double temp;
    ret_code = 0;
    
    count = (int **)malloc(world_stars * sizeof(int *));
    if(count == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    count[0] = (int *)malloc(world_stars * world_stars * sizeof(int));
    if(count[0] == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    for(i=1;i<world_stars;i++)
    {
        count[i] = count[i-1] + world_stars;
    }
    for(i=0;i<world_stars*world_stars;i++)
    {
        count[0][i] = 0;
    }
    for(i=0;i<num_parent;i++)
    {
        for(j=0;j<world_stars-1;j++)
        {
            count[gene[i][j]][gene[i][j+1]]++;
            count[gene[i][j+1]][gene[i][j]]++;
        }
        count[gene[i][0]][gene[i][world_stars-1]]++;
        count[gene[i][world_stars-1]][gene[i][0]]++;
    }
    for(i=0;i<num_parent;i++)
    {
        temp = 0.0;
        for(j=0;j<world_stars-1;j++)
        {
            temp += entropy_table[count[gene[i][j+1]][gene[i][j]]];
        }
        temp += entropy_table[count[gene[i][0]][gene[i][world_stars-1]]];
        entropy[i] = (RealNumber)temp;
    }
    free(count[0]);
    error02:
    free(count);
    error01:
    return ret_code;
}

int gene_regulator(int *gene,int *temp_gene)
{
    int i,j,k,ret_code;
    if(gene == NULL || temp_gene == NULL || gene == temp_gene)
    {
        ret_code = -1;
        goto error01;
    }
    for(i=0;i<world_stars;i++)
    {
        temp_gene[i] = gene[i];
    }
    j = 0;
    k = 0;
    for(i=1;i<world_stars;i++)
    {
        j += i * temp_gene[i];
        k += (world_stars-i) * temp_gene[i];
    }
    if(j < k)
    {
        for(i=1;i<world_stars;i++)
        {
            gene[i] = temp_gene[i];
        }
    }
    else
    {
        for(i=1;i<world_stars;i++)
        {
            gene[world_stars-i] = temp_gene[i];
        }
    }
    error01:
    return ret_code;
}

int genesis(int **gene)
{
    int i,j,k,ret_code;
    int *temp[2];
    uint64_t *urnd;
    ret_code = 0;
    if(gene == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    temp[0] = (int *)malloc(2 * world_stars*sizeof(int));
    if(temp[0] == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    temp[1] = temp[0] + world_stars;
    urnd = (uint64_t *)malloc(world_stars*sizeof(uint64_t));
    if(urnd == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    
    for(i=0;i<num_parent;i++)
    {
        get_random(urnd+1,world_stars-1);
        urnd[0] = 0;
        I_sort_p(urnd,temp[0],temp[1],&j,world_stars);
        gene_regulator(temp[j],temp[j^0x01]);
        for(k=0;k<world_stars;k++)
        {
            gene[i][k] = temp[j][k];
        }
    }
    free(urnd);
    error02:
    free(temp[0]);
    error01:
    return ret_code;
}

int sort(int *index,RealNumber *diff)
{
    int i,j,k,l,m,n,m_max,n_max,ret_code;
    int *idx[2];
    if(index == NULL || diff == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    idx[1] = (int *)malloc(num_parent * sizeof(int));
    if(idx[1] == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    idx[0] = index;

    for(i=0;i<num_parent;i++)
    {
        idx[0x00][i] = i;
    }
    l = 0x00;
    for(i=1;i<num_parent;i*=2)
    {
        for(j=0;j<num_parent;j+=i*2)
        {
            k=j;
            m=j;
            n=j+i;
            m_max = (m+i < num_parent)? m+i : num_parent;
            n_max = (n+i < num_parent)? n+i : num_parent;
            while(m<m_max && n<n_max)
            {
                if(diff[idx[l][m]] < diff[idx[l][n]])
                {
                    idx[l^0x01][k++] = idx[l][m++];
                }
                else
                {
                    idx[l^0x01][k++] = idx[l][n++];
                }
            }
            while(m<m_max)
            {
                idx[l^0x01][k++] = idx[l][m++];
            }
            while(n<n_max)
            {
                idx[l^0x01][k++] = idx[l][n++];
            }
        }
        l ^= 0x01;
    }
    if(l == 0x01)
    {
        for(i=0;i<num_parent;i++)
        {
            idx[0][i] = idx[1][i];
        }
    }
    free(idx[1]);
    error01:
    return ret_code;
}

int progress_generation(int **gene,int *index)
{
    int i,j,k,ret_code;
    uint64_t *urnd;
    int **temp_gene;
    int *temp_list[3];
    ret_code = 0;

    urnd = (uint64_t *)malloc((world_stars + 5) * sizeof(uint64_t));
    if(urnd == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    temp_gene = (int **)malloc(num_parent * sizeof(int *));
    if(temp_gene == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    temp_gene[0] = (int *)malloc(num_parent * world_stars * sizeof(int));
    if(temp_gene[0] == NULL)
    {
        ret_code = -1;
        goto error03;
    }
    for(i=1;i<num_parent;i++)
    {
        temp_gene[i] = temp_gene[i-1] + world_stars;
    }
    i = (num_parent > world_stars)? num_parent : world_stars;
    i = (i > 5)? i : 5;
    temp_list[0] = (int *)malloc(3 * i * sizeof(int));
    if(temp_list[0] == NULL)
    {
        ret_code = -1;
        goto error04;
    }
    temp_list[1] = temp_list[0] + i;
    temp_list[2] = temp_list[1] + i;

    for(i=0;i<5;i++)
    {
        for(j=0;j<world_stars;j++)
        {
            temp_gene[i][j] = gene[index[i]][j];
        }
    }
    for(i=5;i<num_parent;i++)
    {
        mutation_timer++;
        if(mutation_timer == 10)
        {
            mutation_timer = 0;
            roulette_wheel_selector(urnd,temp_list[0],temp_list[1],1,num_parent);
            for(j=0;j<world_stars;j++)
            {
                temp_gene[i][j] = gene[index[urnd[0]]][j];
            }
            get_random(urnd,2);
            urnd[0] %= world_stars - 1;
            urnd[1] %= world_stars - 2;
            if(urnd[0] > urnd[1])
            {
                urnd[2] = urnd[0];
                urnd[0] = urnd[1];
                urnd[1] = urnd[2];
            }
            else if(urnd[0] == urnd[1])
            {
                urnd[1]++;
            }
            urnd[0]++;
            urnd[1]++;
            urnd[2] = urnd[0] + urnd[1];
            for(j=0;j<world_stars;j++)
            {
                if(urnd[0] <= temp_gene[i][j] && temp_gene[i][j] <= urnd[1])
                {
                    temp_gene[i][j] = urnd[2] - temp_gene[i][j];
                }
            }
        }
        else
        {
            roulette_wheel_selector(urnd,temp_list[0],temp_list[1],2,num_parent/2);
            get_random(urnd+2,world_stars);
            for(j=0;j<world_stars;j++)
            {
                temp_list[0][j] = 0;
                temp_list[1][j] = 0;
            }
            for(j=0;j<world_stars;j++)
            {
                if(urnd[3+j] % 2 == 0)
                {
                    temp_gene[i][j] = gene[index[urnd[0]]][j];
                    temp_list[0][gene[index[urnd[0]]][j]] = 1;
                    temp_list[1][j] = 1;
                }
            }
            k = 0;
            while(temp_list[1][k])
            {
                k++;
            }
            j = 0;
            while(j<world_stars)
            {
                if(temp_list[0][gene[index[urnd[1]]][j]] == 0)
                {
                    temp_gene[i][k] = gene[index[urnd[1]]][j];
                    temp_list[0][gene[index[urnd[1]]][j]] = 1;
                    temp_list[1][k] = 1;
                    while(temp_list[1][k])
                    {
                        k++;
                    }
                }
                j++;
            }
        }
    }

    for(i=0;i<num_parent;i++)
    {
        for(j=0;j<world_stars;j++)
        {
            gene[i][j] = temp_gene[i][j];
        }
        gene_regulator(gene[i],temp_gene[i]);
    }

    free(temp_list[0]);
    error04:
    free(temp_gene[0]);
    error03:
    free(temp_gene);
    error02:
    free(urnd);
    error01:
    return ret_code;
}

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
    int **gene,*index;
    RealNumber *diff,*entropy;
    double t1,t2;
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
    world_points = NULL;
    world_stars = 0;
    if(get_world("world01.bin",&world_points,&world_stars))
    {
        ret_code = -1;
        goto error03;
    }
    num_parent = 3 * world_stars;
    num_parent = (num_parent > 100)? num_parent : 100;
    mutation_timer = 0;
    gene = (int **)malloc((num_parent + 1) * sizeof(int *));
    if(gene == NULL)
    {
        ret_code = -1;
        goto error04;
    }
    gene[0] = (int *)malloc((num_parent + 1) * world_stars * sizeof(int));
    if(gene[0] == NULL)
    {
        ret_code = -1;
        goto error05;
    }    
    for(i=1;i<num_parent+1;i++)
    {
        gene[i] = gene[i-1] + world_stars;
    }
    i = (num_parent > world_stars)? num_parent : world_stars;
    index = (int *)malloc(i*sizeof(int));
    if(index == NULL)
    {
        ret_code = -1;
        goto error06;
    }
    diff = (RealNumber *)malloc(num_parent*sizeof(RealNumber));
    if(diff == NULL)
    {
        ret_code = -1;
        goto error07;
    }
    entropy = (RealNumber *)malloc(num_parent*sizeof(RealNumber));
    if(entropy == NULL)
    {
        ret_code = -1;
        goto error08;
    }
    entropy_table = (double *)malloc((num_parent+1)*sizeof(double));
    if(entropy_table == NULL)
    {
        ret_code = -1;
        goto error09;
    }
    distance_table = (RealNumber **)malloc(world_stars*sizeof(double));
    if(distance_table == NULL)
    {
        ret_code = -1;
        goto error10;
    }
    distance_table[0] = (RealNumber *)malloc(world_stars*world_stars*sizeof(double));
    if(distance_table[0] == NULL)
    {
        ret_code = -1;
        goto error11;
    }
    for(i=1;i<world_stars;i++)
    {
        distance_table[i] = distance_table[i-1] + world_stars;
    }

    if(world_stars <= 1)
    {
        goto too_few_input;
    }

    for(i=1;i<=num_parent;i++)
    {
        entropy_table[i] = -log(((double)i) / ((double)num_parent));
    }
    entropy_table[0] = entropy_table[1];
    
    for(i=0;i<world_stars;i++)
    {
        for(j=0;j<world_stars;j++)
        {
            t1 = 0.0;
            for(k=0;k<Dimension;k++)
            {
                t2 = world_points[Dimension*i + k] - world_points[Dimension*j + k];
                t1 += t2*t2;
            }
            distance_table[i][j] = sqrt(t1);
        }
    }
    
    
    for(i=0;i<world_stars;i++)
    {
        index[i] = i;
    }

    for(i=0;i<world_stars-1;i++)
    {
        k = i+1;
        t1 = distance_table[index[k]][index[i]];
        for(j=i+2;j<world_stars;j++)
        {
            t2 = distance_table[index[j]][index[i]];
            if(t1 > t2)
            {
                t1 = t2;
                k = j;
            }
        }
        l = index[i+1];
        index[i+1] = index[k];
        index[k] = l;
    }
    genesis(gene);
    for(i=0;i<world_stars;i++)
    {
        gene[0][i] = index[i];
    }
    gene_regulator(gene[0],gene[num_parent]);
    evaluate_score(gene,diff);
    if(world_stars <= 3)
    {
        for(i=0;i<world_stars;i++)
        {
            printf("%d\n",index[i]);
        }
        printf("%f\n",diff[0]);
        goto too_few_input;
    }
    sort(index,diff);
    t2 = 2.0*(diff[index[(3*num_parent)/4]] - diff[index[num_parent/4]])/entropy_table[1];
    while(t2*entropy_table[1] > 0.02)
    {
        evaluate_score(gene,diff);
        evaluate_entropy(gene,entropy);
        for(j=0;j<num_parent;j++)
        {
            diff[j] -= t2 * entropy[j];
        }
        sort(index,diff);
        progress_generation(gene,index);
        t2 *= 0.9994;
        printf("%f C\n",t2*entropy_table[1]);
    }
    evaluate_score(gene,diff);
    sort(index,diff);
    printf("%f\n",diff[index[0]]);
    
    too_few_input:
    free(distance_table[0]);
    error11:
    free(distance_table);
    error10:
    free(entropy_table);
    error09:
    free(entropy);
    error08:
    free(diff);
    error07:
    free(index);
    error06:
    free(gene[0]);
    error05:
    free(gene);
    error04:
    free(world_points);
    error03:
    free(lp_urandom);
    error02:
    close(urandom_fd);
    error01:
    return ret_code;
}
