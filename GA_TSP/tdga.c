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
int **probability_table;
RealNumber *entropy_table;
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

int U64_sort(uint64_t *figure,int *index_a,int *index_b,int *s,int size)
{
    int i,j,k,l,m,n,m_max,n_max,ret_code;
    int *index[2];
    ret_code = 0;
    if(figure == NULL || index_a == NULL || index_b == NULL || index_a == index_b || size <= 0)
    {
        ret_code = -1;
        goto error;
    }
    index[0] = index_a;
    index[1] = index_b;
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
                if(figure[index[l][m]] < figure[index[l][n]])
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

int I_sort(int *figure,int *index_a,int *index_b,int *s,int size)
{
    int i,j,k,l,m,n,m_max,n_max,ret_code;
    int *index[2];
    ret_code = 0;
    if(figure == NULL || index_a == NULL || index_b == NULL || index_a == index_b || size <= 0)
    {
        ret_code = -1;
        goto error;
    }
    index[0] = index_a;
    index[1] = index_b;
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
                if(figure[index[l][m]] < figure[index[l][n]])
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

int R_sort(RealNumber *figure,int *index_a,int *index_b,int *s,int size)
{
    int i,j,k,l,m,n,m_max,n_max,ret_code;
    int *index[2];
    ret_code = 0;
    if(figure == NULL || index_a == NULL || index_b == NULL || index_a == index_b || size <= 0)
    {
        ret_code = -1;
        goto error;
    }
    index[0] = index_a;
    index[1] = index_b;
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
                if(figure[index[l][m]] < figure[index[l][n]])
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

int initialize_distance_table(void)
{
    int i,j,k,ret_code;
    double t1,t2;
    ret_code = 0;
    if(distance_table == NULL || distance_table[0] == NULL)
    {
        ret_code = -1;
        goto error01;
    }

    for(i=0;i<world_stars;i++)
    {
        for(j=0;j<world_stars;j++)
        {
            t1 = 0.0;
            for(k=0;k<Dimension;k++)
            {
                t2 = (double)(world_points[Dimension*i + k]);
                t2 -= (double)(world_points[Dimension*j + k]);
                t1 += t2*t2;
            }
            distance_table[i][j] = (RealNumber)(sqrt(t1));
        }
    }

    error01:
    return ret_code;
}

int initialize_entropy_table(int num_parent)
{
    int i,ret_code;
    double *temp_v,temp_s;
    ret_code = 0;
    if(entropy_table == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    temp_v = (double *)malloc(num_parent*sizeof(double));
    if(temp_v == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    for(i=1;i<num_parent;i++)
    {
        temp_s = ((double)i) / ((double)num_parent);
        temp_s = temp_s * -log(temp_s);
        temp_v[i] = temp_s * ((double)i);
    }
    temp_v[0] = 0.0;
    temp_v[num_parent] = 0.0;
    for(i=0;i<num_parent;i++)
    {
        entropy_table[i] = temp_v[i];
    }
    entropy_table[num_parent] = 0.0;

    free(temp_v);
    error01:
    return ret_code;
}

int initialize_probability_table(void)
{
    int i,j,ret_code;
    ret_code = 0;
    if(probability_table == NULL || probability_table[0] == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    for(i=0;i<world_stars;i++)
    {
        for(j=0;j<world_stars;j++)
        {
            probability_table[i][j] = 0;
        }
    }
    error01:
    return ret_code;
}

int calc_probability_table(int **gene,int num_parent)
{
    int i,j,ret_code;
    ret_code = 0;
    if(probability_table == NULL || probability_table[0] == NULL || gene == NULL || gene[0] == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    initialize_probability_table();
    for(i=0;i<num_parent;i++)
    {
        for(j=0;j<world_stars-1;j++)
        {
            probability_table[gene[i][j]][gene[i][j+1]]++;
            probability_table[gene[i][j+1]][gene[i][j]]++;
        }
        probability_table[gene[i][0]][gene[i][world_stars-1]]++;
        probability_table[gene[i][world_stars-1]][gene[i][0]]++;
    }
    error01:
    return ret_code;
}

int remove_probability_table(int *gene)
{
    int i,ret_code;
    ret_code = 0;
    if(probability_table == NULL || probability_table[0] == NULL || gene  == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    
    for(i=0;i<world_stars-1;i++)
    {
        probability_table[gene[i]][gene[i+1]]--;
        probability_table[gene[i+1]][gene[i]]--;
    }
    probability_table[gene[0]][gene[world_stars-1]]--;
    probability_table[gene[world_stars-1]][gene[0]]--;
    
    error01:
    return ret_code;
}

int add_probability_table(int *gene)
{
    int i,ret_code;
    ret_code = 0;
    if(probability_table == NULL || probability_table[0] == NULL || gene  == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    
    for(i=0;i<world_stars-1;i++)
    {
        probability_table[gene[i]][gene[i+1]]++;
        probability_table[gene[i+1]][gene[i]]++;
    }
    probability_table[gene[0]][gene[world_stars-1]]++;
    probability_table[gene[world_stars-1]][gene[0]]++;

    error01:
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

int evaluate_distance(int **gene,RealNumber *distance,int size)
{
    int i,j,ret_code;
    ret_code = 0;
    if(gene == NULL || distance == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    for(i=0;i<size;i++)
    {
        distance[i] = 0.0;
        for(j=0;j<world_stars-1;j++)
        {
            distance[i] += distance_table[gene[i][j+1]][gene[i][j]];
        }
        distance[i] += distance_table[gene[i][0]][gene[i][world_stars-1]];
    }
    error01:
    return ret_code;
}

int evaluate_entropy(int **gene,double *entropy,int size)
{
    int i,j,k,ret_code;
    double *temp_v,temp_s;
    ret_code = 0;

    temp_v = (double *)malloc(size*sizeof(double));
    if(temp_v == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    for(i=0;i<size;i++)
    {
        temp_s = ((double)i) / ((double)size);
        temp_s = temp_s * -log(temp_s);
        temp_v[i] = temp_s;
    }

    temp_s = 0.0;
    for(i=0;i<size;i++)
    {
        for(j=0;j<world_stars-1;j++)
        {
            temp_s += temp_v[probability_table[gene[i][j+1]][gene[i][j]]];
        }
        temp_s += temp_v[probability_table[gene[i][0]][gene[i][world_stars-1]]];
    }
    
    *entropy = temp_s;

    free(temp_v);
    error01:
    return ret_code;
}

int evaluate_diff_entropy(int **gene,RealNumber *entropy,int size)
{
    int i,j,k,ret_code;
    RealNumber temp;
    ret_code = 0;

    for(i=0;i<size;i++)
    {
        temp = 0.0;
        for(j=0;j<world_stars-1;j++)
        {
            temp += entropy_table[probability_table[gene[i][j+1]][gene[i][j]]];
        }
        temp += entropy_table[probability_table[gene[i][0]][gene[i][world_stars-1]]];
        entropy[i] = temp;
    }
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

int genesis(int **gene,int size)
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
    
    for(i=0;i<size;i++)
    {
        get_random(urnd+1,world_stars-1);
        urnd[0] = 0;
        U64_sort(urnd,temp[0],temp[1],&j,world_stars);
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


int minimize_free_energy(int **child_gene,RealNumber *child_distance,int child_size,int **parent_gene,RealNumber *parent_distance,int parent_size,double Temperature)
{
    int i,j,k,l,m,n,ret_code;
    int *index[3],**temp_gene,temp_gene_size;
    uint64_t *route;
    RealNumber score_p,score_c,*temp_gene_distance;

    ret_code = 0;
    if(child_gene == NULL || child_gene[0] == NULL || child_distance == NULL || parent_gene == NULL || parent_gene[0] == NULL || parent_distance == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    temp_gene_size = child_size - parent_size;
    temp_gene = (int **)malloc(temp_gene_size*sizeof(int *));
    if(temp_gene == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    j = (2*world_stars > child_size)? 2*world_stars : child_size;
    temp_gene[0] = (int *)malloc((temp_gene_size*world_stars + 3*j)*sizeof(int));
    if(temp_gene[0] == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    for(i=1;i<temp_gene_size;i++)
    {
        temp_gene[i] = temp_gene[i-1] + world_stars;
    }
    index[0] = temp_gene[0] + temp_gene_size*world_stars;
    index[1] = index[0] + j;
    index[2] = index[1] + j;
    route = (uint64_t *)malloc(2*world_stars*sizeof(uint64_t));
    if(route == NULL)
    {
        ret_code = -1;
        goto error03;
    }
    temp_gene_distance = (RealNumber *)malloc(temp_gene_size*sizeof(RealNumber));
    if(temp_gene_distance == NULL)
    {
        ret_code = -1;
        goto error04;
    }

    initialize_probability_table();
    R_sort(child_distance,index[0],index[1],&k,child_size);
    for(i=0;i<parent_size;i++)
    {
        add_probability_table(child_gene[index[k][i]]);
        parent_distance[i] = child_distance[index[k][i]];
        for(j=0;j<world_stars;j++)
        {
            parent_gene[i][j] = child_gene[index[k][i]][j];
            child_gene[index[k][i]][j] = -1;
        }
    }
    k=0;
    for(i=0;i<child_size;i++)
    {
        if(child_gene[i][0] != -1)
        {
            for(j=0;j<world_stars;j++)
            {
                temp_gene[k][j] = child_gene[i][j];
            }
            temp_gene_distance[k] = child_distance[i];
            k++;
        }
    }
    for(i=0;i<parent_size;i++)
    {
        for(j=0;j<temp_gene_size;j++)
        {
            l=0;
            for(k=0;k<world_stars-1;k++)
            {
                if(parent_gene[i][k+1] < parent_gene[i][k])
                {
                    route[l++] = world_stars * parent_gene[i][k+1] + parent_gene[i][k];
                }
                else
                {
                    route[l++] = world_stars * parent_gene[i][k] + parent_gene[i][k+1];
                }
            }
            if(parent_gene[i][0] < parent_gene[i][world_stars-1])
            {
                route[l++] = world_stars * parent_gene[i][0] + parent_gene[i][world_stars-1];
            }
            else
            {
                route[l++] = world_stars * parent_gene[i][world_stars-1] + parent_gene[i][0];
            }
            for(k=0;k<world_stars-1;k++)
            {
                if(temp_gene[j][k+1] < temp_gene[j][k])
                {
                    route[l++] = world_stars * temp_gene[j][k+1] + temp_gene[j][k];
                }
                else
                {
                    route[l++] = world_stars * temp_gene[j][k] + temp_gene[j][k+1];
                }
            }
            if(temp_gene[j][0] < temp_gene[j][world_stars-1])
            {
                route[l++] = world_stars * temp_gene[j][0] + temp_gene[j][world_stars-1];
            }
            else
            {
                route[l++] = world_stars * temp_gene[j][world_stars-1] + temp_gene[j][0];
            }

            U64_sort(route,index[2],index[1],NULL,l);
            m = 0;
            index[0][m] = (int)(route[index[2][0]] / ((uint64_t)world_stars));
            index[1][m] = (int)(route[index[2][0]] % ((uint64_t)world_stars));
            m++;
            for(k=1;k<l;k++)
            {
                if(route[index[2][k]] != route[index[2][k-1]])
                {
                    index[0][m] = (int)(route[index[2][k]] / ((uint64_t)world_stars));
                    index[1][m] = (int)(route[index[2][k]] % ((uint64_t)world_stars));
                    m++;
                }
            }

            score_p = 0.0;
            for(k=0;k<m;k++)
            {
                score_p += entropy_table[probability_table[index[0][k]][index[1][k]]];
            }
            score_p *= -Temperature;
            score_p += parent_distance[i];
            remove_probability_table(parent_gene[i]);
            add_probability_table(temp_gene[j]);
            score_c = 0.0;
            for(k=0;k<m;k++)
            {
                score_c += entropy_table[probability_table[index[0][k]][index[1][k]]];
            }
            score_c *= -Temperature;
            score_c += temp_gene_distance[j];
            if(score_p > score_c)
            {
                score_c = temp_gene_distance[j];
                temp_gene_distance[j] = parent_distance[i];
                parent_distance[i] = score_c;
                for(k=0;k<world_stars;k++)
                {
                    l = parent_gene[i][k];
                    parent_gene[i][k] = temp_gene[j][k];
                    temp_gene[j][k] = l;
                }
            }
            else
            {
                remove_probability_table(temp_gene[j]);
                add_probability_table(parent_gene[i]);
            }
        }
    }
    
    free(temp_gene_distance);
    error04:
    free(route);
    error03:
    free(temp_gene[0]);
    error02:
    free(temp_gene);
    error01:
    return ret_code;
}

int child_generator(int **parent_gene,int *parent_index,int num_parent,int **child_gene,int num_child)
{
    int i,j,k,l,ret_code;
    uint64_t *urnd;
    int *temp[3];
    ret_code = 0;

    urnd = (uint64_t *)malloc((world_stars + 5) * sizeof(uint64_t));
    if(urnd == NULL)
    {
        ret_code = -1;
        goto error01;
    }
    i = (num_parent > world_stars)? num_parent : world_stars;
    i = (i > num_child)? i : num_child;
    i = (i > 5)? i : 5;
    temp[0] = (int *)malloc(3 * i * sizeof(int));
    if(temp[0] == NULL)
    {
        ret_code = -1;
        goto error02;
    }
    temp[1] = temp[0] + i;
    temp[2] = temp[1] + i;

    for(i=0;i<num_child;i++)
    {
        get_random(urnd,1);
        urnd[0] %= 1000;
        if(urnd[0] < 1)
        {
            roulette_wheel_selector(urnd,temp[0],temp[1],2,num_parent);
            get_random(urnd+2,1);
            for(j=0;j<world_stars;j++)
            {
                temp[0][j] = 0;
                temp[1][j] = 0;
            }
            urnd[2+j] %= world_stars - 2;
            urnd[2+j]++;
            for(j=0;j<urnd[2+j];j++)
            {
                child_gene[i][j] = parent_gene[parent_index[urnd[0]]][j];
                temp[0][parent_gene[parent_index[urnd[0]]][j]] = 1;
                temp[1][j] = 1;
            }
            k = 0;
            while(temp[1][k])
            {
                k++;
            }
            j = 0;
            while(j<world_stars)
            {
                if(temp[0][parent_gene[parent_index[urnd[1]]][j]] == 0)
                {
                    child_gene[i][k] = parent_gene[parent_index[urnd[1]]][j];
                    temp[0][parent_gene[parent_index[urnd[1]]][j]] = 1;
                    temp[1][k] = 1;
                    while(temp[1][k])
                    {
                        k++;
                    }
                }
                j++;
            }
        }
        else
        {
            roulette_wheel_selector(urnd,temp[0],temp[1],2,num_parent);
            get_random(urnd+2,world_stars);
            for(j=0;j<world_stars;j++)
            {
                temp[0][j] = 0;
                temp[1][j] = 0;
            }
            for(j=0;j<world_stars;j++)
            {
                if(urnd[2+j] % 2 == 0)
                {
                    child_gene[i][j] = parent_gene[parent_index[urnd[0]]][j];
                    temp[0][parent_gene[parent_index[urnd[0]]][j]] = 1;
                    temp[1][j] = 1;
                }
            }
            k = 0;
            while(temp[1][k])
            {
                k++;
            }
            j = 0;
            while(j<world_stars)
            {
                if(temp[0][parent_gene[parent_index[urnd[1]]][j]] == 0)
                {
                    child_gene[i][k] = parent_gene[parent_index[urnd[1]]][j];
                    temp[0][parent_gene[parent_index[urnd[1]]][j]] = 1;
                    temp[1][k] = 1;
                    while(temp[1][k])
                    {
                        k++;
                    }
                }
                j++;
            }
        }

        get_random(urnd,1);
        urnd[0] %= 1000;
        if(urnd[0] < 30)
        {
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
            if(urnd[0] == 1 && urnd[1] == world_stars - 1)
            {
                get_random(urnd+2,1);
                if(urnd[2] % 2 == 0)
                {
                    urnd[0]++;
                }
                else
                {
                    urnd[1]--;
                }
            }
            urnd[2] = urnd[0] + urnd[1];
            for(j=0;j<world_stars;j++)
            {
                if(urnd[0] <= child_gene[i][j] && child_gene[i][j] <= urnd[1])
                {
                    child_gene[i][j] = urnd[2] - child_gene[i][j];
                }
            }
        }
        else if(urnd[0] < 100)
        {
            roulette_wheel_selector(urnd,temp[0],temp[1],1,num_parent);
            for(j=0;j<world_stars;j++)
            {
                temp[2][j] = child_gene[i][j];
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
            if(urnd[0] == 1 && urnd[1] == world_stars - 1)
            {
                get_random(urnd+2,1);
                if(urnd[2] % 2 == 0)
                {
                    urnd[0]++;
                }
                else
                {
                    urnd[1]--;
                }
            }
            k = urnd[1] - urnd[0] + 1;
            if(k<=0)
            {
                k = 1;
            }
            get_random(urnd+2,k);
            U64_sort(urnd+2,temp[0],temp[1],NULL,k);
            k=0;
            for(j=0;j<world_stars;j++)
            {
                if(urnd[0] <= temp[2][j] && temp[2][j] <= urnd[1])
                {
                    child_gene[i][j] = urnd[0] + temp[0][k++];
                }
                else
                {
                    child_gene[i][j] = temp[2][j];
                }
            }
        }
    }

    for(i=0;i<num_child;i++)
    {
        gene_regulator(child_gene[i],temp[0]);
        for(j=0;j<world_stars;j++)
        {
            temp[0][j] = 0;
        }
        for(j=0;j<world_stars;j++)
        {
            temp[0][child_gene[i][j]] = 1;
        }
        for(j=0;j<world_stars;j++)
        {
            if(temp[0][j] == 0)
            {
                printf("bug\n");
                exit(0);
            }
        }
    }
    
    free(temp[0]);
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
    int i,j,k,l,m,n,ret_code;
    int **parent_gene,**child_gene,*index[2];
    RealNumber *parent_distance,*parent_entropy,*child_distance,*child_entropy,*parent_score,*child_score;
    double Temperature,t1,t2;
    int parent_size,child_size;
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
    parent_size = world_stars;
    parent_size = (parent_size > 100)? parent_size : 100;
    child_size = 4 * parent_size;
    
    child_size = (parent_size > child_size)? parent_size : child_size;
    parent_gene = (int **)malloc((parent_size + child_size)*sizeof(int *));
    if(parent_gene == NULL)
    {
        ret_code = -1;
        goto error04;
    }
    parent_gene[0] = (int *)malloc(world_stars*(parent_size + child_size)*sizeof(int));
    if(parent_gene[0] == NULL)
    {
        ret_code = -1;
        goto error05;
    }
    for(i=1;i<parent_size + child_size;i++)
    {
        parent_gene[i] = parent_gene[i-1] + world_stars;
    }
    child_gene = parent_gene + parent_size;
    
    i = (parent_size > child_size)? parent_size : child_size;
    index[0] = (int *)malloc(2*i*sizeof(int));
    if(index[0] == NULL)
    {
        ret_code = -1;
        goto error06;
    }
    index[1] = index[0] + i;
    parent_distance = (RealNumber *)malloc(3*(parent_size + child_size)*sizeof(RealNumber));
    if(parent_distance == NULL)
    {
        ret_code = -1;
        goto error07;
    }
    parent_entropy = parent_distance + parent_size;
    parent_score = parent_entropy + parent_size;
    child_distance = parent_score + parent_size;
    child_entropy = child_distance + child_size;
    child_score = child_entropy + child_size;

    entropy_table = (RealNumber *)malloc((parent_size+1)*sizeof(RealNumber));
    if(entropy_table == NULL)
    {
        ret_code = -1;
        goto error08;
    }
    distance_table = (RealNumber **)malloc(world_stars*sizeof(RealNumber *));
    if(distance_table == NULL)
    {
        ret_code = -1;
        goto error09;
    }
    distance_table[0] = (RealNumber *)malloc(world_stars*world_stars*sizeof(RealNumber));
    if(distance_table[0] == NULL)
    {
        ret_code = -1;
        goto error10;
    }
    for(i=1;i<world_stars;i++)
    {
        distance_table[i] = distance_table[i-1] + world_stars;
    }
    probability_table = (int **)malloc(world_stars * sizeof(int *));
    if(probability_table == NULL)
    {
        ret_code = -1;
        goto error11;
    }
    probability_table[0] = (int *)malloc(world_stars*world_stars*sizeof(int));
    if(probability_table[0] == NULL)
    {
        ret_code = -1;
        goto error12;
    }
    for(i=1;i<world_stars;i++)
    {
        probability_table[i] = probability_table[i-1] + world_stars;
    }
    
    initialize_entropy_table(parent_size);
    initialize_distance_table();
    genesis(child_gene,child_size);
    initialize_probability_table();
    for(i=0;i<child_size;i++)
    {
        child_distance[i] = 0.0;
    }
    Temperature = 1.0;
    minimize_free_energy(child_gene,child_distance,child_size,parent_gene,parent_distance,parent_size,Temperature);
    evaluate_distance(parent_gene,parent_distance,parent_size);
    evaluate_entropy(parent_gene,&t1,parent_size);
    t2 = 0.0;
    for(i=0;i<parent_size;i++)
    {
        t2 += (double)parent_distance[i];
    }
    Temperature = t2 / t1;

    i=0;
    R_sort(parent_distance,index[0],index[1],&j,parent_size);
    t2 = parent_distance[index[j][0]];
    while(1)
    {
        R_sort(parent_distance,index[0],index[1],&j,parent_size);
        printf("%8d %f %f\n",i++,Temperature,parent_distance[index[j][0]]);
        child_generator(parent_gene,index[j],parent_size,child_gene,child_size);
        evaluate_distance(child_gene,child_distance,child_size);
        minimize_free_energy(child_gene,child_distance,child_size,parent_gene,parent_distance,parent_size,Temperature);
        Temperature *= 0.99;
    }

    free(probability_table[0]);
    error12:
    free(probability_table);
    error11:
    free(distance_table[0]);
    error10:
    free(distance_table);
    error09:
    free(entropy_table);
    error08:
    free(parent_distance);
    error07:
    free(index[0]);
    error06:
    free(parent_gene[0]);
    error05:
    free(parent_gene);
    error04:
    free(world_points);
    error03:
    free(lp_urandom);
    error02:
    close(urandom_fd);
    error01:
    return ret_code;
}
