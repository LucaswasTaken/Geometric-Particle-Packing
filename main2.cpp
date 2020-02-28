#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include "stdc++.h"

using namespace std;
void cell_divide(int n_particle, float *particle_x, float *particle_y, float *particle_r,  int *particle_cell, int *cell,  float  max_r, int n_cellsv)
{
    int cell_now;
    for (int id = 0; id < n_particle; ++id) {
        cell_now = (int)(particle_x[id]/(2*max_r))*n_cellsv + (int)(particle_y[id]/(2*max_r));
        particle_cell[id]= cell_now;
        cell[id]= cell_now;
    }
}

void cell_count(int n_particle, int* particle_cell,int* cell_start,int* cell_end, int n_cells)
{
    for (int id = 0; id < n_particle; ++id) {
        if((id>0)&&(particle_cell[id-1]<particle_cell[id]))
        {
            cell_start[particle_cell[id]] = id;
        }

        if((particle_cell[id]< n_cells-1 )&&(particle_cell[id+1]>particle_cell[id]))
        {
            cell_end[particle_cell[id]] = id;
        }
        if (id == n_particle-1 )
        {
            cell_end[particle_cell[id]] = id;
        }
    }

}

void interactive_solve(int n_particle, float *particle_x,float *particle_y,float *particle_r, int *particle_cell,int *particle_id, int *cell_start, int *cell_end, float *particle_sobx ,float *particle_soby, float  max_r, int n_cells, int n_cellsv,float *max_sob)
{
    float distx;
    float disty;
    float sobx;
    float soby;
    float tang;
    float sen;
    float cosen;
    int j;
    int c;
    int per_cell;
    int check;
    int cell_now;
    for (int id = 0; id < n_particle; ++id) {
        for (c = 0;c<9;c++)
        {
            check = 0;
            if (c<3)
            {
                cell_now = (int)(particle_cell[id]) - n_cellsv -1+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    if (((cell_now > 0) && (cell_end[cell_now]!=0))||((cell_now == 0)))
                        check = 1;

                }
            }
            else if (c<6)
            {
                cell_now = (int)(particle_cell[id]) -4+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    if (((cell_now > 0) && (cell_end[cell_now]!=0))||((cell_now == 0)))
                        check = 1;
                }
            }
            else
            {
                cell_now = (int)(particle_cell[id])  - 7 + n_cellsv + c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    if (((cell_now > 0) && (cell_end[cell_now]!=0))||((cell_now == 0)))
                        check = 1;
                }
            }

            if (check==1)
            {
                for (per_cell = cell_start[cell_now]; per_cell < (cell_end[cell_now]+1); ++per_cell)
                {
                    if(particle_id[per_cell]!=id)
                    {
                        j = particle_id[per_cell];
                        distx = fabs(particle_x[id]-particle_x[j]);
                        disty = fabs(particle_y[id]-particle_y[j]);
                        if (distx>0.05*disty)
                        {
                            tang = disty/distx;
                            sen = (tang/(tang*tang+1));
                            cosen = (1/(tang*tang+1));
                        }

                        else
                        {
                            sen = 1.0;
                            cosen = 0.0;
                        }
                        if(sqrt(distx*distx+disty*disty)<(particle_r[id]+particle_r[j]))
                        {

                            sobx = ((particle_r[id]+particle_r[j]) - (sqrt(distx*distx+disty*disty)))*cosen;
                            soby = ((particle_r[id]+particle_r[j]) - (sqrt(distx*distx+disty*disty)))*sen;
                            if ((particle_x[id]<particle_x[j]))
                            {
                                particle_sobx[id] = particle_sobx[id]-sobx/2;
                            }
                            else if ((particle_x[id]>particle_x[j]))
                            {
                                particle_sobx[id] = particle_sobx[id]+sobx/2;

                            }


                            if ((particle_y[id]<particle_y[j]))
                            {
                                particle_soby[id] = particle_soby[id]-soby/2;
                            }
                            else if((particle_y[id]>particle_y[j]))
                            {
                                particle_soby[id] = particle_soby[id]+soby/2;
                            }

                            else
                            {
                                if ((id < j))
                                {
                                    particle_soby[id] = particle_soby[id]-soby/2;
                                }
                                else
                                {
                                    particle_soby[id] = particle_soby[id]+soby/2;

                                }

                            }
                            if (sqrt(particle_sobx[id]*particle_sobx[id]+particle_soby[id]*particle_soby[id])>0.01*particle_r[id])
                            {
                                if(sqrt(particle_sobx[id]*particle_sobx[id]+particle_soby[id]*particle_soby[id])>max_sob[0])
                                {
                                    max_sob[0] = sqrt(particle_sobx[id]*particle_sobx[id]+particle_soby[id]*particle_soby[id]);
                                }
                            }
                        }

                    }
                }
            }
        }
    }





}

void particle_erase(int n_particle, float* particle_x, float* particle_y, float* particle_sobx, float* particle_soby, int* particle_id , float* particle_r, float altura, float largura)
{
    for (int id = 0; id <n_particle ; ++id) {
        particle_x[id] = particle_x[id] + particle_sobx[id];
        particle_y[id] = particle_y[id] + particle_soby[id];
        particle_sobx[id] = 0;
        particle_soby[id] = 0;
        particle_id[id] = id;
        if ((particle_x[id]>largura-particle_r[id]))
        {
            particle_x[id] = largura-particle_r[id];
        }
        if ((particle_x[id]<0+particle_r[id]))
        {
            particle_x[id] = 0+particle_r[id];
        }
        if ((particle_y[id]<0+particle_r[id]))
        {
            particle_y[id] = 0+particle_r[id];
        }
        if ((particle_y[id]>altura-particle_r[id]))
        {
            particle_y[id] = altura-particle_r[id];
        }
    }
    
}

void cell_erase(int n_particle, int* cell_start, int* cell_end)
{
    for (int id = 0; id < n_particle; ++id) {
        cell_start[id] = 0;
        cell_end[id] = 0;
    }

}

void sob_test(float *max_sob, float min_r)
{
    if(max_sob[0]<0.05*min_r)
    {
        max_sob[1] = 1;
    }
    else
    {
        max_sob[0] = 0;
    }
}

void final_delete(int n_particle, float *particle_x,float *particle_y,float *particle_r, int *particle_cell,int *particle_id, int *cell_start, int *cell_end, float *particle_sobx ,float *particle_soby, float  max_r, int n_cells, int n_cellsv,float *max_sob)
{
    float distx;
    float disty;
    float sobx;
    float soby;
    float tang;
    float sen;
    float cosen;
    int j;
    int c;
    int per_cell;
    int check;
    int cell_now;
    for (int id = 0; id < n_particle ; ++id)
    {
        for (c = 0;c<9;c++)
        {
            check = 0;
            if (c<3)
            {
                cell_now = (int)(particle_cell[id]) - n_cellsv -1+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    if (((cell_now > 0) && (cell_end[cell_now]!=0))||((cell_now == 0)))
                        check = 1;

                }
            }
            else if (c<6)
            {
                cell_now = (int)(particle_cell[id]) -4+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    if (((cell_now > 0) && (cell_end[cell_now]!=0))||((cell_now == 0)))
                        check = 1;
                }
            }
            else
            {
                cell_now = (int)(particle_cell[id])  - 7 + n_cellsv + c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    if (((cell_now > 0) && (cell_end[cell_now]!=0))||((cell_now == 0)))
                        check = 1;
                }
            }

            if (check==1)
            {
                for (per_cell = cell_start[cell_now]; per_cell < (cell_end[cell_now]+1); ++per_cell)
                {
                    if(particle_id[per_cell]!=id)
                    {
                        j = particle_id[per_cell];
                        distx = fabs(particle_x[id]-particle_x[j]);
                        disty = fabs(particle_y[id]-particle_y[j]);
                        if (distx>0.05*disty)
                        {
                            tang = disty/distx;
                            sen = (tang/(tang*tang+1));
                            cosen = (1/(tang*tang+1));
                        }

                        else
                        {
                            sen = 1.0;
                            cosen = 0.0;
                        }
                        if(sqrt(distx*distx+disty*disty)<(particle_r[id]+particle_r[j]))
                        {
                            sobx = ((particle_r[id]+particle_r[j]) - (sqrt(distx*distx+disty*disty)))*cosen;
                            soby = ((particle_r[id]+particle_r[j]) - (sqrt(distx*distx+disty*disty)))*sen;
                            if ((particle_x[id]<particle_x[j]))
                            {
                                particle_sobx[id] = particle_sobx[id]-sobx/2;
                            }
                            else if ((particle_x[id]>particle_x[j]))
                            {
                                particle_sobx[id] = particle_sobx[id]+sobx/2;

                            }


                            if ((particle_y[id]<particle_y[j]))
                            {
                                particle_soby[id] = particle_soby[id]-soby/2;
                            }
                            else if((particle_y[id]>particle_y[j]))
                            {
                                particle_soby[id] = particle_soby[id]+soby/2;
                            }

                            else
                            {
                                if ((id < j))
                                {
                                    particle_soby[id] = particle_soby[id]-soby/2;
                                }
                                else
                                {
                                    particle_soby[id] = particle_soby[id]+soby/2;

                                }

                            }
                            if (sqrt(particle_sobx[id]*particle_sobx[id]+particle_soby[id]*particle_soby[id])>0.05*particle_r[id])
                            {
                                particle_x[id] = 0;
                                particle_y[id] = 0;
                                particle_r[id] = 0;
                            }
                        }

                    }
                }
            }
        }
    }




}

int main(int argc, char** argv) {


    float altura ;
    float largura;
    float porosidade;
    float volume_ocupado;
    int N_granulometria;

    float max_r;
    float min_r;
    float mean_r;
    int n_particle_t = 0;
    float **granulometrias;
    int *n_particle;

    //Dados iniciais
    if(argc<3){
        altura = 250;
        largura = 250;
        porosidade = 0.2;
        volume_ocupado = (1 - porosidade)*(altura*largura);
        N_granulometria = 1;
        granulometrias =  (float**)malloc(N_granulometria * sizeof(float*));
        for (int index=0;index<N_granulometria;++index)
        {
            granulometrias[index] = (float*)malloc(2 * sizeof(float));
        }
        max_r = 0;
        min_r = 100;
        mean_r = 0;

        for (int i = 0; i < 2 * N_granulometria; i += 2) {
            granulometrias[i / 2][0] = 1;
            //cout<<granulometrias[i/2][0]<<endl;
            granulometrias[i / 2][1] = 1;
            //cout<<granulometrias[i/2][1]<<endl;
            if (granulometrias[i / 2][0] > max_r)
                max_r = granulometrias[i / 2][0];
            if (granulometrias[i / 2][0] < min_r)
                min_r = granulometrias[i / 2][0];
            mean_r = mean_r + granulometrias[i / 2][0] * granulometrias[i / 2][1];
        }
        //Geracao aleatoria das particulas

        n_particle = (int *)malloc(N_granulometria * sizeof(int));

        for (int i = 0; i < N_granulometria; ++i) {
            n_particle[i] = int(floor(volume_ocupado*granulometrias[i][1] / (M_PI*granulometrias[i][0] * granulometrias[i][0])));
            n_particle_t += n_particle[i];
        }
    }
    else
        {
            N_granulometria = 1;
            granulometrias =  (float**)malloc(N_granulometria * sizeof(float*));
            for (int index=0;index<N_granulometria;++index)
            {
                granulometrias[index] = (float*)malloc(2 * sizeof(float));
            }
            porosidade = 0.2;
            n_particle_t = stoi(argv[2]);
            altura = sqrt(M_PI*n_particle_t/(1-porosidade));
            largura = altura;
            volume_ocupado = (1 - porosidade)*(altura*largura);
            max_r = 0;
            min_r = 100;
            mean_r = 0;
            for (int i = 0; i < 2 * N_granulometria; i += 2) {
                granulometrias[i / 2][0] = 1;
                //cout<<granulometrias[i/2][0]<<endl;
                granulometrias[i / 2][1] = 1;
                //cout<<granulometrias[i/2][1]<<endl;
                if (granulometrias[i / 2][0] > max_r)
                    max_r = granulometrias[i / 2][0];
                if (granulometrias[i / 2][0] < min_r)
                    min_r = granulometrias[i / 2][0];
                mean_r = mean_r + granulometrias[i / 2][0] * granulometrias[i / 2][1];
            }
            n_particle = (int *)malloc(N_granulometria * sizeof(int));
            n_particle[0] = n_particle_t;
    }
    int * particle_cell = new  int[n_particle_t];
    int * cell = new  int[n_particle_t];
    int * particle_id = new int[n_particle_t];
    for (int i = 0; i < n_particle_t; i++)
    {
        //m_hInput[i] = m_N - i;			// Use this for debugging. Use 1 or i or similar
        particle_cell[i] = 0;
        cell[i] = 0;
        particle_id[i] = i;
    }
    int		*cell_start;
    int		*cell_end;
    float				*particle_x;
    float				*particle_y;
    float				*particle_sobx;
    float				*particle_soby;
    float				*particle_r;
    float               max_sob[2];
    int		n_cellsv;
    int		n_cells;

    n_cellsv = (int(altura / (2 * max_r))) + 1;
    n_cells = n_cellsv * (int(largura / (2 * max_r)) + 1);


    /*Criacao dos data_arrays de particulas:*/

    particle_x = (float *)malloc((n_particle_t) * sizeof(float));
    particle_y = (float *)malloc((n_particle_t) * sizeof(float));
    particle_sobx = (float *)malloc((n_particle_t) * sizeof(float));
    particle_soby = (float *)malloc((n_particle_t) * sizeof(float));
    particle_r = (float *)malloc((n_particle_t) * sizeof(float));


    //Inicializa��o das part�culas
    //cout<<n_particle[0]<<endl;
    int aux_particle_alloc = 0;
    for (int i = 0; i < N_granulometria; ++i) {
        for (int j = 0; j < n_particle[i]; ++j) {
            particle_x[aux_particle_alloc] = (rand() % (1000))*largura / 1000;
            //cout<<particle_x[aux_particle_alloc]<<endl;
            particle_sobx[aux_particle_alloc] = 0;
            particle_y[aux_particle_alloc] = (rand() % (1000))*altura / 1000;
            particle_soby[aux_particle_alloc] = 0;
            particle_r[aux_particle_alloc] = granulometrias[i][0];
            aux_particle_alloc += 1;
        }
    }

    //ini cell_number
    cell_start = new int[n_cells];
    cell_end = new int[n_cells];
    max_sob[0] = 0;
    max_sob[1] = 0;

    //run selected task

    int n_int = 0;
    //cl_int ret;
    clock_t timer = clock();
    clock_t timer_sort = 0;
    clock_t timer_sortaux;
    vector< pair <int,int> > cell_ind_vet;



    while (n_int < 2000)
    {
        particle_erase(n_particle_t, particle_x, particle_y, particle_sobx, particle_soby, particle_id, particle_r, altura, largura);
        cell_erase(n_particle_t, cell_start,cell_end);
        cell_divide(n_particle_t, particle_x, particle_y, particle_r,  particle_cell, cell,  max_r, n_cellsv);
        timer_sortaux = clock();
        for (int i = 0; i <n_particle_t; ++i)
        {
            particle_id[i] = i;
            cell_ind_vet.emplace_back(make_pair(particle_cell[i],particle_id[i]));
        }
        sort(cell_ind_vet.begin(),cell_ind_vet.end());
        for (int i = 0; i <n_particle_t; ++i)
        {
            particle_id[i] = cell_ind_vet[i].second;
            particle_cell[i] = cell_ind_vet[i].first;
        }
        timer_sortaux = clock() - timer_sortaux;
        timer_sort = timer_sort + timer_sortaux;
        cell_count(n_particle_t, particle_cell,cell_start,cell_end, n_cells);
        interactive_solve(n_particle_t, particle_x, particle_y, particle_r, particle_cell,particle_id, cell_start, cell_end, particle_sobx, particle_soby,  max_r,  n_cells,  n_cellsv, max_sob);
        sob_test(max_sob, min_r);
        if(int(max_sob[1])==1)
            break;
        n_int += 1;



    }
    particle_erase(n_particle_t, particle_x, particle_y, particle_sobx, particle_soby, particle_id, particle_r, altura, largura);
    cell_erase(n_particle_t, cell_start,cell_end);
    cell_divide(n_particle_t, particle_x, particle_y, particle_r,  particle_cell, cell,  max_r, n_cellsv);
    timer_sortaux = clock();
    for (int i = 0; i <n_particle_t; ++i)
    {
        particle_id[i] = i;
        cell_ind_vet.emplace_back(make_pair(particle_cell[i],particle_id[i]));
    }
    sort(cell_ind_vet.begin(),cell_ind_vet.end());
    for (int i = 0; i <n_particle_t; ++i)
    {
        particle_id[i] = cell_ind_vet[i].second;
        particle_cell[i] = cell_ind_vet[i].first;
    }
    timer_sortaux = clock() - timer_sortaux;
    timer_sort = timer_sort + timer_sortaux;
    cell_count(n_particle_t, particle_cell,cell_start,cell_end, n_cells);
    final_delete(n_particle_t, particle_x, particle_y, particle_r, particle_cell,particle_id, cell_start, cell_end, particle_sobx, particle_soby,  max_r,  n_cells,  n_cellsv, max_sob);;
    sob_test(max_sob, min_r);
    float timer_sortfloat = float(timer_sort);
    float timer_solvefloat;
    timer = clock() - timer;
    float timer_float = float(timer);
    timer_float = float(timer_float/CLOCKS_PER_SEC);
    timer_sortfloat = float(timer_sortfloat/CLOCKS_PER_SEC);
    timer_solvefloat = timer_float - timer_sortfloat;
    cout << "Time: " << timer_float << endl;
    //read back the results synchronously.
}
