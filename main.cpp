//Includes
//http://arkanis.de/weblog/2014-11-25-minimal-opencl-development-on-windows
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include "stdc++.h"

using namespace std;

//Determinação do número de partículas máximo por célula
int cell_max(float *particles_x, float *particles_y, int n_particles, int n_cells, float max_r, int n_cellsv)
{

    int cells[n_cells];
    int cell_now;


    for (int i = 0; i < n_cells; ++i)
    {
        cells[i] = 0;
    }

    for (int i = 0; i < n_particles; ++i)
    {
        cell_now = (int)(particles_y[i]/(2*max_r)) + (int)(particles_x[i]/(2*max_r))*n_cellsv;
        cells[cell_now]+=1;
    }

    int max_cell = 0;
    for (int i = 0; i < n_cells; ++i)
    {
        if (max_cell<cells[i])
        {
            max_cell = cells[i];
        }
    }

    return max_cell;
}


void particle_divide(float *particles_x,float *particles_y,float *particles_r, int *particles_cell, float *particles_sobx ,float *particles_soby, int n_cellsv, int n_particles_t, float  max_r, int max_cell, float altura,float largura)
{
    for (int id = 0; id < n_particles_t; ++id)
    {

        int cell_now;
        particles_x[id] = particles_x[id] + particles_sobx[id];
        particles_y[id] = particles_y[id] + particles_soby[id];
        particles_sobx[id] = 0;
        particles_soby[id] = 0;

        if ((particles_x[id]>largura-particles_r[id]))
        {
            particles_x[id] = largura-particles_r[id];
        }
        if ((particles_x[id]<0+particles_r[id]))
        {
            particles_x[id] = 0+particles_r[id];
        }
        if ((particles_y[id]<0+particles_r[id]))
        {
            particles_y[id] = 0+particles_r[id];
        }
        if ((particles_y[id]>altura-particles_r[id]))
        {
            particles_y[id] = altura-particles_r[id];
        }

        cell_now = (int)(particles_x[id]/(2*max_r))*n_cellsv + (int)(particles_y[id]/(2*max_r));
        particles_cell[id]= (float)(cell_now);

    }
}



void cell_divide(int* particles_cell, int* cells, int max_cell, int n_particles_t )
{
    int cell_now;
    int n_particles_in_cell;
    int cell_id;
    for (int i = 0; i <n_particles_t; ++i)
    {
        cell_now  = (int)(particles_cell[i]);
        cell_id = max_cell*cell_now;
        cells[cell_id]+=1;
        n_particles_in_cell = cells[cell_id];
        cells[cell_id + n_particles_in_cell] = i;
    }
}



void cell_erase(int* cells, int max_cell, int n_cells, float *cell_PArea, int * id_smaller_particle, int *aux_id_cell)
{

    for (int i = 0; i <n_cells*max_cell; ++i)
    {
        cells[i]=0;
    }

    for (int i = 0; i <n_cells; ++i)
    {
        cell_PArea[i]=0;
        aux_id_cell[i] = i;
        id_smaller_particle[i] = 0;
    }

}



void interactive_solve(float *particles_x,float *particles_y,float *particles_r, int *particles_cell, float *particles_sobx ,float *particles_soby,int *cells,float *max_overlap, int n_particles_t, int n_cellsv, int n_cells, int max_cell, int *max_sob_id , float *max_sob_value)
{
    for (int id = 0; id < n_particles_t; ++id)
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

        for (c = 0;c<9;c++)
        {
            check = 0;
            if (c<3)
            {
                cell_now = (int)(particles_cell[id]) - n_cellsv -1+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    check = 1;
                }
            }
            else if (c<6)
            {
                cell_now = (int)(particles_cell[id]) -4+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    check = 1;
                }
            }
            else
            {
                cell_now = (int)(particles_cell[id])  - 7 + n_cellsv + c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    check = 1;
                }
            }

            if (check==1)
            {
                for (per_cell = 1; per_cell < cells[cell_now*max_cell]+1; ++per_cell) {
                    if(cells[cell_now*max_cell+per_cell]!=id)
                    {
                        j = cells[cell_now*max_cell+per_cell];
                        distx = fabs(particles_x[id]-particles_x[j]);
                        disty = fabs(particles_y[id]-particles_y[j]);
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
                        if(sqrt(distx*distx+disty*disty)<(particles_r[id]+particles_r[j]+0.005*(particles_r[id]+particles_r[j])))
                        {

                            sobx = ((particles_r[id]+particles_r[j]) - (sqrt(distx*distx+disty*disty)))*cosen;
                            soby = ((particles_r[id]+particles_r[j]) - (sqrt(distx*distx+disty*disty)))*sen;
                            if ((particles_x[id]<particles_x[j]))
                            {
                                particles_sobx[id] = particles_sobx[id]-sobx/2;
                            }
                            else if ((particles_x[id]>particles_x[j]))
                            {
                                particles_sobx[id] = particles_sobx[id]+sobx/2;

                            }
                            else
                            {
                                if ((id < j))
                                {
                                    particles_sobx[id] = particles_sobx[id]-sobx/2;
                                }
                                else
                                {
                                    particles_sobx[id] = particles_sobx[id]+sobx/2;

                                }

                            }

                            if ((particles_y[id]<particles_y[j]))
                            {
                                particles_soby[id] = particles_soby[id]-soby/2;
                            }
                            else if ((particles_y[id]>particles_y[j]))
                            {
                                particles_soby[id] = particles_soby[id]+soby/2;
                            }
                            else
                            {
                                if ((id < j))
                                {
                                    particles_soby[id] = particles_soby[id]-soby/2;
                                }
                                else
                                {
                                    particles_soby[id] = particles_soby[id]+soby/2;
                                }
                            }
                            if((particles_r[id]+particles_r[j]-sqrt(distx*distx+disty*disty))>max_sob_value[id])
                            {
                                max_sob_value[id] = (particles_r[id]+particles_r[j]-sqrt(distx*distx+disty*disty))/particles_r[id];
                                max_sob_id[id] = id;
                            }


                            if((particles_r[id]+particles_r[j]-sqrt(distx*distx+disty*disty))>max_overlap[1])
                            {
                                max_overlap[1] = (particles_r[id]+particles_r[j]-sqrt(distx*distx+disty*disty))/particles_r[id];
                                if (max_overlap[1]>0.05)
                                {
                                    max_overlap[0]=1;
                                }
                            }
                        }

                    }
                }
            }
        }
    }
}

//

void area_cell_compute(float *particles_r, int *particles_cell, int n_cellsv,int n_particles_t, int n_cells, float *cell_PArea, int * id_smaller_particle, float max_r )

{
    int cell_now;
    float area_now;
    float * area_before = (float *)malloc(n_cells * sizeof(float));
    for (int k = 0; k < n_cells; ++k) {
        area_before[k] = (float)M_PI*max_r*max_r;
    }
    for (int i = 0; i < n_particles_t; ++i) {
        cell_now  = (int)(particles_cell[i]);
        area_now = (float)((particles_r[i])*(particles_r[i])*M_PI);
        cell_PArea[cell_now] = cell_PArea[cell_now] + area_now;
        if(area_now<area_before[cell_now])
        {
            id_smaller_particle[cell_now] = i;
            area_before[cell_now] = area_now;
        }
    }

    int c, check, count_area;
    float mean_area;
    for (int j = 0; j < n_cells ; ++j)

    {
        mean_area = 0;
        count_area = 0;
        for (c = 0;c<9;c++)
        {
            check = 0;
            if (c<3)
            {
                cell_now = j - n_cellsv -1+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    check = 1;
                }
            }
            else if (c<6)
            {
                cell_now = j -4+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    check = 1;
                }
            }
            else
            {
                cell_now = j  - 7 + n_cellsv + c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    check = 1;
                }
            }

            if (check==1)
            {
                mean_area = cell_PArea[cell_now] + mean_area;
                count_area = count_area +1;
            }
        }
        if (count_area>0)
        {
            cell_PArea[j] = cell_PArea[j] + 2*mean_area/count_area;
        }
    }

}


void realloc_area(float *particles_x,float *particles_y,float *particles_r, int *particles_cell,  int n_cells, int n_cellsv, float max_r, vector< pair <float,int> > cell_PArea_vec, int * id_smaller_particle, float mean_r)
{
    int particle_now,cell_now, cell_after;
    int aux_count = 0;
    float r_now;
    for (int i = 0; i < int(0.01*n_cells+1); ++i) {

        cell_after = cell_PArea_vec[aux_count].second;
        cell_now = cell_PArea_vec[n_cells-aux_count-1].second;
        particle_now = id_smaller_particle[cell_now];
        r_now = particles_r[particle_now];

        if (r_now <= mean_r)
        {
            particles_cell[particle_now] =cell_after;
            particles_x[particle_now] = int(floor(cell_after/n_cellsv)*2*max_r) + (rand()%int(2*max_r));
            particles_y[particle_now] = int(floor(cell_after%n_cellsv)*2*max_r) + (rand()%int(2*max_r));
            cout<<"teteeeeee\n\n\n\n"<<r_now<<endl;

        }
        else
        {
            i=i-1;
            aux_count = aux_count +1;
            continue;
        }
    }
}



void p_remove(float *particles_x,float *particles_y,float *particles_r, int *particles_cell, int n_particles_t, int n_cells, float max_r, vector< pair <float,int> > sob_vet, int* insert_index, float* insert_r, float mean_r)
{
    int particle_now;
    int aux_count = 0;
    float r_now;
    for (int i = 0; i < int(0.01*n_cells+1); ++i) {

        particle_now = sob_vet[n_particles_t- aux_count -1].second;
        r_now = particles_r[particle_now];


        if (r_now <= mean_r)
        {
            insert_index[0] = insert_index[0] +1;
            insert_index[insert_index[0]] = particle_now;
            insert_r[0] = insert_r[0]+1;
            insert_r[int(insert_r[0])] = int(particles_r[particle_now]);
            particles_cell[particle_now] = 0;
            particles_x[particle_now] = 0;
            particles_y[particle_now] = 0;
            particles_r[particle_now] = 0;
            cout<<"teteeeeee\n\n\n\n"<<r_now<<endl;

        }
        else
        {
            i=i-1;
            aux_count = aux_count +1;
            continue;
        }



    }
}


void p_insert(float *particles_x,float *particles_y,float *particles_r, int *particles_cell, int n_cellsv, float max_r, vector< pair <float,int> > cell_PArea_vec, int* insert_index, float* insert_r)
{
    int particle_now;
    int  cell_after;
    int aux_count = 0;
    float r_now;

    for (int i = 1; i < (insert_index[0]+1); ++i) {

        cell_after = cell_PArea_vec[aux_count].second;
        particle_now = insert_index[i];
        r_now = insert_r[i];

        particles_cell[particle_now] =cell_after;
        particles_x[particle_now] = int(floor(cell_after/n_cellsv)*2*max_r) + (rand()%int(2*max_r));
        particles_y[particle_now] = int(floor(cell_after%n_cellsv)*2*max_r) + (rand()%int(2*max_r));
        particles_r[particle_now] = r_now;


    }

}

void grain_print(string str_grain, int N_granulometria, float granulometrias[][2], int *n_particles, vector< pair <float,int> > cell_PArea_vec, int insert_aux, float* raios_r)
{
    int * granulometria_after = (int *)malloc(N_granulometria * sizeof(int));
    for (int j = 0; j <N_granulometria ; ++j) {
        granulometria_after[j] = n_particles[j] - raios_r[j];
    }

    ofstream grain_poro;
    grain_poro.open(str_grain,std::ofstream::out | std::ofstream::app);
    grain_poro << "granulometria inicial ----> granulometria final\n";
    for (int k = 0; k <N_granulometria ; ++k) {
        grain_poro << granulometrias[k][0] <<"---->"<<n_particles[k]<<"---->"<<raios_r[k] <<endl;
    }

    grain_poro.close();
}

int main(int argc, char** argv)
{

    clock_t timer = clock();
    double timetotal;
    clock_t timer2;
    if ( argc < 6 )
    {
        printf("uso: ./particlegeneration  [Altura do Domínio] [Comprimento do Domínio] [Porosidade] "
                       "[Numero de Granumetrias] [raio 1] [Granulometria 1] [raio 2] [Granulometria 2] ....\n");
        return -1;
    }
//    Declaracao das Variaveis de Auxilio

    float altura = stof(argv[2],NULL);
    float largura = stof(argv[1],NULL);
    float porosidade = stof(argv[3],NULL);
    float volume_ocupado = (1-porosidade)*(altura*largura);
    int N_granulometria = stoi(argv[4],NULL);
    float granulometrias[N_granulometria][2];
    float raios_r[N_granulometria];
//    float altura = 5;
//    float largura = 5;
//    float porosidade = 0.2;
//    float volume_ocupado = (1-porosidade)*(altura*largura);
//    int N_granulometria = 1;
//    float granulometrias[N_granulometria][2];

    float max_r = 0;
    float min_r = 100;
    float mean_r = 0;

    for (int i = 0; i < 2*N_granulometria ; i+=2) {
        granulometrias[i/2][0] = stof(argv[5+i],NULL);
        //cout<<granulometrias[i/2][0]<<endl;
        granulometrias[i/2][1] = stof(argv[6+i],NULL);
        //cout<<granulometrias[i/2][1]<<endl;
        if (granulometrias[i/2][0]>max_r)
            max_r = granulometrias[i/2][0];
        if (granulometrias[i/2][0]<min_r)
            min_r = granulometrias[i/2][0];
        mean_r = mean_r + granulometrias[i/2][0]*granulometrias[i/2][1];
    }


    //Geracao aleatoria das particulas

    int *n_particles = (int *)malloc(N_granulometria * sizeof(int));
    int n_particles_t = 0;
    for (int i = 0; i < N_granulometria ; ++i) {
        n_particles[i] = int(floor(volume_ocupado*granulometrias[i][1]/(M_PI*granulometrias[i][0]*granulometrias[i][0])));
        n_particles_t+=n_particles[i];
        raios_r[i] = 0;
    }


    /*Criacao dos data_arrays de particulas:*/

    int count_granu =0;
    float *particles_x = (float *)malloc((n_particles_t)* sizeof(float));
    float *particles_y = (float *)malloc((n_particles_t)* sizeof(float));
    float *particles_r = (float *)malloc((n_particles_t)* sizeof(float));
    float *particles_sobx = (float *)malloc((n_particles_t)* sizeof(float));
    float *particles_soby = (float *)malloc((n_particles_t)* sizeof(float));
    int *particles_cell = (int *)malloc((n_particles_t)* sizeof(int));





    //Inicialização das partículas


    int aux_particle_alloc = 0;
    for (int i = 0; i < N_granulometria ; ++i) {
        for (int j = 0; j < n_particles[i]; ++j) {
            particles_x[aux_particle_alloc] =  (rand()%(1000))*largura/1000;
            particles_y[aux_particle_alloc] = (rand()%(1000))*altura/1000;
            particles_r[aux_particle_alloc] = granulometrias[i][0];
            particles_sobx[aux_particle_alloc] = 0;
            particles_soby[aux_particle_alloc] = 0;
            aux_particle_alloc+=1;
        }
    }
/*
    //Print OFS inicial
    ofstream myfile;
    string str;
    str.append("generationi");
    str.append(argv[2]);
    str.append("x");
    str.append(argv[1]);
    str.append("x");
    string n_particles_ts = to_string(n_particles_t);
    str.append(n_particles_ts);
    str.append(".ofs");
    myfile.open (str);
    myfile << "%DEM.MATERIAL.COLOR\n1\n1 0.725 0.478 0.341 1.000\n%DEM.PARTICLE\n"<<n_particles_t<<endl<<"%DEM.PARTICLE.CIRCLE\n"<<n_particles_t<<endl;



    for (int i = 0; i < n_particles_t; ++i) {
        myfile <<i+1<<" 1"<<" "<<particles_r[i]<<" "<<particles_x[i]<<" "<<particles_y[i]<<" 0.00"<<endl;
    }
    myfile.close();
    */

    //Auxiliares de célula

    int  n_int = 0;
    int n_cellsv = (int(altura/(2*max_r))) + 1;
    int n_cells = n_cellsv*(int(largura/(2*max_r))+1);



    //Estimando o numero de particulas maximo por celula
    int max_cell = int(1.4*cell_max(particles_x, particles_y, n_particles_t, n_cells, max_r, n_cellsv));
    //Criando as cells na host
    int *cells_host = (int *)malloc(max_cell*n_cells * sizeof(int));

    //Criando estruturas auxiliares de indice
    int *aux_id_cell = (int *)malloc(n_cells * sizeof(int));
    float *cell_PArea = (float *)malloc(n_cells * sizeof(float));
    int * id_smaller_particle = (int *)malloc(n_cells * sizeof(int));
    int *insert_index = (int *)malloc(int(n_particles_t+1) * sizeof(int));
    float *insert_r = (float *)malloc(int(n_particles_t+1) * sizeof(float));
    int insert_trials = 0;
    insert_index[0] = 0;
    insert_r[0] = 0;
    int index_aux = 0;

    //Criação das variáveis auxiliares de sobreposição

    int *max_sob_id = (int *)malloc((n_particles_t)* sizeof(int));
    float *max_sob_value = (float *)malloc((n_particles_t)* sizeof(float));
    float *max_overlap = (float *)malloc((2)* sizeof(float));
    max_overlap[0] = 1.0;
    max_overlap[1] = 10*max_r;

    //inicialização das células
    for (int i = 0; i < n_cells*max_cell; ++i)
    {
        cells_host[i] = 0;
    }

    vector< pair <float,int> > cell_PArea_vet;
    vector< pair <float,int> > sob_vet;


    //auxiliares de realocação e remoção

    for (int i = 0; i <n_cells; ++i)
    {
        cell_PArea[i] = 0.0;
        aux_id_cell[i] = i;
        cell_PArea_vet.emplace_back(make_pair(cell_PArea[i],aux_id_cell[i]));
    }

    for (int i = 0; i <n_particles_t; ++i)
    {
        max_sob_id[i] = i;
        max_sob_value[i] = 0;
        sob_vet.emplace_back(make_pair(max_sob_value[i],max_sob_id[i]));
    }


    //Loop iterativo

    clock_t timer_aux, timer_celld= clock() - clock()
    , timer_partd= clock() - clock(), timer_areac= clock() - clock()
    ,timer_sort= clock() - clock(),timer_insert= clock() - clock(),timer_realloc= clock() - clock()
    , timer_solve= clock() - clock(), timer_remove= clock() - clock();
    while (((max_overlap[1]>0.05)||(int(max_overlap[0])==1)||((insert_index[0]>0)&&(insert_trials<10))))
    // (((max_overlap[1]>0.05)||(int(max_overlap[0])==1)))
    {

        //inicialização das variáveis de sobreposição
        max_overlap[0] = 0;
        max_overlap[1] = 0;

        //divisão em células
        timer_aux = clock();
        particle_divide(particles_x,particles_y,particles_r, particles_cell, particles_sobx ,particles_soby, n_cellsv, n_particles_t, max_r, max_cell, altura,largura);
        timer_partd = timer_partd + timer_aux - clock();
        //Métodos de realocação, remoção e inserção

        if(((n_int>10)&&(n_int%100==0))) {

            cell_PArea_vet.erase(cell_PArea_vet.begin(),cell_PArea_vet.end());
            timer_aux = clock();
            area_cell_compute(particles_r, particles_cell,  n_cellsv, n_particles_t,  n_cells, cell_PArea,  id_smaller_particle, max_r);
            timer_areac = timer_areac+ timer_aux - clock();
            for (int i = 0; i <n_cells; ++i)
            {
                cell_PArea_vet.emplace_back(make_pair(cell_PArea[i],aux_id_cell[i]));

            }

            for (int i = 0; i <n_particles_t; ++i)
            {
                sob_vet.emplace_back(make_pair(max_sob_value[i],max_sob_id[i]));

            }
            timer_aux = clock();
            sort(cell_PArea_vet.begin(),cell_PArea_vet.end());
            sort(sob_vet.begin(),sob_vet.end());
            timer_sort = timer_sort + timer_aux - clock();
            if (((n_int%1000==0)&&(insert_index[0]>0)&&(insert_trials<10))||(max_overlap[1]<0.05))
            {
                timer_aux = clock();
                p_insert(particles_x, particles_y, particles_r, particles_cell,  n_cellsv,  max_r, cell_PArea_vet,  insert_index,  insert_r);
                insert_trials+=1;
                insert_index[0] = 0;
                insert_r[0] = 0;
                timer_insert = timer_insert + timer_aux - clock();
            }

            else if(n_int%500==0)
            {
                timer_aux = clock();
                if (insert_trials>=10) {
                    index_aux = index_aux + insert_index[0];
                    for (int i = 1; i < (insert_index[0]+1); ++i) {
                        for (int j = 0; j <N_granulometria ; ++j) {
                            if (insert_r[i]==granulometrias[j][0])
                            {
                                raios_r[j] = raios_r[j] + 1;
                            }
                        }
                    }

                }
                insert_index[0] = 0;
                insert_r[0] = 0;
                p_remove(particles_x,particles_y,particles_r, particles_cell, n_particles_t,  n_cells,  max_r, sob_vet,  insert_index,  insert_r, mean_r);
                timer_remove= timer_remove + timer_aux - clock();
            }

            else{
                timer_aux = clock();
                realloc_area(particles_x,particles_y,particles_r, particles_cell,  n_cells,  n_cellsv, max_r, cell_PArea_vet,  id_smaller_particle, mean_r);
                timer_realloc = timer_realloc + timer_aux - clock();
            }

        }


        //reinicialização das células
        timer_aux = clock();
        cell_erase(cells_host, max_cell, n_cells,  cell_PArea, id_smaller_particle, aux_id_cell);
        cell_divide(particles_cell, cells_host, max_cell, n_particles_t);
        timer_celld = timer_celld + timer_aux - clock();
        //calculo das sobreposições
        timer_aux = clock();
        interactive_solve(particles_x, particles_y, particles_r, particles_cell, particles_sobx, particles_soby, cells_host,max_overlap, n_particles_t, n_cellsv, n_cells, max_cell, max_sob_id ,max_sob_value);
        n_int=n_int+1;
        timer_solve = timer_solve + timer_aux - clock();
        // print the results
        printf("nint e max_overlap: %d %f\n",n_int,max_overlap[1]);
        if(n_int>2000)
            break;
    }



    timer2 = clock();
    timetotal = ((float)timer2-(float)timer)/CLOCKS_PER_SEC;
    //PrintOFS2
    /*ofstream myfile2;
    string str2, str_grain;
    str2.append("generationf");
    str2.append(argv[2]);
    str2.append("x");
    str2.append(argv[1]);
    str2.append("_(");
    string nints = to_string(n_int);
    string timeS = to_string(timetotal);
    string sobfs = to_string(max_overlap[1]);
    str2.append(n_particles_ts);
    str2.append("p)");
    str2.append("_(");
    str2.append(nints);
    str2.append("int)");
    str2.append("_(");
    str2.append(timeS);
    str2.append("s)");
    str2.append("_(");
    str2.append(sobfs);
    str2.append(")");
    str_grain = str2;
    str2.append(".ofs");
    str_grain.append("_granulometry.grain");
    myfile2.open (str2);

    myfile2 << "%DEM.MATERIAL.COLOR\n1\n1 0.725 0.478 0.341 1.000\n%DEM.PARTICLE\n"<<n_particles_t<<endl<<"%DEM.PARTICLE.CIRCLE\n"<<n_particles_t<<endl;


    float area_total = 0;
    for (int i = 0; i < n_particles_t; ++i) {
        myfile2 <<i+1<<" 1"<<" "<<particles_r[i]<<" "<<particles_x[i]<<" "<<particles_y[i]<<" 0.00"<<endl;
        area_total = area_total + (float)M_PI*particles_r[i]*particles_r[i];
    }
    myfile2.close();
    grain_print(str_grain, N_granulometria, granulometrias, n_particles,  cell_PArea_vet,  index_aux, raios_r);
*/
    ofstream myfile;
    myfile.open ("time_record.grain", std::ios_base::app);
    float area_total = 0;
    for (int i = 0; i < n_particles_t; ++i) {
        area_total = area_total + (float)M_PI*particles_r[i]*particles_r[i];
    }
    float new_p = int(100*(area_total/(altura*largura)));

    myfile<<" Numero de particulas: "<< n_particles_t<< "N_int: "<< n_int<< "Tempo total: "<< timetotal<<endl;
    /*
    myfile<<new_p<<endl;
    myfile<<"time solve: "<<timer_solve/CLOCKS_PER_SEC<<endl;
    myfile<<"time sort: "<<timer_sort/CLOCKS_PER_SEC<<endl;
    myfile<<"time realloc: "<<timer_realloc/CLOCKS_PER_SEC<<endl;
    myfile<<"time remove: "<<timer_remove/CLOCKS_PER_SEC<<endl;
    myfile<<"time area: "<<timer_areac/CLOCKS_PER_SEC<<endl;
    myfile<<"time insert: "<<timer_insert/CLOCKS_PER_SEC<<endl;
    myfile<<"time cell_d: "<<timer_celld</CLOCKS_PER_SEC<endl;
    myfile<<"time part_d: "<<timer_partd/CLOCKS_PER_SEC<<endl;*/
    myfile.close();

    return int(1000 - 1000*(area_total/(altura*largura)));

}

//void realloc(float *particles, int n_particles_t, int n_cellsv, float max_r, vector< pair <int,int> > cells_id_par, vector< pair <int,int> > part_cont_par)
//{
//    cout<<part_cont_par[0].second<<endl;
//    for (int i = 0; i < 10; ++i) {
//        particles[4*n_particles_t + part_cont_par[i].second] = cells_id_par[i].second;
//        particles[n_particles_t + part_cont_par[i].second] = int(floor(cells_id_par[i].second/n_cellsv)*2*max_r) + (rand()%int(2*max_r));
//        particles[2*n_particles_t + part_cont_par[i].second] = int(floor(cells_id_par[i].second%n_cellsv)*2*max_r) + (rand()%int(2*max_r));
//    }
//}
//
//void p_remove(float *particles, int n_particles_t, int n_cellsv, float max_r, vector< pair <int,int> > cells_id_par, vector< pair <int,int> > part_cont_par)
//{
//    cout<<part_cont_par[0].second<<endl;
//    for (int i = 0; i < int(0.001*n_particles_t); ++i) {
//        particles[4*n_particles_t + part_cont_par[i].second] = n_cellsv-1;
//        particles[n_particles_t + part_cont_par[i].second] = 0;
//        particles[2*n_particles_t + part_cont_par[i].second] = 0;
//        particles[3*n_particles_t + part_cont_par[i].second] = 0;
//    }
//}