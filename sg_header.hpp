//
// Created by lucas_omena on 18/03/2019.
//


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "readfile.hpp"
#include <math.h>
#include "stdc++.h"


int cell_max(float *particles, int n_particles, int n_cells, float max_r, int n_cellsv)
{
    cout<<"teste1"<<endl;

    int cells[n_cells];
    int cell_now;
    cout<<"teste2"<<endl;


    for (int i = 0; i < n_cells; ++i)
    {
        cells[i] = 0;
    }
    cout<<"teste3"<<endl;

    for (int i = 0; i < n_particles; ++i)
    {
        cell_now = (int)(particles[i+2*n_particles]/(2*max_r)) + (int)(particles[i+n_particles]/(2*max_r))*n_cellsv;
        cells[cell_now]+=1;
    }
    cout<<"teste4"<<endl;

    int max_cell = 0;
    for (int i = 0; i < n_cells; ++i)
    {
        if (max_cell<cells[i])
        {
            max_cell = cells[i];
        }
    }
    cout<<"teste5"<<endl;

    return max_cell;
}


void particle_divide(float *particles, int n_cellsv, int idp, float  max_r, int max_cell, float altura,float largura)
{
    for (int id = 0; id < idp; ++id)
    {

        int cell_now;
        particles[id+idp] = particles[id+idp]+particles[5*idp+id];
        particles[id+2*idp] = particles[id+2*idp]+particles[6*idp+id];
        particles[5*idp+id] = 0;
        particles[6*idp+id] = 0;
        particles[7*idp] = 0;
        particles[7*idp+1] = 0;
        if ((particles[id+idp]>largura-particles[id+3*idp]))
        {
            particles[id+idp] = largura-particles[id+3*idp];
        }
        if ((particles[id+idp]<0+particles[id+3*idp]))
        {
            particles[id+idp] = 0+particles[id+3*idp];
        }
        if ((particles[id+2*idp]<0+particles[id+3*idp]))
        {
            particles[id+2*idp] = 0+particles[id+3*idp];
        }
        if ((particles[id+2*idp]>altura-particles[id+3*idp]))
        {
            particles[id+2*idp] = altura-particles[id+3*idp];
        }

        cell_now = (int)(particles[id+idp]/(2*max_r))*n_cellsv + (int)(particles[id+2*idp]/(2*max_r));
        particles[id+4*idp]= (float)(cell_now);

    }
}




void cell_divide(float* particles, int* cells, int max_cell, int idp )
{
    int cell_now;
    int n_particles_in_cell;
    int cell_id;
    for (int i = 0; i <idp; ++i)
    {
        cell_now  = (int)(particles[i+4*idp]);
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



void interactive_solve(float *particles,int *cells,float *max_overlap, int idp, int n_cellsv, int n_cells, int max_cell, int *max_sob_id , float *max_sob_value, float min_r  )
{
    for (int id = 0; id < idp; ++id)
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
                cell_now = (int)(particles[id+4*idp]) - n_cellsv -1+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    check = 1;
                }
            }
            else if (c<6)
            {
                cell_now = (int)(particles[id+4*idp]) -4+c;
                if ((cell_now>=0)&&(cell_now<n_cells))
                {
                    check = 1;
                }
            }
            else
            {
                cell_now = (int)(particles[id+4*idp])  - 7 + n_cellsv + c;
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
                        distx = fabs(particles[id+idp]-particles[j+idp]);
                        disty = fabs(particles[id+2*idp]-particles[j+2*idp]);
                        if (distx>0.01)
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
                        if(sqrt(distx*distx+disty*disty)<(particles[id+3*idp]+particles[j+3*idp]+0.005*(particles[id+3*idp]+particles[j+3*idp])))
                        {

                            sobx = ((particles[id+3*idp]+particles[j+3*idp]) - (sqrt(distx*distx+disty*disty)))*cosen;
                            soby = ((particles[id+3*idp]+particles[j+3*idp]) - (sqrt(distx*distx+disty*disty)))*sen;
                            if ((particles[id+idp]<particles[j+idp]))
                            {
                                particles[id+5*idp] = particles[id+5*idp]-sobx/2;
                            }
                            else if ((particles[id+idp]>particles[j+idp]))
                            {
                                particles[id+5*idp] = particles[id+5*idp]+sobx/2;

                            }
                            else
                            {
                                if ((particles[id]<particles[j]))
                                {
                                    particles[id+5*idp] = particles[id+5*idp]-sobx/2;
                                }
                                else
                                {
                                    particles[id+5*idp] = particles[id+5*idp]+sobx/2;

                                }

                            }

                            if ((particles[id+2*idp]<particles[j+2*idp]))
                            {
                                particles[id+6*idp] = particles[id+6*idp]-soby/2;
                            }
                            else if ((particles[id+2*idp]>particles[j+2*idp]))
                            {
                                particles[id+6*idp] = particles[id+6*idp]+soby/2;
                            }
                            else
                            {
                                if ((particles[id]<particles[j]))
                                {
                                    particles[id+6*idp] = particles[id+6*idp]-soby/2;
                                }
                                else
                                {
                                    particles[id+6*idp] = particles[id+6*idp]+soby/2;
                                }
                            }
                            if((particles[id+3*idp]+particles[j+3*idp]-sqrt(distx*distx+disty*disty))>max_sob_value[id])
                            {
                                max_sob_value[id] = (particles[id+3*idp]+particles[j+3*idp]-sqrt(distx*distx+disty*disty))/particles[id+3*idp];
                                max_sob_id[id] = id;
                            }


                            if((particles[id+3*idp]+particles[j+3*idp]-sqrt(distx*distx+disty*disty))>max_overlap[1])
                            {
                                max_overlap[1] = (particles[id+3*idp]+particles[j+3*idp]-sqrt(distx*distx+disty*disty))/particles[id+3*idp];
                                if (max_overlap[1]>0.05)
                                {
                                    max_overlap[0]=1;
                                }
                            }
                            particles[7*idp+1] = particles[7*idp+1]+particles[id+3*idp]+particles[j+3*idp]-sqrt(distx*distx+disty*disty);
                        }

                    }
                }
            }
        }
    }
}

//

void area_cell_compute(float *particles,  int n_cellsv,int idp, int n_cells, float *cell_PArea, int * id_smaller_particle, float max_r )

{
    int cell_now;
    float area_now;
    float * area_before = (float *)malloc(n_cells * sizeof(float));
    for (int k = 0; k < n_cells; ++k) {
        area_before[k] = (float)M_PI*max_r*max_r;
    }
    for (int i = 0; i < idp; ++i) {
        cell_now  = (int)(particles[i+4*idp]);
        area_now = (float)((particles[i+3*idp])*(particles[i+3*idp])*M_PI);
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


void realloc_area(float *particles, int idp, int n_cells, int n_cellsv, float max_r, vector< pair <float,int> > cell_PArea_vec, int * id_smaller_particle)
{
    int particle_now,cell_now, cell_after;
    int aux_count = 0;
    float r_now;
    for (int i = 0; i < int(0.01*n_cells+1); ++i) {

        cell_after = cell_PArea_vec[aux_count].second;
        cell_now = cell_PArea_vec[n_cells-aux_count-1].second;
        particle_now = id_smaller_particle[cell_now];
        r_now = particles[particle_now+3*idp];

        if (r_now<max_r/2)
        {
            particles[4*idp + particle_now] =cell_after;
            particles[idp + particle_now] = int(floor(cell_after/n_cellsv)*2*max_r) + (rand()%int(2*max_r));
            particles[2*idp +particle_now] = int(floor(cell_after%n_cellsv)*2*max_r) + (rand()%int(2*max_r));
        }
        else
        {
            i=i-1;
            aux_count = aux_count +1;
            continue;
        }
    }
}

void p_remove(float *particles, int idp, int n_cells, float max_r, vector< pair <float,int> > sob_vet, int* insert_index, float* insert_r)
{
    int particle_now;
    int aux_count = 0;
    float r_now;
    for (int i = 0; i < int(0.01*n_cells+1); ++i) {

        particle_now = sob_vet[idp- aux_count -1].second;
        r_now = particles[particle_now+3*idp];


        if (r_now<max_r/2)
        {
            insert_index[0] = insert_index[0] +1;
            insert_index[insert_index[0]] = int(particles[particle_now]);
            insert_r[0] = insert_r[0]+1;
            insert_r[int(insert_r[0])] = int(particles[particle_now+3*idp]);
            particles[4*idp + particle_now] = 0;
            particles[idp + particle_now] = 0;
            particles[2*idp +particle_now] = 0;
            particles[3*idp +particle_now] = 0;
        }
        else
        {
            i=i-1;
            aux_count = aux_count +1;
            continue;
        }



    }
}

void p_insert(float *particles, int idp, int n_cellsv, float max_r, vector< pair <float,int> > cell_PArea_vec, int* insert_index, float* insert_r)
{
    int particle_now;
    int  cell_after;
    int aux_count = 0;
    float r_now;

    for (int i = 1; i < (insert_index[0]+1); ++i) {

        cell_after = cell_PArea_vec[aux_count].second;
        particle_now = insert_index[i];
        r_now = insert_r[i];

        particles[4*idp + particle_now] =cell_after;
        particles[idp + particle_now] = int(floor(cell_after/n_cellsv)*2*max_r) + (rand()%int(2*max_r));
        particles[2*idp +particle_now] = int(floor(cell_after%n_cellsv)*2*max_r) + (rand()%int(2*max_r));
        particles[3*idp +particle_now] = r_now;


    }

}
