//
// Created by bruno on 1/29/17.
//

#include <math.h>
#include "memory.h"
#include "definitions.h"


parts_t * create_parts_obj(int num_particles){
    int i = 0;
    parts_t * parts = NULL;
    parts = mem_alloc(sizeof(parts_t), true);

    /* compute hash table size */
    parts->hash_table_size = (num_particles);
    parts->neighbours = mem_alloc( num_particles * sizeof(list_t), true);
    parts->rji = mem_alloc(num_particles * sizeof(real_p*), true);
    parts->rji[0] = mem_alloc(num_particles * num_particles * sizeof(real_p), true);

    for(i = 1; i < num_particles; i++)
        parts->rji[i] = parts->rji[0] + i*num_particles;

    parts->r        = mem_alloc(sizeof(f_t), true);
    parts->v        = mem_alloc(sizeof(f_t), true);
    parts->n        = mem_alloc(sizeof(f_t), true);

    parts->dv_dt    = mem_alloc(sizeof(f_t), true);
    parts->dvis_dt  = mem_alloc(sizeof(f_t), true);
    parts->xsph_v   = mem_alloc(sizeof(f_t), true);

    /* reserve space for the positions */
    parts->r->x = mem_alloc(num_particles * sizeof(real_p), true);
    parts->r->y = mem_alloc(num_particles * sizeof(real_p), true);
    parts->r->z = mem_alloc(num_particles * sizeof(real_p), true);

    /* reserve space for thevelocities */
    parts->v->x = mem_alloc(num_particles * sizeof(real_p), true);
    parts->v->y = mem_alloc(num_particles * sizeof(real_p), true);
    parts->v->z = mem_alloc(num_particles * sizeof(real_p), true);


    /* reserve space for the fluid surface normal */
    parts->n->x = mem_alloc(num_particles * sizeof(real_p), true);
    parts->n->y = mem_alloc(num_particles * sizeof(real_p), true);
    parts->n->z = mem_alloc(num_particles * sizeof(real_p), true);

    /* reserve space for accelerations */
    parts->dv_dt->x = mem_alloc(num_particles * sizeof(real_p), true);
    parts->dv_dt->y = mem_alloc(num_particles * sizeof(real_p), true);
    parts->dv_dt->z = mem_alloc(num_particles * sizeof(real_p), true);

    /* reserve space for viscosity derivative */
    parts->dvis_dt->x = mem_alloc(num_particles * sizeof(real_p), true);
    parts->dvis_dt->y = mem_alloc(num_particles * sizeof(real_p), true);
    parts->dvis_dt->z = mem_alloc(num_particles * sizeof(real_p), true);

    parts->xsph_v->x = mem_alloc(num_particles * sizeof(real_p), true);
    parts->xsph_v->y = mem_alloc(num_particles * sizeof(real_p), true);
    parts->xsph_v->z = mem_alloc(num_particles * sizeof(real_p), true);

    /* reserve space for the density */
    parts->rho = mem_alloc(num_particles * sizeof(real_p), true);

    /* reserve space for the density time derivative */
    parts->drho_dt = mem_alloc(num_particles * sizeof(real_p), true);

    /* reserve space for the mass */
    parts->mass = mem_alloc(num_particles * sizeof(real_p), true);

    /* reserve space for the pressure */
    parts->P = mem_alloc(num_particles * sizeof(real_p), true);

    /* stiffness */
    parts->k = mem_alloc(num_particles * sizeof(real_p), true);

    /* smoothing length */
    parts->h = mem_alloc(num_particles * sizeof(real_p), true);

    /* color field */
    parts->c = mem_alloc(num_particles * sizeof(real_p), true);

    /* marching cubes draw particle state */
    parts->march_cubes_vis_state = mem_alloc(num_particles * sizeof(int), true);

    parts->inside_outside = mem_alloc(num_particles * sizeof(int), true);

    parts->num_particles = num_particles;

    return parts;
}



void destroy_parts_obj(parts_t * parts){

    free_mem(parts->inside_outside);

    /* particle visualization state (marching cubes) memory */
    free_mem(parts->march_cubes_vis_state);

    /* color field */
    free_mem(parts->c);

    /* smoothing length */
    free_mem( parts->h );

    /* stiffness */
    free_mem( parts->k );

    /* free pressure memory space */
    free_mem( parts->P );

    /* free mass memory space */
    free_mem( parts->mass );

    /* free density time derivative memory space */
    free_mem( parts->drho_dt );

    /* free density memory space */
    free_mem( parts->rho );

    free_mem( parts->xsph_v->z );
    free_mem( parts->xsph_v->y );
    free_mem( parts->xsph_v->x );

    /* free viscosity derivative memory space */
    free_mem( parts->dvis_dt->z );
    free_mem( parts->dvis_dt->y );
    free_mem( parts->dvis_dt->x );

    /* free accelerations memory space */
    free_mem( parts->dv_dt->z );
    free_mem( parts->dv_dt->y );
    free_mem( parts->dv_dt->x );

    /* free surface normal memory space */
    free_mem( parts->n->x );
    free_mem( parts->n->y );
    free_mem( parts->n->z );

    /* free velocities memory space */
    free_mem( parts->v->z );
    free_mem( parts->v->y );
    free_mem( parts->v->x );

    /* free positions memory space*/
    free_mem( parts->r->z );
    free_mem( parts->r->y );
    free_mem( parts->r->x );

    free_mem( parts->xsph_v );
    free_mem( parts->dvis_dt );
    free_mem( parts->dv_dt );

    free_mem( parts->n   );
    free_mem( parts->v );
    free_mem( parts->r );

    free_mem( parts->rji[0] );
    free_mem( parts->rji );

    free_mem(parts->neighbours);
    free_mem(parts);
}

/* # size bytes memory allocation and initialization function */
void * mem_alloc(size_t size, bool zero_set)
{
    void * p = malloc(size);
    if(p == NULL){
        printf("could not allocate memory\n");
        exit(0);
    }

    if(zero_set == true){
        memset(p, 0, size);
    }
    return p;
}

void * mem_realloc(void * ptr, int new_sz)
{
    void * new_ptr = realloc( ptr, new_sz );

    if(new_ptr == NULL)
    {
        printf("could not reallocate memory!\n");
        exit(0);
    }
    return new_ptr;
}


/* memory deallocation function */
void free_mem(void * ptr)
{
    free( ptr );
}