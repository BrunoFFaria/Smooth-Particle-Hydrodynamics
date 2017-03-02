//
// Created by bruno on 1/29/17.
//

#include "definitions.h"
#include "dynamics.h.h"
#include "marching_cubes.h"

#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


real_p kernel_pressure(real_p r, real_p h, kern_t kernel_operation){

    real_p value = 0, h6 = h*h*h*h*h*h, hr=h-r;
    if(kernel_operation == EVAL_DERIV2) {
        if(r >= 0 && r <= h){
            value = (real_p)( -90/(PI*h6) ) * (1/r) * (h - r) * (h - 2 * r);
        }
    }else if(kernel_operation == EVAL_DERIV) {
        if(r >= 0 && r <= h){
            value = (real_p)( (-45*hr*hr)/(PI*h6) );
        }
    }else{
        if(r >= 0 && r <= h){
            value = (real_p)( (15*hr * hr * hr)/(PI*h6) );
        }
    }
    return value;
}

real_p kernel_viscosity(real_p r, real_p h, kern_t kernel_operation){

    real_p value = 0, h2 = h*h, h3 = h * h * h, r3 = r*r*r,h9=h3*h3, r2 = r*r;
    if(kernel_operation == EVAL_DERIV2) {
        if(r >= 0 && r <= h){
            value = (real_p)( (-45/(PI*h9)) * (h - r));
        }
    }else if(kernel_operation == EVAL_DERIV) {
        if(r >= 0 && r <= h){
            value = (real_p)( 15/(2*PI*h3 ) * (-1.5*(r/h3) + (2/h2) - 0.5*(h/r3)));
        }
    }else{
        if(r >= 0 && r <= h){
            value = (real_p)( 15/(2*PI*h3 ) * ( (-0.5*(r3/h3)) + (r2/h2) + 0.5*(h/r) - 1));
        }
    }
    return value;
}


real_p kernel_normal(real_p r, real_p h, kern_t kernel_operation){

    real_p value = 0, h3 = h*h*h, h9 = h3*h3*h3, h2 = h*h, r2 = r*r;
    if(kernel_operation == EVAL_DERIV2) {
        if(r >= 0 && r <= h){
            value = (real_p)( -(945 * (h2 - r2) * (3*h2 - 7 * r2))/(32*PI*h9) );
        }
    }else if(kernel_operation == EVAL_DERIV) {
        if(r >= 0 && r <= h){
            value = (real_p)( (-945*(h2 - r2) * (h2 - r2))/(32*PI*h9) );
        }
    }else{
        if(r >= 0 && r <= h){
            value = (real_p)( (315 * (h2 - r2) * (h2 - r2) * (h2 - r2))/(64*PI*h9)) ;
        }
    }
    return value;
}

real_p kernel_cubic(real_p r, real_p h, kern_t kernel_operation){
    real_p value = 0, h3 = h*h*h, x = r/h;
    if(kernel_operation == EVAL_DERIV2) {
        if(r >= 0 && r <= h){
            value = (real_p)(1);
        }
    }else if(kernel_operation == EVAL_DERIV){
        if(x >= 0 && x <= 1){
            value = (real_p)( (1/(PI*h3*h)) * ((2.25) * x * x - 3 * x));
        }else if(x > 1 && x < 2){
            value = (real_p)( (1/(PI*h3*h)) * (-0.75) * ((2 - x) * (2-x)));
        }

    }else{
        if(x >= 0 && x <= 1){
            value = (real_p)( (1/(PI*h3)) * (1 - 1.5f * x * x + 0.75f * x * x * x));
        }else if(x > 1 && x < 2){
            value = (real_p)( (1/(PI*h3)) * 0.25f * (2-x) * (2-x) * (2-x) );
        }
    }
    return value;
}

#define DOT(a, b, i, j) (a->x[i]*b->x[j] + a->y[i]*b->y[j] + a->z[i]*b->z[j])

#define NORM(a, b, i, j) sqrtf((a->x[i]-b->x[j]) * (a->x[i]-b->x[j]) + \
                               (a->y[i]-b->y[j]) * (a->y[i]-b->y[j]) + \
                               (a->z[i]-b->z[j]) * (a->z[i]-b->z[j]))



void acceleration_computations(parts_t * parts, pol_t * pol){
    int i = 0, j = 0, k = 0;
    real_p cell_x, cell_y, cell_z;
    real_p cell_x_st, cell_y_st, cell_z_st;
    real_p rji = 0.0f, rji_m_o = 0.0f, nabbla2 = 0.0f, nabbla = 0.0f,P_rho_i2 = 0.0f,P_rho_j2 = 0.0f, surf_mag = 0.0f, rho_i_mo = 0.0f;
    vector_t neighbour_list;

    /* start by computing the derivative of the density and compute the pressure for all particles */
    #pragma omp parallel for private(i, j,cell_x, cell_y, cell_z, cell_x_st, cell_y_st, cell_z_st, rji, rji_m_o, nabbla2, nabbla, P_rho_i2, P_rho_j2, surf_mag, neighbour_list) shared(parts, pol)
    for(j = 0; j < parts->num_particles; j++){
        //parts->drho_dt[j] = 0.0;
        parts->rho[j] = 0.0f;

        for(cell_x = -parts->default_h; cell_x <= parts->default_h; cell_x+=parts->default_h){
            for(cell_y = -parts->default_h; cell_y <= parts->default_h; cell_y+=parts->default_h){
                for(cell_z = -parts->default_h; cell_z <= parts->default_h; cell_z+=parts->default_h){

                    /* query neighbour cell from hash table */
                    cell_x_st = (real_p)floor( (parts->r->x[j] + cell_x) / parts->default_h);
                    cell_y_st = (real_p)floor( (parts->r->y[j] + cell_y) / parts->default_h);
                    cell_z_st = (real_p)floor( (parts->r->z[j] + cell_z) / parts->default_h);

                    neighbour_list = query_hash_neighbours_table(parts, pol, parts->default_h, cell_x_st, cell_y_st, cell_z_st);

                    if(neighbour_list.num_items == 0){
                        continue;
                    }

                    /* compute density with cells from the list */
                    for(k = 0; k < neighbour_list.num_items; k++){
                        i = neighbour_list.lp[k].key;
                        if(i != j){
                            /* compute norm the rji */
                            rji = NORM(parts->r, parts->r, j, i);
                            parts->rji[j][i] = rji;
                            if(rji > parts->default_h){ continue; }

                            //parts->drho_dt[j] += parts->mass[i] * DOT(parts->v, parts->r, j, i) *
                            //                     kernel_cubic(rji, parts->default_h, EVAL_DERIV)/rji;
                            parts->rho[j] += parts->mass[i] * kernel_normal(rji, parts->default_h, EVAL_FUNCT);
                       }
                    }
                }
            }

        }

        //parts->rho[j] +=  parts->dt *  parts->drho_dt[j];
        /* Compute new pressure */
        parts->P[j] = parts->k[j] * (parts->rho[j] - parts->rho_0);
        //parts->P[j] = parts->k[j] * parts->rho[j];

    }

    /* now compute the viscous term of the navier-stokes equation and
     * the impact of the pressure on paricles (acceleration)*/
    #pragma omp parallel for private(i, j,cell_x, cell_y, cell_z, cell_x_st, cell_y_st, cell_z_st, rji, rji_m_o, nabbla2, nabbla, P_rho_i2, P_rho_j2, surf_mag, neighbour_list) shared(parts, pol)
    for(j = 0; j < parts->num_particles; j++) {
        parts->dvis_dt->x[j] = 0.0f; parts->dvis_dt->y[j] = 0.0f; parts->dvis_dt->z[j] = 0.0f;
        parts->xsph_v->x[j] = 0.0f;  parts->xsph_v->y[j] = 0.0f;  parts->xsph_v->z[j] = 0.0f;
        parts->dv_dt->x[j] = 0.0f;   parts->dv_dt->y[j] = 0.0f;   parts->dv_dt->z[j] = 0.0f;
        parts->n->x[j] = 0.0f;       parts->n->y[j] = 0.0f;       parts->n->z[j] = 0.0f;
        parts->c[j] = 0.0f;

        for(cell_x = -parts->default_h; cell_x <= parts->default_h; cell_x+=parts->default_h) {
            for (cell_y = -parts->default_h; cell_y <= parts->default_h; cell_y += parts->default_h) {
                for (cell_z = -parts->default_h; cell_z <= parts->default_h; cell_z += parts->default_h) {

                    /* query neighbour cell from hash table */
                    cell_x_st = (real_p)floor( (parts->r->x[j] + cell_x) / parts->default_h);
                    cell_y_st = (real_p)floor( (parts->r->y[j] + cell_y) / parts->default_h);
                    cell_z_st = (real_p)floor( (parts->r->z[j] + cell_z) / parts->default_h);

                    neighbour_list = query_hash_neighbours_table(parts, pol, parts->default_h, cell_x_st, cell_y_st, cell_z_st);

                    if(neighbour_list.num_items == 0){
                        continue;
                    }

                    /* compute particle */
                    for(k = 0; k < neighbour_list.num_items; k++){
                        i = neighbour_list.lp[k].key;


                        if (i != j) {
                            /* compute norm the rji */
                            //rji = NORM(parts->r, parts->r, j, i);
                            rji = parts->rji[j][i];
                            if(rji > parts->default_h){
                                continue;
                            }
                            rji_m_o = 1/rji;
                            rho_i_mo = 1/parts->rho[i];
                            /* nabbla2 */
                            nabbla2 = kernel_viscosity(rji, parts->default_h, EVAL_DERIV2);

                            /* viscous term */
                            parts->dvis_dt->x[j] +=
                                    parts->mu_0 * parts->mass[i] * (parts->v->x[j] - parts->v->x[i]) * nabbla2 * rho_i_mo;

                            parts->dvis_dt->y[j] +=
                                    parts->mu_0 * parts->mass[i] * (parts->v->y[j] - parts->v->y[i]) * nabbla2 * rho_i_mo;
                            parts->dvis_dt->z[j] +=
                                    parts->mu_0 * parts->mass[i] * (parts->v->z[j] - parts->v->z[i]) * nabbla2 * rho_i_mo;

                            /* acceleration due to pressure term */

                            /* nabla term */
                            nabbla = kernel_pressure(rji, parts->default_h, EVAL_DERIV);

                            P_rho_i2 = parts->P[i] / (parts->rho[i] * parts->rho[i]);
                            P_rho_j2 = parts->P[j] / (parts->rho[j] * parts->rho[j]);

                            parts->dv_dt->x[j] +=
                                    parts->mass[i] * (P_rho_j2 + P_rho_i2) * nabbla *
                                    (parts->r->x[j] - parts->r->x[i]) * rji_m_o;
                            parts->dv_dt->y[j] +=
                                    parts->mass[i] * (P_rho_j2 + P_rho_i2) * nabbla *
                                    (parts->r->y[j] - parts->r->y[i]) * rji_m_o;
                            parts->dv_dt->z[j] +=
                                    parts->mass[i] * (P_rho_j2 + P_rho_i2) * nabbla *
                                    (parts->r->z[j] - parts->r->z[i]) * rji_m_o;

                            /* evaluate surface color even if I don't use it */
                            parts->n->x[j] += parts->mass[i] * (parts->r->x[j] - parts->r->x[i]) * kernel_normal(rji, parts->default_h, EVAL_DERIV) / (parts->rho[i] * rji);
                            parts->n->y[j] += parts->mass[i] * (parts->r->y[j] - parts->r->y[i]) * kernel_normal(rji, parts->default_h, EVAL_DERIV) / (parts->rho[i] * rji);
                            parts->n->z[j] += parts->mass[i] * (parts->r->z[j] - parts->r->z[i]) * kernel_normal(rji, parts->default_h, EVAL_DERIV) / (parts->rho[i] * rji);

                            parts->c[j] += parts->mass[i] * kernel_normal(rji, parts->default_h, EVAL_DERIV2) * rho_i_mo;

                            /* evaluate xsph correction */
                            parts->xsph_v->x[j] += parts->mass[i] * (parts->v->x[i] - parts->v->x[j]) * kernel_pressure(rji, parts->default_h, EVAL_FUNCT) / (0.5f * (parts->rho[i]+parts->rho[j]));
                            parts->xsph_v->y[j] += parts->mass[i] * (parts->v->y[i] - parts->v->y[j]) * kernel_pressure(rji, parts->default_h, EVAL_FUNCT) / (0.5f * (parts->rho[i]+parts->rho[j]));
                            parts->xsph_v->z[j] += parts->mass[i] * (parts->v->z[i] - parts->v->z[j]) * kernel_pressure(rji, parts->default_h, EVAL_FUNCT) / (0.5f * (parts->rho[i]+parts->rho[j]));
                        }
                    }
                }
            }
        }

        /* compute surface tension */
        surf_mag = DOT(parts->n,parts->n,j,j);
        if( surf_mag > (parts->surface_threshold * parts->surface_threshold)){
            surf_mag = (real_p)sqrt(surf_mag);
            parts->n->x[j] = -parts->surface_tension * parts->n->x[j] * parts->c[j] / surf_mag ;
            parts->n->y[j] = -parts->surface_tension * parts->n->y[j] * parts->c[j] / surf_mag ;
            parts->n->z[j] = -parts->surface_tension * parts->n->z[j] * parts->c[j] / surf_mag ;
        }else{
            parts->n->x[j] = 0.0f;
            parts->n->y[j] = 0.0f;
            parts->n->z[j] = 0.0f;
        }

        /* sum all force components to get dv_dt derivative */
        parts->dv_dt->x[j] = (real_p)( -parts->dv_dt->x[j] + (parts->n->x[j] + parts->dvis_dt->x[j] ) / parts->rho[j] + 0.00 );
        parts->dv_dt->y[j] = (real_p)( -parts->dv_dt->y[j] + (parts->n->y[j] + parts->dvis_dt->y[j] ) / parts->rho[j] + 0.00 );
        parts->dv_dt->z[j] = (real_p)( -parts->dv_dt->z[j] + (parts->n->z[j] + parts->dvis_dt->z[j] ) / parts->rho[j] - 9.81 );
    }
}


int reflect_particle(v_t * r, v_t  * v, v_t * a, v_t plane_point, v_t plane_norm, real_p elastic_coeff, real_p dt)
{
    real_p d = 0.0f, norm = 0.0f, d_real_space = 0.0f, l_dot_n = 0.0f, dot = 0;
    v_t  part_dir, contact_point;
    int reflected = 0;
    /* compute particle direction */
    part_dir.x = v->x; part_dir.y = v->y; part_dir.z = v->z;

    norm = (real_p)sqrt(part_dir.x * part_dir.x + part_dir.y * part_dir.y + part_dir.z * part_dir.z);

    if(norm == 0){ return reflected; }

    part_dir.x /= norm; part_dir.y /= norm; part_dir.z /= norm;

    d_real_space = (plane_point.x - r->x)  * plane_norm.x +
                   (plane_point.y - r->y)  * plane_norm.y +
                   (plane_point.z - r->z)  * plane_norm.z;

    l_dot_n = (part_dir.x * plane_norm.x + part_dir.y * plane_norm.y + part_dir.z * plane_norm.z);
    if(l_dot_n == 0.0f || d_real_space == 0.0f){ return reflected; }

    d_real_space /= l_dot_n;

    contact_point.x = d_real_space/l_dot_n * part_dir.x + r->x;
    contact_point.y = d_real_space/l_dot_n * part_dir.y + r->y;
    contact_point.z = d_real_space/l_dot_n * part_dir.z + r->z;

    /* compute distance travelled away from the contact point */
    d = (real_p)sqrt( (r->x - contact_point.x) * (r->x - contact_point.x) +
                      (r->y - contact_point.y) * (r->y - contact_point.y) +
                      (r->z - contact_point.z) * (r->z - contact_point.z));

    if(d_real_space*l_dot_n > 0.0) {

        /* compute velocity normal */
        //norm = ()sqrt(v->x * v->x + v->y * v->y + v->z * v->z);

        r->x += (d_real_space * l_dot_n + 0.00001 * (rand()/(double)RAND_MAX)) * plane_norm.x;
        r->y += (d_real_space * l_dot_n + 0.00001 * (rand()/(double)RAND_MAX)) * plane_norm.y;
        r->z += (d_real_space * l_dot_n + 0.00001 * (rand()/(double)RAND_MAX)) * plane_norm.z;
        dot = (v->x * plane_norm.x + v->y * plane_norm.y + v->z * plane_norm.z);
        v->x -= dot * (1+elastic_coeff) * plane_norm.x;
        v->y -= dot * (1+elastic_coeff) * plane_norm.y;
        v->z -= dot * (1+elastic_coeff) * plane_norm.z;

        reflected = 1;
    }

    return reflected;
}

void boundary_collision(parts_t * parts, pol_t * pol)
{
    int j = 0;

    v_t plane_point,  plane_norm, r, v, a;

    /* Perform collision detection with boundaries (walls, floor, etc) */
    /* for now I will use a simple polygon */
    #pragma omp parallel for private( j, r, v, a, plane_point, plane_norm)
    for (j = 0; j < parts->num_particles; j++) {

        r.x = parts->r->x[j];
        r.y = parts->r->y[j];
        r.z = parts->r->z[j];
        v.x = parts->v->x[j];
        v.y = parts->v->y[j];
        v.z = parts->v->z[j];
        a.x = parts->dv_dt->x[j];
        a.y = parts->dv_dt->y[j];
        a.z = parts->dv_dt->z[j];

        plane_point.x = pol->xmin;
        plane_point.y = (pol->ymin + 0.5f * (pol->ymax - pol->ymin))*1;
        plane_point.z = (pol->zmin + 0.5f * (pol->zmax - pol->zmin))*1;
        plane_norm.x = 1;
        plane_norm.y = 0;
        plane_norm.z = 0;

        reflect_particle(&r, &v, &a, plane_point, plane_norm, parts->elastic_coeff, parts->dt);

        plane_point.x = pol->xmax;
        plane_point.y = (pol->ymin + 0.5f * (pol->ymax - pol->ymin))*1;
        plane_point.z = (pol->zmin + 0.5f * (pol->zmax - pol->zmin))*1;
        plane_norm.x = -1;
        plane_norm.y = 0;
        plane_norm.z = 0;
        reflect_particle(&r, &v, &a, plane_point, plane_norm, parts->elastic_coeff, parts->dt);


        plane_point.x = (pol->xmin + 0.5f * (pol->xmax - pol->xmin))*1;
        plane_point.y = pol->ymin;
        plane_point.z = (pol->zmin + 0.5f * (pol->zmax - pol->zmin))*1;
        plane_norm.x = 0;
        plane_norm.y = 1;
        plane_norm.z = 0;
        reflect_particle(&r, &v, &a, plane_point, plane_norm, parts->elastic_coeff, parts->dt);

        plane_point.x = (pol->xmin + 0.5f * (pol->xmax - pol->xmin)) * 1;
        plane_point.y = pol->ymax;
        plane_point.z = (pol->zmin + 0.5f * (pol->zmax - pol->zmin)) * 1;
        plane_norm.x = 0;
        plane_norm.y = -1;
        plane_norm.z = 0;
        reflect_particle(&r, &v, &a, plane_point, plane_norm, parts->elastic_coeff, parts->dt);

        plane_point.x = (pol->xmin + 0.5f * (pol->xmax - pol->xmin))*1;
        plane_point.y = (pol->ymin + 0.5f * (pol->ymax - pol->ymin))*1;
        plane_point.z = pol->zmin;
        plane_norm.x = 0;
        plane_norm.y = 0;
        plane_norm.z = 1;
        reflect_particle(&r, &v, &a, plane_point, plane_norm, parts->elastic_coeff, parts->dt);


        plane_point.x = (pol->xmin + 0.5f * (pol->xmax - pol->xmin))*1;
        plane_point.y = (pol->ymin + 0.5f * (pol->ymax - pol->ymin))*1;
        plane_point.z = pol->zmax;
        plane_norm.x = 0;
        plane_norm.y = 0;
        plane_norm.z = -1;
        reflect_particle(&r, &v, &a, plane_point, plane_norm, parts->elastic_coeff, parts->dt);


        parts->r->x[j] = r.x;
        parts->r->y[j] = r.y;
        parts->r->z[j] = r.z;
        parts->v->x[j] = v.x;
        parts->v->y[j] = v.y;
        parts->v->z[j] = v.z;

        parts->dv_dt->x[j] = a.x;
        parts->dv_dt->y[j] = a.y;
        parts->dv_dt->z[j] = a.z;

    }
}




int compare_list(const void * a, const void *b)
{
    if( (*(list_t*)a).value < (*(list_t *)b).value){
        return -1;
    }else if( (*(list_t*)a).value == (*(list_t *)b).value){
        return  0;
    }else{
        return  1;
    }
}

void hash_neighbours(parts_t * parts, pol_t * pol, real_p cube_length){
    int j = 0;
    unsigned int p1 = 73856093, p2 = 19349663, p3 = 83492791;
    unsigned int x1 = 0, x2 = 0, x3 = 0, hash = 0;

    /* clear previous hash table */
    memset(parts->neighbours, 0, parts->num_particles * sizeof(list_t));

    /* fill new hash table information */
    #pragma omp parallel for private( j, x1, x2, x3, hash )
    for(j = 0; j < parts->hash_table_size; j++){

        x1 = (unsigned int) (floor(parts->r->x[j] / cube_length));
        x2 = (unsigned int) (floor(parts->r->y[j] / cube_length));
        x3 = (unsigned int) (floor(parts->r->z[j] / cube_length));

        parts->neighbours[j].key=j;
        hash = (unsigned int)(x3 * (2.0f*(pol->xmax - pol->xmin)/cube_length) * 2.0f *((pol->ymax - pol->ymin)/cube_length) +
                              x2 * 2.0f *((pol->xmax - pol->xmin)/cube_length) + x1);

        parts->neighbours[j].value = hash;
    }


    /* sort hash table */
    qsort(parts->neighbours, (size_t)parts->hash_table_size, sizeof(list_t), compare_list);
}

vector_t query_hash_neighbours_table(parts_t * parts, pol_t * pol, real_p cube_length, real_p x, real_p y, real_p z){

    unsigned int p1 = 73856093, p2 = 19349663, p3 = 83492791, hash = 0;
    int num_items = 0, i = 0;
    list_t * lt_ptr = NULL;
    vector_t vect;

    hash = (unsigned int)(z * (2.0f*(pol->xmax - pol->xmin)/cube_length) * 2.0f *((pol->ymax - pol->ymin)/cube_length) +
                          y * 2.0f *((pol->xmax - pol->xmin)/cube_length) + x);


    /* perform a binary search for this hash */
    lt_ptr = bsearch(&hash, parts->neighbours, (size_t)parts->hash_table_size, sizeof(list_t), compare_list );
    if(lt_ptr == NULL){
        num_items = 0;
        vect.lp = NULL;
        vect.num_items = 0;
    }else{
        /* search for the farthest left item */
        while(lt_ptr > parts->neighbours ){
            if(lt_ptr->value != hash){
                if( (lt_ptr + 1) < (parts->neighbours + parts->hash_table_size))
                lt_ptr++;
                break;
            }else{
                lt_ptr--;
            }
        }
        vect.lp=lt_ptr;

        /* search for the right most item */
        while(lt_ptr < (parts->neighbours + parts->hash_table_size)){
            if(lt_ptr->value != hash){
                break;
            }else{
                lt_ptr++;
            }
        }

        /* number of items */
        vect.num_items = (int)(lt_ptr - vect.lp);
    }

    return vect;
}

void initialize_integrator(parts_t * parts, pol_t * pol, params_t params)
{
    real_p default_mass = 0.0f;
    real_p default_h = 0.0f, rji;
    int j = 0, i = 0;
    default_mass = params.rho_0 * params.Vwater/(real_p)(parts->num_particles);
    default_h = pow((3 * params.Vwater * (real_p)params.x_particles / (4*PI*(real_p)(parts->num_particles))),1.0f/3.0f);

    parts->default_h = default_h;
    parts->rho_0 = params.rho_0;
    parts->elastic_coeff = params.elastic_coeff;
    parts->dt = params.dt;
    parts->mu_0 = params.mu;
    parts->surface_tension = params.surface_tension;
    parts->surface_threshold = params.surface_threshold;
    parts->xsph_corr = params.xsph_corr;

    for( j = 0; j < parts->num_particles; j++){
        parts->mass[ j ] = default_mass;
        parts->k[ j ] = params.stiffness;
        parts->rho[ j ] = params.rho_0;
        parts->P[ j ] = 0;
    }

    /* define initial positions */
    (*pol) = params.define_system_init_state(parts);

}

real_p get_iso_value(parts_t * parts, real_p x0, real_p y0, real_p z0){
    int i = 0, j = 0, k = 0;
    real_p r = 0.0f, rsqr = 0.0f, iso = 0.0f, a = parts->default_h/3, b = parts->default_h;

    for(j = 0; j < parts->num_particles; j++){
        r = (x0 - parts->r->x[j]) * (x0 - parts->r->x[j]) +
            (y0 - parts->r->y[j]) * (y0 - parts->r->y[j]) +
            (z0 - parts->r->z[j]) * (z0 - parts->r->z[j]);

            iso += ((parts->default_h/20) * (parts->default_h/20)) / (r);
    }

    return iso;
}

real_p get_iso_interp(real_p val1, real_p val2, real_p iso_thresh){
    if(val2 - val1  == 0.0){
        return 0.5f;
    }else{
        return  (iso_thresh - val1)/(val2-val1);
    }
}
xyz get_normal(parts_t * parts, vector_t * neighbours_list, real_p x, real_p y, real_p z){
    xyz n;
    real_p norm = 0.0f;
    n.x = get_iso_value(parts,  x - 0.01, y, z) - get_iso_value(parts,  x + 0.01, y, z);
    n.y = get_iso_value(parts,  x, y - 0.01, z) - get_iso_value(parts,  x, y + 0.01, z);
    n.z = get_iso_value(parts, x, y, z - 0.01) - get_iso_value(parts,  x, y, z + 0.01);

    norm = sqrt(n.x * n.x + n.y*n.y + n.z * n.z);
    n.x /= norm;
    n.y /= norm;
    n.z /= norm;
    return n;
}




int compute_surface2(parts_t * parts, pol_t * pol, real_p resolution, real_p iso_thresh, vertex_t * vertexes, vertex_t * normals)
{

    real_p cell_x, cell_y, cell_z;
    real_p interped_iso = 0.0f;

    xyz vertex_list[12];
    xyz normal_list[12];

    extern char  tri_table[256][16];
    extern int   edge_table[256];
    extern int   edge_conn[12][2];
    extern float edge_dir[12][3];

    int j = 0, i = 0, k = 0, cube_index = 0, edges_flags = 0, ivertex = 0;
    vector_t neighbours_list[27];
    grid_cell_t grid;
    triangles_t triangles[5];
    int ntrig = 0, num_total_triangles = 0;
    real_p neighbours[8][3] = { {0,0,0},{resolution,0,0},
                                { resolution, resolution,0},{0,resolution,0},
                                {0,0,resolution},{resolution,0,resolution},
                                { resolution, resolution,resolution},{0,resolution,resolution}};

    num_total_triangles = 0;

    /* first clear visited state, then start drawing */
    memset(parts->march_cubes_vis_state, 0, parts->num_particles * sizeof(int));


    for (cell_x = pol->xmin; cell_x <= pol->xmax; cell_x += resolution) {
        for (cell_y = pol->ymin; cell_y <= pol->ymax; cell_y += resolution) {
            for (cell_z = pol->zmin; cell_z <= pol->zmax; cell_z += resolution) {

                /* get isovalues */
                for (i = 0; i < 8; i++) {
                    grid.p[i].x = cell_x + neighbours[i][0];
                    grid.p[i].y = cell_y + neighbours[i][1];
                    grid.p[i].z = cell_z + neighbours[i][2];
                    grid.val[i] = get_iso_value(parts, grid.p[i].x, grid.p[i].y, grid.p[i].z);

                    if (grid.val[i] < iso_thresh) {
                        cube_index |= 1 << i;
                    }
                }
                edges_flags = edge_table[cube_index];
                /* check whether the cube is inside or outside */

                /* if totally inside or outside move on to the next one */
                if (edges_flags == 0) { continue; }

                /* go over each edge and interpolate */
                for (i = 0; i < 12; i++) {

                    /* is there any intersection? */
                    if (edges_flags & (1 << i)) {

                        /* interpolate the iso value on this edge */
                        vertex_list[i] = VertexInterp(iso_thresh, grid.p[edge_conn[i][0]], grid.p[edge_conn[i][1]],
                                                      grid.val[edge_conn[i][0]], grid.val[edge_conn[i][1]]);

                        /* now get the normal at this vertex */
                        normal_list[i] = get_normal(parts, neighbours_list, vertex_list[i].x, vertex_list[i].y,
                                                    vertex_list[i].z);
                    }
                }
                /* add found the triangles and normals to the list */
                for(i=0;tri_table[cube_index][i]!=-1;i+=3) {
                    vertexes[num_total_triangles].x = (GLfloat)vertex_list[tri_table[cube_index][i  ]].x;
                    vertexes[num_total_triangles].y = (GLfloat)vertex_list[tri_table[cube_index][i+1]].y;
                    vertexes[num_total_triangles].z = (GLfloat)vertex_list[tri_table[cube_index][i+2]].z;

                    /* compute surface normal */
                    normals[num_total_triangles].x = (GLfloat)normal_list[tri_table[cube_index][i  ]].x;
                    normals[num_total_triangles].y = (GLfloat)normal_list[tri_table[cube_index][i+1]].y;
                    normals[num_total_triangles].z = (GLfloat)normal_list[tri_table[cube_index][i+2]].z;

                    num_total_triangles++;
                }
            }
        }
    }

    return num_total_triangles;
}



void integrator_step(parts_t * parts, pol_t * pol)
{
    int j = 0;
    v_t new_pos,new_vel;

    /* start by hashing neighbours of the integer step (positions and velocities used for force computation) */
    hash_neighbours(parts, pol, parts->default_h);

    /* start by computing the forces at t+1/2 using r and v at time t */
    acceleration_computations( parts, pol );
    #pragma omp parallel for private( j, new_pos, new_vel )
    for(j = 0; j < parts->num_particles; j++){
        /* compute new positions at the next time step */
        new_pos.x = parts->r->x[j] + (parts->v->x[j] + parts->xsph_corr * parts->xsph_v->x[j]) * parts->dt + 0.5f * parts->dv_dt->x[j] * parts->dt * parts->dt;
        new_pos.y = parts->r->y[j] + (parts->v->y[j] + parts->xsph_corr * parts->xsph_v->y[j]) * parts->dt + 0.5f * parts->dv_dt->y[j] * parts->dt * parts->dt;
        new_pos.z = parts->r->z[j] + (parts->v->z[j] + parts->xsph_corr * parts->xsph_v->z[j]) * parts->dt + 0.5f * parts->dv_dt->z[j] * parts->dt * parts->dt;

        /* compute the velocity at t+1 */
        new_vel.x = (new_pos.x - parts->r->x[j]) / (parts->dt);
        new_vel.y = (new_pos.y - parts->r->y[j]) / (parts->dt);
        new_vel.z = (new_pos.z - parts->r->z[j]) / (parts->dt);

        parts->r->x[j] = new_pos.x;
        parts->r->y[j] = new_pos.y;
        parts->r->z[j] = new_pos.z;
        parts->v->x[j] = new_vel.x;
        parts->v->y[j] = new_vel.y;
        parts->v->z[j] = new_vel.z;
    }

    /* check boundary conditions */
    boundary_collision( parts, pol);

}
