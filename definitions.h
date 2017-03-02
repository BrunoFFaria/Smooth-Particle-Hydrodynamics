//
// Created by bruno on 1/29/17.
//

#ifndef SPHSIMULATION_DEFINITIONS_H
    #define SPHSIMULATION_DEFINITIONS_H

    #include <stdbool.h>
    typedef float real_p;

    typedef struct{
        real_p * x;          /* components of the position */
        real_p * y;
        real_p * z;
    }f_t;

    typedef struct{
        real_p xmin;
        real_p xmax;
        real_p ymin;
        real_p ymax;
        real_p zmin;
        real_p zmax;
    }pol_t;



    typedef struct{
        real_p x;
        real_p y;
        real_p z;
    }v_t;

    typedef struct list{                        /* defines a list container general purpose */
        int value;
        int key;
    }list_t;

    typedef struct vector{
        list_t * lp;
        int num_items;
    }vector_t;

    typedef struct{
        list_t  * neighbours;
        real_p ** rji;         /* table with the distances (to avoid unnecessary calculations)*/

        f_t  * r;
        f_t  * v;
        f_t  * n;           /* surface normal: required for fluid visualization */
        f_t  * dv_dt;       /* */
        f_t  * dvis_dt;     /*  */
        f_t  * xsph_v;      /* velocity correction */
        real_p * rho;       /* density */
        real_p * drho_dt;

        real_p * mass;      /* mass of the ith particle */

        real_p * P;         /* pressure at particle i*/

        real_p * k;          /* stiffness */
        real_p * h;          /* variable smoothing length */

        real_p * c;          /* color field used for the martching cubes algorithm */

        real_p default_h;    /* default smoothing length */
        real_p rho_0;
        real_p mu_0;
        real_p elastic_coeff;
        real_p dt;
        real_p surface_threshold;
        real_p surface_tension;
        real_p xsph_corr;
        int     num_particles;
        int     num_bound_particles;
        int hash_table_size;
        int * march_cubes_vis_state;    /* holds the particles that have already been visited */
        int * inside_outside;           /* should we consider this particle for rendering ? */
    }parts_t;

    typedef struct{
        int x_particles;
        real_p Vwater;
        real_p rho_0;
        real_p elastic_coeff;
        real_p dt;
        real_p mu;
        real_p stiffness;
        real_p surface_threshold;
        real_p surface_tension;
        real_p xsph_corr;
        pol_t  (*define_system_init_state)(parts_t *);
    }params_t;

#define PI  3.141592653589793

typedef enum { EVAL_FUNCT, EVAL_DERIV, EVAL_DERIV2 }kern_t;

#endif //SPHSIMULATION_DEFINITIONS_H
