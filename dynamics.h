//
// Created by bruno on 2/3/17.
//

#ifndef SPHSIMULATION_DYNAMICS_H_H
#define SPHSIMULATION_DYNAMICS_H_H

    #include "definitions.h"
    #include "marching_cubes.h"

    void initialize_integrator(parts_t * parts, pol_t * pol, params_t params);
    void integrator_step(parts_t * parts, pol_t * pol);
    vector_t query_hash_neighbours_table(parts_t * parts, pol_t * pol, real_p cube_length, real_p x, real_p y, real_p z);
    int compute_surface(parts_t * parts, pol_t * pol, real_p resolution, real_p iso_thresh, vertex_t * vertexes, vertex_t * normals);
    int compute_surface2(parts_t * parts, pol_t * pol, real_p resolution, real_p iso_thresh, vertex_t * vertexes, vertex_t * normals);
#endif //SPHSIMULATION_DYNAMICS_H_H
