//
// Created by bruno on 1/29/17.
//

#ifndef SPHSIMULATION_MEMORY_H
#define SPHSIMULATION_MEMORY_H
    #include <stdlib.h>
    #include <stdio.h>
    #include <stdbool.h>
    #include <string.h>
    #include "definitions.h"

    parts_t * create_parts_obj(int num_particles);
    void destroy_parts_obj(parts_t * parts);
    void * mem_alloc(size_t size, bool zero_set);
    void * mem_realloc(void * ptr, int new_sz);
    void free_mem(void * ptr);
#endif //SPHSIMULATION_MEMORY_H
