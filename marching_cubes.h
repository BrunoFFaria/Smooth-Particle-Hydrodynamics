//
// Created by bruno on 2/9/17.
//

#ifndef SPHSIMULATION_MARCHING_CUBES_H
#define SPHSIMULATION_MARCHING_CUBES_H
    #include <GL/glut.h>

    typedef struct {
        double x,y,z;
    }xyz;

    typedef struct{
        GLfloat x,y,z;
    }vertex_t;

    typedef struct{
        GLfloat r,g,b,w;
    }color_t;

    typedef struct {
        xyz p[8];		//position of each corner of the grid in world space
        double val[8];	//value of the function at this grid corner
        xyz n[8];       // normal at this position
    }grid_cell_t;

    typedef struct {
        xyz p[3];         /* Vertices */
        xyz c;            /* Centroid */
        xyz n[3];         /* Normal   */
    }triangles_t;

    int Polygonise(grid_cell_t * grid, triangles_t * triangles, double iso);
    xyz VertexInterp(double isolevel, xyz p1,xyz p2,double valp1,double valp2);

#endif //SPHSIMULATION_MARCHING_CUBES_H
