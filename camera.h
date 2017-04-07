//
// Created by bruno on 3/26/17.
//

#ifndef SPHSIMULATION_CAMERA_H
#define SPHSIMULATION_CAMERA_H
    #include <GL/glut.h>
    #include <math.h>
    void init_camera(int w,int h);
    void camera_set_pos(float x, float y, float z);
    void camera_set_rot(float alpha, float beta, float gamma);
    void camera_render_view();
    void camera_mouse(int button, int state, int x, int y);
    void camera_reshape_window(int w, int h);
    void camera_motion(int x, int y);
    float * camera_get_view();
    float * camera_get_projection();
    float camera_get_fov();
#endif //SPHSIMULATION_CAMERA_H
