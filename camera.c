//
// Created by bruno on 3/26/17.
//

#include "renderer.h"
#include "camera.h"
#include "definitions.h"
#include <stdio.h>
static int w_width;
static int w_height;
static int camera_button_state = 0;
static int camera_old_x = 0;
static int camera_old_y = 0;
static float camera_z_near = 0.1f;
static float camera_z_far  = 5.0f;
static float camera_fov = 60.0;
static float camera_trans[] = {0, 0, -3};
static float camera_rot[]   = {0, 0, 0};
static float camera_trans_lag[] = {0, 0, -3};
static float camera_rot_lag[] = {0, 0, 0};
static float camera_inertia = 1;
static float camera_model_view_matrix[16];
static float camera_projection_matrix[16];

void init_camera(int w,int h){
    glutReshapeFunc(camera_reshape_window);
    glutMouseFunc(camera_mouse);
    glutMotionFunc(camera_motion);
    camera_reshape_window(w, h);
}

void camera_set_pos(float x, float y, float z){
    camera_trans[0] = x;
    camera_trans[1] = y;
    camera_trans[2] = z;
    camera_trans_lag[0] = x;
    camera_trans_lag[1] = y;
    camera_trans_lag[2] = z;
}

void camera_set_rot(float alpha, float beta, float gamma){
    camera_rot[0] = alpha;
    camera_rot[1] = beta;
    camera_rot[2] = gamma;
    camera_rot_lag[0] = alpha;
    camera_rot_lag[1] = beta;
    camera_rot_lag[2] = gamma;
}

void camera_render_view(){
    int c = 0;
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    for (c = 0; c < 3; ++c)
    {
        camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]) * camera_inertia;
        camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]) * camera_inertia;
    }
    glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2]);
    glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
    glRotatef(camera_rot_lag[2], 0.0, 0.0, 1.0);

    glGetFloatv(GL_MODELVIEW_MATRIX,  camera_model_view_matrix);
    glGetFloatv(GL_PROJECTION_MATRIX, camera_projection_matrix);
}

float * camera_get_view(){
    return camera_model_view_matrix;
}

float * camera_get_projection(){
    return camera_projection_matrix;
}

float camera_get_fov(){
    return camera_fov;
}

void camera_mouse(int button, int state, int x, int y) {
    if (state == GLUT_DOWN) {
        switch (button) {
            case GLUT_LEFT_BUTTON:
                camera_button_state = 1; // rotation
                break;
            case GLUT_MIDDLE_BUTTON:
                camera_button_state = 2; // zoom
                break;
            case GLUT_RIGHT_BUTTON:
                camera_button_state = 3; // movement
                break;
        }
    }else if (state == GLUT_UP){
       camera_button_state = 0;
    }

    camera_old_x = x;
    camera_old_y = y;

    glutPostRedisplay();
}

void camera_reshape_window(int w, int h){
    extern renderer_t renderer;
    extern parts_t * particles;
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective(camera_fov,((float)(w) / (float) (h)), camera_z_near, camera_z_far);

    glMatrixMode( GL_MODELVIEW );
    glViewport(0, 0, w, h);
    w_width = w;
    w_height = h;

    /* restart renderer here */
    destroy_renderer(renderer);
    renderer = init_water_renderer(w, h, particles->num_particles);
}

void camera_motion(int x, int y){
    float dx = 0, dy = 0;
    dx = (float)(x - camera_old_x);
    dy = (float)(y - camera_old_y);

    if(camera_button_state == 2){ // zoom
        camera_trans[2] += (dy/100.0f) * 0.5 * fabs(camera_trans[2]);
    }else if(camera_button_state == 3){
        camera_trans[0] += dx / 100.0f;
        camera_trans[1] -= dy / 100.0f;
    }else if(camera_button_state == 1){
        camera_rot[0] += dy / 5.0f;
        camera_rot[2] += dx / 5.0f;
    }
    camera_old_x = x;
    camera_old_y = y;
    glutPostRedisplay();
}