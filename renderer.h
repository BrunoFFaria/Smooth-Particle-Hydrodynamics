//
// Created by bruno on 3/19/17.
//

#ifndef SPHSIMULATION_RENDERER_H
#define SPHSIMULATION_RENDERER_H
    #include <GL/glew.h>
    #include <GL/glut.h>

    #include <stdio.h>
    #include "definitions.h"
    #define GLERROR                                                \
    {                                                              \
        GLenum code = glGetError();                                \
        while (code!=GL_NO_ERROR)                                  \
        {                                                          \
            printf("%s\n",(char *) gluErrorString(code));          \
                code = glGetError();                               \
        }                                                          \
    }
    typedef enum {PLANE,DEPTH, BLUR, THICKNESS, FLUIDFINAL, FINAL} shader_type_t ;
    typedef struct {
        GLuint fbo;
        GLuint tex;
        GLuint tex2;
        GLuint vao;
        GLuint vbo;
        GLuint ebo;
        /* blur */
        GLuint fboV;
        GLuint fboH;
        GLuint texV;
        GLuint texH;

        GLuint program;
    }shader_t;

    typedef struct{
        shader_t plane_shader;
        shader_t depth_shader;
        shader_t blur_shader;
        shader_t thickness_shader;
        shader_t fluidfinal_shader;
        /* used to cast shadows */
        shader_t shadow_depth_shader;
        shader_t shadow_map_shader;
        shader_t final_shader;

        GLuint   position_vbo;
    }renderer_t;

    GLuint compile_shader(const char * vsource, const char * fsource);
    shader_t init_shader(const char * vsource, const char * fsource, shader_type_t shader_type, int width, int height);
    renderer_t init_water_renderer(int width, int height, int num_particles);
    void render_water(renderer_t renderer, int width, int height, double _scale, double particle_radius, int num_particles, vertex_t * particle_vertices);
    void destroy_renderer(renderer_t renderer);
#endif //SPHSIMULATION_RENDERER_H
