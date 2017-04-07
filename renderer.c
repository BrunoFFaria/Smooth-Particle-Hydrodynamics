//
// Created by bruno on 3/19/17.
//

#include "renderer.h"
#include "shaders.h"
#include <math.h>
#include "camera.h"
#include "definitions.h"

GLuint compile_shader(const char * vsource, const char * fsource){
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint frament_shader = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint program = 0;
    GLint success = 0;
    char temp[256];
    GLERROR;
    glShaderSource(vertex_shader, 1, &vsource, 0);
    glShaderSource(frament_shader, 1, &fsource, 0);

    glCompileShader(vertex_shader);
    glCompileShader(frament_shader);

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, frament_shader);
    glLinkProgram(program);

    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if(!success){
        glGetProgramInfoLog(program, 256, 0, temp);
        printf("Failed to link program:\n%s\n",temp);
        glDeleteProgram(program);
        program = 0;
    }
    glDeleteShader(vertex_shader);
    glDeleteShader(frament_shader);
    return program;
}

/* init shader (default) */
shader_t init_shader(const char * vsource, const char * fsource, shader_type_t shader_type, int width, int height){
    shader_t shader;

    /* start by compiling the shader */
    shader.program = compile_shader(vsource, fsource);

    glGenBuffers(1, &shader.ebo);
    glGenBuffers(1, &shader.vbo);
    glGenVertexArrays(1, &shader.vao);

    /* initialize textures */
    switch(shader_type){
        case DEPTH:
            /* initialize framebuffer */
            glGenFramebuffers(1,&shader.fbo);
            glBindFramebuffer(GL_FRAMEBUFFER, shader.fbo);

            glGenTextures(1, &shader.tex);
            glBindTexture(GL_TEXTURE_2D, shader.tex);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
            glBindTexture(GL_TEXTURE_2D, 0);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shader.tex, 0);
            glBindFramebuffer(GL_FRAMEBUFFER,0);
             break;
        case PLANE:
        case THICKNESS:
        case FLUIDFINAL:
        case FINAL:
            /* initialize framebuffer */
            glGenFramebuffers(1,&shader.fbo);
            glBindFramebuffer(GL_FRAMEBUFFER, shader.fbo);

            glGenTextures(1, &shader.tex);
            glBindTexture(GL_TEXTURE_2D, shader.tex);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
            glBindTexture(GL_TEXTURE_2D, 0);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, shader.tex, 0);
            glBindFramebuffer(GL_FRAMEBUFFER,0);
            break;
        case BLUR:


            /* initialize vertical framebuffer */
            glGenFramebuffers(1,&shader.fboV);
            glBindFramebuffer(GL_FRAMEBUFFER, shader.fboV);

            glGenTextures(1, &shader.texV);
            glBindTexture(GL_TEXTURE_2D, shader.texV);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

            glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
            glBindTexture(GL_TEXTURE_2D, 0);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shader.texV, 0);
            glBindFramebuffer(GL_FRAMEBUFFER,0);

            /* initialize horizontal framebuffer */
            glGenFramebuffers(1,&shader.fboH);
            glBindFramebuffer(GL_FRAMEBUFFER, shader.fboH);

            glGenTextures(1, &shader.texH);
            glBindTexture(GL_TEXTURE_2D, shader.texH);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

            glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
            glBindTexture(GL_TEXTURE_2D, 0);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shader.texH, 0);
            glBindFramebuffer(GL_FRAMEBUFFER,0);
            break;

        default:
            break;
    }
    GLint status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    /* check framebuffer status */
    if ( status != GL_FRAMEBUFFER_COMPLETE) {
        printf("frame buffer is not complete: %d\n", shader_type);
    }


    return shader;
}
void destroy_shader(shader_t shader, shader_type_t shader_type){
    /* start by destroying the shader */
    glDeleteProgram(shader.program);
    /* clear buffers */
    glDeleteBuffers(1, &shader.ebo);
    glDeleteBuffers(1, &shader.vbo);
    glDeleteBuffers(1, &shader.vao);
    switch(shader_type){
        case DEPTH:
        case PLANE:
        case THICKNESS:
        case FLUIDFINAL:
        case FINAL:
            glDeleteFramebuffers(1,&shader.fbo);
            glDeleteTextures(1,&shader.tex);
            break;
        case BLUR:
            glDeleteFramebuffers(1,&shader.fboH);
            glDeleteTextures(1,&shader.texH);
            glDeleteFramebuffers(1,&shader.fboV);
            glDeleteTextures(1,&shader.texV);
            break;
        default:
            break;
    }

}
void destroy_renderer(renderer_t renderer){

    destroy_shader(renderer.plane_shader, PLANE);
    destroy_shader(renderer.depth_shader, DEPTH);
    destroy_shader(renderer.blur_shader, BLUR);
    destroy_shader(renderer.thickness_shader, THICKNESS);
    destroy_shader(renderer.fluidfinal_shader, FLUIDFINAL);
    destroy_shader(renderer.final_shader, FINAL);
    glDeleteBuffers(1, &renderer.position_vbo);
}

renderer_t init_water_renderer(int width, int height, int num_particles){
    renderer_t renderer;

    /* start by initializing the shaders */
    renderer.plane_shader       = init_shader(plane_vertex_shader,        plane_fragment_shader,       PLANE,      width, height); /* depth     shader */
    renderer.depth_shader       = init_shader(water_depth_vertex_shader,  water_depth_fragment_shader, DEPTH,      width, height); /* depth     shader */
    renderer.blur_shader        = init_shader(blur_vertex_shader,         blur_fragment_shader,        BLUR,       width, height); /* blur      shader */
    renderer.thickness_shader   = init_shader(water_depth_vertex_shader,  thickness_fragment_shader,   THICKNESS,  width, height); /* thickness shader */
    renderer.fluidfinal_shader  = init_shader(fluidfinal_vertex_shader,   fluidfinal_fragment_shader,  FLUIDFINAL, width, height); /* thickness shader */
    renderer.final_shader       = init_shader(final_vertex_shader,        final_fragment_shader,       FINAL,      width, height); /* thickness shader */

    /* init VBO */
    glGenBuffers(1, &renderer.position_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, renderer.position_vbo);
    glBufferData(GL_ARRAY_BUFFER, num_particles * 3 * sizeof(float), 0, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    return renderer;
}


void render_water(renderer_t renderer, int width, int height, double fov, double particle_radius, int num_particles, vertex_t * particle_vertices) {
    GLfloat plane_vertices[] = {  1.0f,  1.0f, -0.25f-(float)particle_radius/2,
                                  1.0f, -1.0f, -0.25f-(float)particle_radius/2,
                                 -1.0f, -1.0f, -0.25f-(float)particle_radius/2,
                                 -1.0f,  1.0f, -0.25f-(float)particle_radius/2};

    GLfloat vertices[] = {1.0f, 1.0f, 1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 1.0f};
    GLint indexes[] = {0, 1, 3, 1, 2, 3};


    float * view_matrix;
    float * projection_matrix;
    int i = 0, j = 0;

    view_matrix = camera_get_view();
    projection_matrix = camera_get_projection();

    // ---- Draw the plane ----

    glUseProgram(renderer.plane_shader.program);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, renderer.plane_shader.fbo);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindVertexArray(renderer.plane_shader.vao);

    glBindBuffer(GL_ARRAY_BUFFER, renderer.plane_shader.vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(plane_vertices), plane_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, renderer.plane_shader.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indexes), indexes, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glUniformMatrix4fv(glGetUniformLocation(renderer.plane_shader.program, "projection"), 1, GL_FALSE, projection_matrix);
    glUniformMatrix4fv(glGetUniformLocation(renderer.plane_shader.program, "mView"), 1, GL_FALSE, view_matrix);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


    // ---- start with depth buffer ----
    glUseProgram(renderer.depth_shader.program);
    glBindFramebuffer(GL_FRAMEBUFFER, renderer.depth_shader.fbo);
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);
    glClear(GL_DEPTH_BUFFER_BIT);

    // ---- bind position VBO ----
    glBindVertexArray(renderer.depth_shader.vao);
    glBindBuffer(GL_ARRAY_BUFFER, renderer.position_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_t) * num_particles, particle_vertices, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glUniform1f( glGetUniformLocation(renderer.depth_shader.program, "pointScale"),(float) (height / (tan(fov * 0.5 * PI/180.0))));
    glUniform1f( glGetUniformLocation(renderer.depth_shader.program, "pointRadius"), (float) particle_radius/2.0f );//0.0125f*0.8f
    glUniformMatrix4fv(glGetUniformLocation(renderer.depth_shader.program, "projection"), 1, GL_FALSE, projection_matrix);
    glUniformMatrix4fv(glGetUniformLocation(renderer.depth_shader.program, "mView"), 1, GL_FALSE, view_matrix);

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SPRITE);
    glDrawArrays(GL_POINTS, 0, num_particles);

    // ---- blurr it ----
    glUseProgram(renderer.blur_shader.program);

    // vertical blur
    glBindFramebuffer(GL_FRAMEBUFFER, renderer.blur_shader.fboV);
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);

    glClear(GL_DEPTH_BUFFER_BIT);

    // ---- SHADER VAO QUAD ----
    glBindVertexArray(renderer.blur_shader.vao);
    glBindBuffer(GL_ARRAY_BUFFER, renderer.blur_shader.vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, renderer.blur_shader.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indexes), indexes, GL_STATIC_DRAW);

    glVertexAttribPointer(0,2,GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
    glEnableVertexAttribArray(0);

    // ---- depth_map ----
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, renderer.depth_shader.tex);
    glUniform1i(glGetUniformLocation(renderer.blur_shader.program,"depthMap"), 0);
    glUniformMatrix4fv(glGetUniformLocation(renderer.blur_shader.program, "projection"), 1, GL_FALSE, projection_matrix);
    glUniform2f(glGetUniformLocation(renderer.blur_shader.program,"screenSize"), (float)width, (float)height);
    glUniform2f(glGetUniformLocation(renderer.blur_shader.program,"blurDir"), 0.0f, (1.0f/(float)height));
    glUniform1f(glGetUniformLocation(renderer.blur_shader.program,"filterRadius"), 3.0f);
    glUniform1f(glGetUniformLocation(renderer.blur_shader.program,"blurScale"), 0.01f);

    glEnable(GL_DEPTH_TEST);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // ---- horizontal blur ----
    glBindFramebuffer(GL_FRAMEBUFFER, renderer.blur_shader.fboH);
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);

    glClear(GL_DEPTH_BUFFER_BIT);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, renderer.blur_shader.texV);
    glUniform1i(glGetUniformLocation(renderer.blur_shader.program, "depthMap"), 0);
    glUniform2f(glGetUniformLocation(renderer.blur_shader.program,"screenSize"), (float)width, (float)height);
    glUniformMatrix4fv(glGetUniformLocation(renderer.blur_shader.program, "projection"), 1, GL_FALSE, projection_matrix);
    glUniform2f(glGetUniformLocation(renderer.blur_shader.program,"blurDir"), (1.0f/(float)width), 0.0f);
    glUniform1f(glGetUniformLocation(renderer.blur_shader.program,"filterRadius"), 3.0f);
    glUniform1f(glGetUniformLocation(renderer.blur_shader.program,"blurScale"), 0.01f);

    glEnable(GL_DEPTH_TEST);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glDisable(GL_DEPTH_TEST);

    // ---- Particle Thickness ----
    glUseProgram(renderer.thickness_shader.program);
    glBindFramebuffer(GL_FRAMEBUFFER, renderer.thickness_shader.fbo);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // ---- bind position VBO ----
    glBindVertexArray(renderer.thickness_shader.vao);
    glBindBuffer(GL_ARRAY_BUFFER, renderer.position_vbo);
    // ---- copy points to GPU ----
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_t) * num_particles, particle_vertices, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glUniformMatrix4fv(glGetUniformLocation(renderer.thickness_shader.program, "projection"), 1, GL_FALSE, projection_matrix);
    glUniformMatrix4fv(glGetUniformLocation(renderer.thickness_shader.program, "mView"), 1, GL_FALSE, view_matrix);
    glUniform1f( glGetUniformLocation(renderer.thickness_shader.program, "pointScale"), (float) (height / (tan(fov * 0.5 * PI/180.0))));
    glUniform1f( glGetUniformLocation(renderer.thickness_shader.program, "pointRadius"), (float)(particle_radius/2));//0.0125f*0.8f

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);
    glBlendEquation(GL_FUNC_ADD);
    glDepthMask(GL_FALSE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    glDrawArrays(GL_POINTS, 0, (GLsizei)num_particles);

    glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    // ---- Particle fluidFinal -----
    glUseProgram(renderer.fluidfinal_shader.program);
    glBindFramebuffer(GL_FRAMEBUFFER, renderer.fluidfinal_shader.fbo);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // ---- SHADER VAO QUAD ----
    glBindVertexArray(renderer.fluidfinal_shader.vao);
    glBindBuffer(GL_ARRAY_BUFFER, renderer.fluidfinal_shader.vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, renderer.fluidfinal_shader.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indexes), indexes, GL_STATIC_DRAW);

    glVertexAttribPointer(0,2,GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, renderer.blur_shader.texH);
    glUniform1i(glGetUniformLocation(renderer.fluidfinal_shader.program, "depthMap"), 0);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, renderer.thickness_shader.tex);
    glUniform1i(glGetUniformLocation(renderer.fluidfinal_shader.program, "thicknessMap"), 1);


    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, renderer.plane_shader.tex);
    glUniform1i(glGetUniformLocation(renderer.fluidfinal_shader.program, "sceneMap"), 2);

    glUniformMatrix4fv(glGetUniformLocation(renderer.fluidfinal_shader.program, "projection"), 1, GL_FALSE, projection_matrix);
    glUniformMatrix4fv(glGetUniformLocation(renderer.fluidfinal_shader.program, "mView"), 1, GL_FALSE, view_matrix);
    glUniform4f(glGetUniformLocation(renderer.fluidfinal_shader.program, "color"),0.275f,0.65f,0.85f,0.9f);
    //glUniform4f(glGetUniformLocation(renderer.fluidfinal_shader.program, "color"),0.1451f,0.4471f,0.6118f,0.9f);
    glUniform2f(glGetUniformLocation(renderer.fluidfinal_shader.program,"invTexScale"), (float)(1.0f/(float)(width)), (float)(1.0f/(float)(height)));

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glDisable(GL_DEPTH_TEST);

    // in this pass render everything into screen

    // final pass

    glUseProgram(renderer.final_shader.program);
    glBindFramebuffer(GL_FRAMEBUFFER,0);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindVertexArray(renderer.final_shader.vao);
    glBindBuffer(GL_ARRAY_BUFFER, renderer.final_shader.vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, renderer.final_shader.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indexes), indexes, GL_STATIC_DRAW);

    glVertexAttribPointer(0,2,GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, renderer.fluidfinal_shader.tex);
    glUniform1i(glGetUniformLocation(renderer.final_shader.program, "fluidMap"), 0);

    //glUniform2f(glGetUniformLocation(renderer.final_shader.program,"screenSize"), (float)width, (float)height);

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glUseProgram(0);
    glDisableVertexAttribArray(0);

}