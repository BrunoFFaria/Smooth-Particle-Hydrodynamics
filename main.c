
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <GL/glew.h>
#include <math.h>
#include "definitions.h"
#include "memory.h"
#include "dynamics.h"
#include "marching_cubes.h"
#include "zpr.h"
#include "shaders.h"



/* define window parameters */
#define WIDTH					800
#define	HEIGHT					800

/* define update parameters */
#define FRAMERATE				60.0f
#define TIMEWAIT			1/FRAMERATE

/* define simulation parameters*/
#define NUMPARTICLES          25000
#define XPARTICLES              20
#define WATERVOLUME            0.1
#define RESTDENSITY            998
#define ELASTICCOEFF           0.0
#define TIMESTEP               0.004
#define VISCOSITYCOEFF         3.5
#define STIFFNESS              3.0
#define SURFACETHRESHOLD     7.065
#define SURFACETENSION      0.0728
#define XSPHCORRECTION         0.0f
#define GLERROR                                                    \
    {                                                              \
        GLenum code = glGetError();                                \
        while (code!=GL_NO_ERROR)                                  \
        {                                                          \
            printf("%s\n",(char *) gluErrorString(code));          \
                code = glGetError();                               \
        }                                                          \
    }
/* define particles as a global array */
parts_t * particles = NULL;
pol_t polygon = { 4*0.43,4*0.43 + 4*0.43, 8*0.43f, 8*0.43f + 8*0.43f, 4*0.43, 4*0.43 + 4*0.43 };
pol_t pol = {-0.25f, 0.25f, -0.5f, 0.5f, -0.25f, 0.25};


vertex_t vertices[20 * 3 * NUMPARTICLES];
vertex_t normals[20 * 3 * NUMPARTICLES];
color_t colors[20 * 3 * NUMPARTICLES];
int num_vertices = NUMPARTICLES;

/* define functions prototypes */
void drawBoundaries(void);
void pick(GLint name);
void display(void);
static void idle(void);
void reshape(int w, int h);
void keyboard(unsigned char key, int x, int y);
pol_t particles_initial_positions(parts_t * parts);

GLuint mprogram;
GLuint m_vbo;
GLuint m_color_vbo;
GLuint compile_shader(const char * vsource, const char * fsource);
double max_z = 0.0f;
double min_z = 0.0f;

/* */
static const GLfloat afAmbientWhite [] = {0.75, 0.75, 0.75, 1.00};
static const GLfloat afAmbientRed   [] = {0.25, 0.00, 0.00, 1.00};
static const GLfloat afAmbientGreen [] = {0.00, 0.25, 0.00, 1.00};
static const GLfloat afAmbientBlue  [] = {0.00, 0.00, 0.25, 1.00};
static const GLfloat afDiffuseWhite [] = {0.75, 0.75, 0.75, 1.00};
static const GLfloat afDiffuseRed   [] = {0.75, 0.00, 0.00, 1.00};
static const GLfloat afDiffuseGreen [] = {0.00, 0.75, 0.00, 1.00};
static const GLfloat afDiffuseBlue  [] = {0.00, 0.00, 0.75, 1.00};
static const GLfloat afSpecularWhite[] = {1.00, 1.00, 1.00, 1.00};
static const GLfloat afSpecularRed  [] = {1.00, 0.25, 0.25, 1.00};
static const GLfloat afSpecularGreen[] = {0.25, 1.00, 0.25, 1.00};
static const GLfloat afSpecularBlue [] = {0.25, 0.25, 1.00, 1.00};

float l_interp(float a, float b, float t){
    return a+t*(b-a);
}

void colorRamp(float t, float * r){
    const int ncolors = 7;
    float c[8][3] = {{1.0,0.0,0.0}, {1.0,0.5,0.0}, {1.0,1.0,0.0}, {0.0,1.0,0.0},
                           {0.0,1.0,1.0}, {0.0,0.0,1.0}, {1.0,0.0,1.0} };
    t = t * (ncolors-1);
    int i = (int)t;
    float u = t-floorf(t);
    r[0] = l_interp(c[i][0], c[i+1][0], u);
    r[1] = l_interp(c[i][1], c[i+1][1], u);
    r[2] = l_interp(c[i][2], c[i+1][2], u);
}

extern int w_width = WIDTH;
extern int w_height = HEIGHT;

int main(int argc, char * argv[]) {
    float *data = NULL;
    int i = 0;
    params_t params;
    GLfloat afPropertiesAmbient [] = {0.75, 0.75, 0.50, 1.00};
    GLfloat afPropertiesDiffuse [] = {0.75, 0.75, 0.75, 1.00};
    GLfloat afPropertiesSpecular[] = {1.00, 1.00, 1.00, 1.00};
    GLfloat lightPos[] = {0.1, 0.1f, -1.0f, 0.0f};
    glewExperimental=GL_FALSE;
    GLenum Err;

    params.x_particles = XPARTICLES;
    params.Vwater = WATERVOLUME;
    params.rho_0 = RESTDENSITY;
    params.elastic_coeff = ELASTICCOEFF;
    params.dt = TIMESTEP;
    params.mu = VISCOSITYCOEFF;
    params.stiffness = STIFFNESS;
    params.surface_threshold = SURFACETHRESHOLD;
    params.surface_tension = SURFACETENSION;
    params.xsph_corr = XSPHCORRECTION;
    params.define_system_init_state = &particles_initial_positions;

    /* allocate memory for the particles */
    particles = create_parts_obj( NUMPARTICLES );

    /* initialize integrator */
    initialize_integrator(particles, &polygon, params);

    /* initialise glut */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(WIDTH,HEIGHT);
    glutCreateWindow("SPH 1.0");

    Err = glewInit();

    if(Err!=GLEW_OK){
        fprintf(stderr,"Glew failed to initialize: %s", glewGetErrorString(Err));
        exit(-1);
    }
    if(!glewIsSupported("GL_VERSION_4_0 GL_VERSION_1_5 GL_ARB_multitexture GL_ARB_vertex_buffer_object")){
        fprintf(stderr,"Required OpenGL extensions missing.\n");
        exit(-1);
    }

    glClearColor( 0.0, 0.0, 0.0, 1.0 );
    glClearDepth( 1.0 );

    /* Configure GLUT callback functions */
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    GLERROR;

    glEnable(GL_DEPTH_TEST);


    /* compile shaders */
    mprogram = compile_shader(vertex_shader, sphere_shader);

    glClampColorARB( GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
    glClampColorARB( GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);

    /* enable VBO (positions)*/

    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo );
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_t) * NUMPARTICLES, vertices, GL_DYNAMIC_DRAW);


    /* enable vbo colours */
    glGenBuffers(1, &m_color_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_color_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GL_FLOAT) * 4 * NUMPARTICLES, 0, GL_DYNAMIC_DRAW);

    /* fill color buffer */
    glBindBufferARB(GL_ARRAY_BUFFER, m_color_vbo);

    /* get write address */
    data = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for(i = 0; i < NUMPARTICLES; i++){

        colorRamp(i/(float)(NUMPARTICLES),data);
        data+=3;
        *data++=1.0f;
    }
    glUnmapBufferARB(GL_ARRAY_BUFFER);


    /* another method for drawing */
    /*
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    glNormalPointer(GL_FLOAT, 0, normals);
    glVertexPointer(3, GL_FLOAT, 0, vertices);
    */

    GLERROR;
    glLineWidth(1.4);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glViewport(0, 0, (GLint) WIDTH, (GLint) HEIGHT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    GLERROR;
    glMatrixMode(GL_MODELVIEW);
    gluPerspective(60, (float)(w_width)/(float)(w_height), -1, 10);
    gluLookAt(0.45,0.80,0.25,0,0,0,0,0,1);

    GLERROR;


    /* Configure ZPR module */
    zprInit();
    zprSelectionFunc(drawBoundaries);       /* Selection mode draw function */
    zprPickFunc(pick);                      /* Pick event client callback   */
    GLERROR;

    /* Enter GLUT event loop */
    glutMainLoop();

    /* Simulation finished */
    destroy_parts_obj( particles );

    glDeleteBuffers(1, &m_vbo);
    glDeleteBuffers(1, &m_color_vbo);
    return 0;
}
int add_layer(real_p * x, real_p * y, real_p * z, real_p xmin, real_p xmax, real_p ymin, real_p ymax, real_p zmin, real_p zmax, real_p step){
    real_p xpos  =0.0f, ypos = 0.0f, zpos = 0.0f;
    int ct = 0;
    for(xpos = xmin; xpos < xmax; xpos+=step){
        for(ypos = ymin; ypos < ymax; ypos+=step){
            for(zpos = zmin; zpos < zmax; zpos+=step){
                /* */
                x[ct] = xpos;
                y[ct] = ypos;
                z[ct] = zpos;
                ct++;
            }
        }
    }
    return ct;
}
pol_t particles_initial_positions(parts_t * parts){
    int  j = 0, num_bound_particles = 0, num_particles = 0;
    real_p frac_bound = 4.00f, frac_part = 2.0f;
    real_p x_start = 0.0f, y_start = 0.0f, z_start = 0.0f;
    real_p x_end = 0.0f, y_end = 0.0f, z_end = 0.0f;
    real_p xpos = 0.0f, ypos = 0.0f, zpos = 0.0f;

    pol_t pol;

    /* start by defining the polygon */
    /* initialize polygon */
    pol.xmin = 10 * parts->default_h;
    pol.xmax = pol.xmin + 22 * parts->default_h;

    pol.ymin = 20 * parts->default_h;
    pol.ymax = pol.ymin + 32 * parts->default_h;

    pol.zmin = 10 * parts->default_h;
    pol.zmax = pol.zmin + 22 * parts->default_h;

    frac_bound = 4.0f;

    /* now define particle positions */
    frac_part = 2.0f;

    /* lets define two blocks of water with num_particles/2 in each*/
    x_start = pol.xmin + parts->default_h; x_end = pol.xmax - 4.0f * (parts->default_h/frac_bound);
    y_start = pol.ymin + parts->default_h; y_end = y_start + 9  * parts->default_h/frac_part;
    z_start = pol.zmin + parts->default_h/20; z_end = z_start + 33 * parts->default_h/frac_part;

    xpos = x_start; ypos = y_start; zpos = z_start;

    /* start by defining particle positions */
    for(j = 0; j < parts->num_particles; j++){
        if( (j % 42) == 0 && j > 0 ){
            xpos = x_start;
            ypos += parts->default_h/frac_part;
        }

        if( (j % 378) == 0 && j > 0 ){
            xpos = x_start;
            ypos = y_start + 35*parts->default_h/frac_part;
        }

        if( (j % 756) == 0 && j > 0 ){
            xpos = x_start;
            ypos = y_start;
            zpos += parts->default_h/frac_part;
        }

        parts->r->x[j] = xpos +  0.0f*rand()/(real_p)(RAND_MAX);;
        parts->r->y[j] = ypos +  0.0f*rand()/(real_p)(RAND_MAX);;
        parts->r->z[j] = zpos +  0.0f*rand()/(real_p)(RAND_MAX);;
        parts->v->x[j] = 0.001f*rand()/(real_p)(RAND_MAX);
        parts->v->y[j] = 0.001f*rand()/(real_p)(RAND_MAX);
        parts->v->z[j] = 0.001f*rand()/(real_p)(RAND_MAX);
        xpos += parts->default_h/frac_part;
    }

    return pol;
}

/*
void reshape(int w, int h){
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (float)w/((float)h),-1,10.0 );
    glMatrixMode(GL_MODELVIEW);
    glViewport(0,0,w,h);
    w_width = w;
    w_height = h;
}
*/

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
    return program;
}


real_p mapminmax(real_p xmin, real_p xmax, real_p ymin, real_p ymax, real_p x){
    return (ymax-ymin)*(x-xmin)/(xmax-xmin) + ymin;
}

static double ftime(void) {
    struct timeval t;
    gettimeofday(&t, NULL);

    return 1.0*t.tv_sec + 1e-6*t.tv_usec;
}

double last_T;

static void idle(void)
{
    int i = 0;
    float color = 0.0;
    const double now_T = ftime();
    const double delta_T = now_T - last_T;

    extern double _scale;

    if(delta_T > TIMEWAIT )
    {
        integrator_step(particles, &polygon);
        /*
        num_vertices = compute_surface2(particles, &polygon, particles->default_h/4, 1.0f, vertices, normals);
        for(i = 0; i < num_vertices; i++){
            if(i == 0){
                max_z = vertices[i].z;
                min_z = vertices[i].z;
            }else{
                if(vertices[i].z>max_z){
                    max_z = vertices[i].z;
                }
                if(vertices[i].z<min_z){
                    min_z = vertices[i].z;
                }
            }
            vertices[i].x = (GLfloat)mapminmax(polygon.xmin, polygon.xmax, pol.xmin, pol.xmax, vertices[i].x);
            vertices[i].y = (GLfloat)mapminmax(polygon.ymin, polygon.ymax, pol.ymin, pol.ymax, vertices[i].y);
            vertices[i].z = (GLfloat)mapminmax(polygon.zmin, polygon.zmax, pol.zmin, pol.zmax, vertices[i].z);

        }
        */

        for(i = 0; i < particles->num_particles; i++) {
            vertices[i].x = (GLfloat)mapminmax(polygon.xmin, polygon.xmax, pol.xmin, pol.xmax, particles->r->x[i]);
            vertices[i].y = (GLfloat)mapminmax(polygon.ymin, polygon.ymax, pol.ymin, pol.ymax, particles->r->y[i]);
            vertices[i].z = (GLfloat)mapminmax(polygon.zmin, polygon.zmax, pol.zmin, pol.zmax, particles->r->z[i]);
            if(i == 0){
                max_z = vertices[i].z;
                min_z = vertices[i].z;
            }else{
                if(vertices[i].z>max_z){
                    max_z = vertices[i].z;
                }
                if(vertices[i].z<min_z){
                    min_z = vertices[i].z;
                }
            }
        }

        for(i = 0; i < particles->num_particles; i++){
            color = (float)((vertices[i].z - min_z)/(real_p)(max_z-min_z));
            if(isnan(color)){
                color = 0.0f;
            }
            colorRamp(color, (float *)&colors[i]);
            colors[i].w=1.0f;
        }

        printf("%f - %f\n",1/delta_T,_scale);
        last_T = now_T;
        glutPostRedisplay();
    }
}

void keyboard(unsigned char key, int x, int y){
    FILE * fp = NULL;
    int i = 0;
    if(key == (int)'w'){
        fp = fopen("data.txt","w+");
        if(fp == NULL){ exit(-1); }
        for(i = 0; i < NUMPARTICLES; i++){
            fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", particles->r->x[i],particles->r->y[i],particles->r->z[i],
                    particles->v->x[i],particles->v->y[i],particles->v->z[i],
                    particles->n->x[i],particles->n->y[i],particles->n->z[i]
            );
        }
        fclose(fp);
    }
}

/* draw Boundaries */
void drawBoundaries(void)
{
    int i = 0;
    extern double _scale;

    /* Name-stack manipulation for the purpose of
       selection hit processing when mouse button
       is pressed.  Names are ignored in normal
       OpenGL rendering mode.                    */
    /* Render animation */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    glEnable(GL_POINT_SPRITE_ARB);
    glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);

    glUseProgram(mprogram);

    glUniform1f( glGetUniformLocation(mprogram, "point_scale"), (float)(w_height/(tan(60.0f * 0.5 * PI/180.0))));
    glUniform1f( glGetUniformLocation(mprogram, "point_radius"), (float)(particles->default_h));//0.0125f*0.8f
    glUniform1f( glGetUniformLocation(mprogram, "far"), 0.9f);//0.0125f*0.8f
    glUniform1f( glGetUniformLocation(mprogram, "near"), 0.1f);//0.0125f*0.8f
    glColor3d(1,1,1);

    /* time to draw the points*/

    glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_t) * particles->num_particles, vertices, GL_DYNAMIC_DRAW);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);

    /* Enable colours */

    glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_color_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(color_t) * particles->num_particles, colors, GL_DYNAMIC_DRAW);
    glColorPointer(4, GL_FLOAT, 0, 0);
    glEnableClientState(GL_COLOR_ARRAY);

    glDrawArrays(GL_POINTS, 0, particles->num_particles);

    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);

    glUseProgram(0);
    glDisable(GL_POINT_SPRITE_ARB);
    glDisableClientState(GL_VERTEX_ARRAY);

    /*
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDrawArrays(GL_TRIANGLES,0, num_vertices);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    */
    #define axis(xe,ye,ze,xs,ys,zs)		glVertex3d(xs,ys,zs);	\
		                                glVertex3d(xe,ye,ze);

    // draw polygon boundaries
    glBegin( GL_LINES );
    glColor3f(1.0,1.0,1.0);

    axis(pol.xmin-particles->default_h/2,pol.ymax+particles->default_h/2,pol.zmin-particles->default_h/2, pol.xmin-particles->default_h/2, pol.ymin-particles->default_h/2, pol.zmin-particles->default_h/2);
    axis(pol.xmax+particles->default_h/2,pol.ymax+particles->default_h/2,pol.zmin-particles->default_h/2, pol.xmax+particles->default_h/2, pol.ymin-particles->default_h/2, pol.zmin-particles->default_h/2);

    axis(pol.xmin-particles->default_h/2,pol.ymax+particles->default_h/2,pol.zmax+particles->default_h/2, pol.xmin-particles->default_h/2, pol.ymin-particles->default_h/2, pol.zmax+particles->default_h/2);
    axis(pol.xmax+particles->default_h/2,pol.ymax+particles->default_h/2,pol.zmax+particles->default_h/2, pol.xmax+particles->default_h/2, pol.ymin-particles->default_h/2, pol.zmax+particles->default_h/2);

    axis(pol.xmax+particles->default_h/2,pol.ymin-particles->default_h/2,pol.zmin-particles->default_h/2, pol.xmin-particles->default_h/2, pol.ymin-particles->default_h/2, pol.zmin-particles->default_h/2);
    axis(pol.xmax+particles->default_h/2,pol.ymin-particles->default_h/2,pol.zmax+particles->default_h/2, pol.xmin-particles->default_h/2, pol.ymin-particles->default_h/2, pol.zmax+particles->default_h/2);

    axis(pol.xmax+particles->default_h/2,pol.ymax+particles->default_h/2,pol.zmin-particles->default_h/2, pol.xmin-particles->default_h/2, pol.ymax+particles->default_h/2, pol.zmin-particles->default_h/2);
    axis(pol.xmax+particles->default_h/2,pol.ymax+particles->default_h/2,pol.zmax+particles->default_h/2, pol.xmin-particles->default_h/2, pol.ymax+particles->default_h/2, pol.zmax+particles->default_h/2);

    axis(pol.xmin-particles->default_h/2,pol.ymin-particles->default_h/2,pol.zmax+particles->default_h/2, pol.xmin-particles->default_h/2, pol.ymin-particles->default_h/2, pol.zmin-particles->default_h/2);
    axis(pol.xmax+particles->default_h/2,pol.ymin-particles->default_h/2,pol.zmax+particles->default_h/2, pol.xmax+particles->default_h/2, pol.ymin-particles->default_h/2, pol.zmin-particles->default_h/2);

    axis(pol.xmin-particles->default_h/2,pol.ymax+particles->default_h/2,pol.zmax+particles->default_h/2, pol.xmin-particles->default_h/2, pol.ymax+particles->default_h/2, pol.zmin-particles->default_h/2);
    axis(pol.xmax+particles->default_h/2,pol.ymax+particles->default_h/2,pol.zmax+particles->default_h/2, pol.xmax+particles->default_h/2, pol.ymax+particles->default_h/2, pol.zmin-particles->default_h/2);

    glEnd();

}


/* Callback function for drawing */
void display(void)
{
    GLERROR;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    drawBoundaries();
    glutSwapBuffers();
    glFlush();
    GLERROR;
}

/* Callback function for pick-event handling from ZPR */
void pick(GLint name)
{
    fflush(stdout);
}
