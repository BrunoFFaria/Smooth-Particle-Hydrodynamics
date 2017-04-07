
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "renderer.h"

#include "definitions.h"
#include "memory.h"
#include "dynamics.h"
//#include "zpr.h"
#include "camera.h"
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
#define ELASTICCOEFF           1
#define TIMESTEP               0.0045
#define VISCOSITYCOEFF         3.5
#define STIFFNESS              3.0
#define SURFACETHRESHOLD     7.065
#define SURFACETENSION      0.0728
#define XSPHCORRECTION         0.5f

/* define particles as a global array */
parts_t * particles = NULL;
pol_t polygon = { 4*0.43f,4*0.43f + 4*0.43f, 8*0.43f, 8*0.43f + 8*0.43f, 4*0.43f, 4*0.43f + 4*0.43f};
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
void print(float x, float y, const char *string);

GLuint depth_program;


/* VBOs */
GLuint position_vbo;
GLuint color_vbo;

double max_z = 0.0f;
double min_z = 0.0f;

/* */
int visualization_mode = 1;
renderer_t renderer;

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

    glClearColor( 0.4, 0.4, 0.4, 1.0 );
    glClearDepth( 1.0 );

    /* Configure GLUT callback functions */
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    GLERROR;

    glEnable(GL_DEPTH_TEST);

    /* compile shaders */
    depth_program = compile_shader(depth_vertex_shader, depth_fragment_shader);

    glClampColor(GL_CLAMP_VERTEX_COLOR, GL_FALSE);
    glClampColor(GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);

    /* enable VBOs (positions)*/
    glGenBuffers(1, &position_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, position_vbo );
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_t) * NUMPARTICLES, vertices, GL_DYNAMIC_DRAW);

    /* enable vbo colours */
    glGenBuffers(1, &color_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, color_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GL_FLOAT) * 4 * NUMPARTICLES, 0, GL_DYNAMIC_DRAW);

    /* fill color buffer */
    glBindBuffer(GL_ARRAY_BUFFER, color_vbo);

    /* get write address */
    data = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for(i = 0; i < NUMPARTICLES; i++){

        colorRamp(i/(float)(NUMPARTICLES),data);
        data+=3;
        *data++=1.0f;
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);
    GLERROR;

    /* init water renderer */
    renderer = init_water_renderer(w_width, w_height, NUMPARTICLES);

    GLERROR;
    /* init camera */
    init_camera(w_width, w_height);
    camera_set_pos(-0.07f, 0.06f, -1.12f);
    camera_set_rot(-60.0f,0,60.0f);

    /* Enter GLUT event loop */
    glutMainLoop();

    /* Simulation finished */
    destroy_parts_obj( particles );

    glDeleteBuffers(1, &position_vbo);
    glDeleteBuffers(1, &color_vbo);
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
	static unsigned int enabled_disabled = 1;
	int r_start = (int)(30.0f * (double)rand()/((double)RAND_MAX));
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
        
		if(enabled_disabled){
	        if( (j % 378) == 0 && j > 0 ){
    	        xpos = x_start;
    	        ypos = y_start + (20 + r_start)*parts->default_h/frac_part;
    	    }
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
    
	/* switch method */
	enabled_disabled ^= 1;
	
    return pol;
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
double current_timestep = 0;
static void idle(void)
{
    int i = 0;
    float color = 0.0;
    const double now_T = ftime();
    const double delta_T = now_T - last_T;

    extern double _scale;

    if(delta_T > TIMEWAIT )
    {	
    	current_timestep+=TIMESTEP;
        integrator_step(particles, &polygon);
        
        if(current_timestep >= 3.0f){
        	/* restart simulation */
        	current_timestep = 0.0f;
			polygon = particles_initial_positions( particles );
        }


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

        //printf("%f - %f\n",1/delta_T,_scale);
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
    }else if(key == (int)'m'){
        if(visualization_mode == 1){
            visualization_mode = 0;
        }else{
            visualization_mode = 1;
        }
    }
}

/* draw Boundaries */
void drawBoundaries(void)
{
    int i = 0;
    float * view_matrix;
    float * projection_matrix;

    double fov = 0.0f;
	char buffer[128];
    GLfloat plane_vertices[] = {1.0f, 1.0f, -0.25f - (float) particles->default_h / 2,
                                1.0f, -1.0f, -0.25f - (float) particles->default_h / 2,
                                -1.0f, -1.0f, -0.25f - (float) particles->default_h / 2,
                                -1.0f, 1.0f, -0.25f - (float) particles->default_h / 2};


    GLint indexes[] = {0, 1, 3, 1, 2, 3};
    /* Name-stack manipulation for the purpose of
       selection hit processing when mouse button
       is pressed.  Names are ignored in normal
       OpenGL rendering mode.                    */
    /* Render animation */

    camera_render_view();
    fov = camera_get_fov();


    view_matrix = camera_get_view();
    projection_matrix = camera_get_projection();

    if(visualization_mode == 0) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        //glClearColor( 0.0, 0.0, 0.0, 1.0 );
        /* render plane */
        glDepthMask(GL_TRUE);
        glEnable(GL_DEPTH_TEST);

        glUseProgram(renderer.plane_shader.program);

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
        glDisableVertexAttribArray(0);

        glEnable(GL_POINT_SPRITE);
        glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);

        glUseProgram(depth_program);

        //glUniform1f( glGetUniformLocation(mprogram, "point_scale"), (float)(w_height/(tan(60.0f * 0.5 * PI/180.0))));
        glUniform1f(glGetUniformLocation(depth_program, "point_scale"),
                    (float) (w_height / (tan(fov * 0.5 * PI/180.0))));
        glUniform1f(glGetUniformLocation(depth_program, "point_radius"),
                    (float) (particles->default_h / 4));//0.0125f*0.8f


        /* time to draw the points*/
        glBindBuffer(GL_ARRAY_BUFFER, position_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_t) * particles->num_particles, vertices, GL_DYNAMIC_DRAW);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        glEnableClientState(GL_VERTEX_ARRAY);

        /* Enable colours */
        glBindBuffer(GL_ARRAY_BUFFER, color_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(color_t) * particles->num_particles, colors, GL_DYNAMIC_DRAW);
        glColorPointer(4, GL_FLOAT, 0, 0);
        glEnableClientState(GL_COLOR_ARRAY);

        glDrawArrays(GL_POINTS, 0, particles->num_particles);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);

        glUseProgram(0);
        glDisable(GL_POINT_SPRITE);
        glDisableClientState(GL_VERTEX_ARRAY);

        glLineWidth(1.4);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

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

    }else{
        render_water(renderer, w_width, w_height, fov, (particles->default_h ), particles->num_particles, vertices);
    }

	print(10, 20, "Developed by: Bruno Faria");
	print(10, 46, "Department of Physics");
	print(10, 72, "University of Aveiro");
	sprintf(buffer,"Time instant:%f seconds",current_timestep);
	print(10, 98, buffer);
    glDisable(GL_BLEND);
    glDisable(GL_LINE_SMOOTH);

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


void print(float x, float y, const char *string)
{
	int i = 0;
    //Assume we are in MODEL_VIEW already
	glPushMatrix ();
	glLoadIdentity ();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix ();
	glLoadIdentity();

	GLint viewport [4];
	glGetIntegerv (GL_VIEWPORT, viewport);
	gluOrtho2D (0,viewport[2], viewport[3], 0);

	glDepthFunc (GL_ALWAYS);
	glColor3f (0.888,0.888,0.888);
	glRasterPos2f(x, y);

	for (i = 0; string[i]!= '\0'; ++i)
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, string[i]);

	glDepthFunc (GL_LESS);
	glPopMatrix ();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix ();

}
