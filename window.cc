// open gl
#include "lib/open_gl.hh"
#include "enums.hh"

// tensor files
#include "lib/file_ops/Parameters.hh"
#include <time.h>
#include <fstream>
#include "lib/std/FatalError.hh"
#include "lib/streamlines/Streamline.hh"
#include "lib/volume/VolumeRenderer.hh"
#include "lib/std/AllocArray.hh"
#include "lib/glyphs/VectorFieldRenderer.hh"
#include <iostream>

#define MAX_LIGHTS  8

// callbacks
extern void display();
extern void keyboard(unsigned char key, int x, int y);
extern void mouse(int button, int state, int x, int y);
extern void motion(int x, int y);
extern void mouseMotion(int x, int y);
extern void wheel(int wheel, int direction, int x, int y);
extern void readData();
extern void reshape(int w, int h);
extern void read_spec(char *fname);

extern InteractionMode CurrIM; // from callbacks.cc

extern Streamline			streamlines ; // from callbacks.cc					
extern VolumeRenderer		volRenderer ; // from callbacks.cc					
extern VectorFieldRenderer		vectorField ; // from callbacks.cc
extern Vector ***velocityField; // from callbacks.cc
extern Vector ***eigenValues; // from callbacks.cc
extern Vector 	***eigenVectors; 

extern int slices;
extern int rows;
extern int cols;

extern std::string parameters;

extern int num_lights;

typedef struct LITE{
	GLfloat amb[4];
	GLfloat diff[4];
	GLfloat spec[4];
	GLfloat pos[4];
	GLfloat dir[3];
	GLfloat angle;
}LITE;

extern struct LITE my_lights[MAX_LIGHTS];


float mag = 1.0;

// GLUT specific variables
unsigned int window_width = 512;
unsigned int window_height = 512;

unsigned int timer = 0; // a timer for FPS calculations

// Forward declarations of GL functionality
bool initGL(int argc, char** argv);

void lighting_setup () {
	int i;
	GLfloat globalAmb[]     = {.1, .1, .1, .1};
	
	//enable lighting
//	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);
	
	// reflective propoerites -- global ambiant light
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);
	
	// create flashlight
	GLfloat amb[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat dif[] = {0.8, 0.8, 0.8, 1.0};
	GLfloat spec[] = {4.0, 4.0, 0.0, 1.0};

	GLfloat dir[4];
	dir[2] = -1; 
	
	GLfloat at[] = {0.0,1.0,0.0,0.0};
	
	glLightfv(GL_LIGHT0, GL_POSITION, at);
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);
	
	glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
	glLightfv(GL_LIGHT0, GL_SPECULAR, spec);
	glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 200.0);
	
	glEnable(GL_LIGHT0);
	
	// setup properties of lighting
	for (i=1; i<num_lights; i++) {
		glEnable(GL_LIGHT0+i);
		glLightfv(GL_LIGHT0+i, GL_AMBIENT, my_lights[i].amb);
		glLightfv(GL_LIGHT0+i, GL_DIFFUSE, my_lights[i].diff);
		glLightfv(GL_LIGHT0+i, GL_SPECULAR, my_lights[i].spec);
		glLightfv(GL_LIGHT0+i, GL_POSITION, my_lights[i].pos);
		if ((my_lights[i].dir[0] > 0) ||  (my_lights[i].dir[1] > 0) ||  (my_lights[i].dir[2] > 0)) {
			glLightf(GL_LIGHT0+i, GL_SPOT_CUTOFF, my_lights[i].angle);
			glLightfv(GL_LIGHT0+i, GL_SPOT_DIRECTION, my_lights[i].dir);
		}
	}
	
}

bool createWindow(int argc, char ** argv)
{
	if (false == initGL(argc, argv)) {
		return false;
	}
}

void startApplication(int argc, char ** argv)
{
	readData() ;			// read a vector field of size slices, rows AND cols
	read_spec("/Users/onorinbejasus/Dropbox/Vis/Projects/CombustionProject/Tensor");
	lighting_setup();
	
	double boxLenX = volRenderer.params.CtFileVoxelSpacing.VectorX ;
	double boxLenY = volRenderer.params.CtFileVoxelSpacing.VectorY ;
	double boxLenZ = volRenderer.params.CtFileVoxelSpacing.VectorZ ;
	
	// setup the streamline object
	streamlines.SetVectorField(velocityField, slices, rows, cols) ;
	
	streamlines.Configure( volRenderer.params.SeedPointLoc, 54 ) ; 
	streamlines.SetSpacing(boxLenX, boxLenY, boxLenZ);
	if( streamlines.bReadyToProcess )							// generate streamlines and save them in display lists
		streamlines.CreateStreamlines() ;
	else{
		printf("We have problem. Streamlines can't be generated.\n") ;
		exit(1) ;
	}	
	
	vectorField.Configure( eigenVectors, eigenValues, &streamlines, slices, rows, cols) ;	
	vectorField.SetSpacing(boxLenX, boxLenY, boxLenZ);
		
	// setting up gloabal status variables
	CurrIM = ExploreMode;   // initial interaction mode is explore
	
	vectorField.RenderStreamlines(mag) ;	// render vector field	
	vectorField.bVisible = false;
	
	// start rendering mainloop
	glutMainLoop();

	// clean up
	vectorField.cleanup();
	vectorField.cleanupzoom();
	
	free(velocityField);
	free(eigenValues); 
}

bool initGL(int argc, char **argv)
{
	//Steps 1-2: create a window and GL context (also register callbacks)
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_STENCIL | GLUT_MULTISAMPLE);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("Tensor");
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutPassiveMotionFunc(mouseMotion);
	
	// check for necessary OpenGL extensions
	glewInit();
	if (! glewIsSupported( "GL_VERSION_2_0 " ) ) {
		fprintf(stderr, "ERROR: Support for necessary OpenGL extensions missing.\n");
		return false;
	}
		
  	glClearColor (0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_SMOOTH);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
	//glEnable(GL_BLEND) ;
	glEnable(GL_LINE_SMOOTH) ;
	glEnable(GL_POINT_SMOOTH) ;
	glEnable(GL_MULTISAMPLE);
	
	volRenderer.Initialize(parameters, window_width, window_height) ; 
	
	reshape(window_width, window_height);
			
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity() ;  // init modelview to identity

	return true;
}
