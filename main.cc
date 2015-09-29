#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// tensor files
#include <time.h>
#include <fstream>

#include "lib/std/FatalError.hh"
#include "lib/streamlines/Streamline.hh"
#include "lib/std/AllocArray.hh"
#include "lib/glyphs/VectorFieldRenderer.hh"
#include "lib/file_ops/Parameters.hh"
#include "lib/volume/VolumeRenderer.hh"

// open gl
#include "window.hh"
#include "lib/open_gl.hh"

// #include <cuda_runtime.h>
#include <iostream>
#include <cmath>

#define PI 3.14159265

extern VolumeRenderer volRenderer;

std::string parameters;

void reshape(int w, int h) {
		
	glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
	   glMatrixMode (GL_PROJECTION);
	   glLoadIdentity ();
	   float fProjectionPlanePixelSize = volRenderer.projector.fDetectorPixelSize ;

	   glFrustum (	-1.0 * w/2 * fProjectionPlanePixelSize,
					1.0	 * w/2 * fProjectionPlanePixelSize,
					-1.0 * h/2 * fProjectionPlanePixelSize,
					1.0  * h/2 * fProjectionPlanePixelSize, 
					volRenderer.projector.fSourceToDetector, 15000);
					
	   glMatrixMode (GL_MODELVIEW);
}

int main(int argc, char **argv) {

	if (argc < 2){
		printf("Need to specify parameter file\n");
		printf("./main [parameter.csv]\n");
		exit(1);
	}
		
	parameters = argv[1];

	createWindow(argc, argv);
	startApplication(argc, argv);
	return 0;
}


