#include "enums.hh"

// tensor files
#include <time.h>
#include <fstream>

#include "lib/std/FatalError.hh"
#include "lib/streamlines/Streamline.hh"
#include "lib/volume/VolumeRenderer.hh"
#include "lib/std/AllocArray.hh"
#include "lib/glyphs/VectorFieldRenderer.hh"
#include "lib/file_ops/Parameters.hh"
#include "lib/vector/Vector.hh"
#include "lib/vector/Point.hh"
#include "lib/open_gl.hh"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

extern void DrawAxisAndBox(double, double, double);
extern unsigned int window_width; 
extern unsigned int window_height;
extern float mag;

static float zoom = 3.0;
static float radius = 0.1;	

static bool shift = false;

GLfloat vertices_axes[][4] = {
	{0.0, 0.0, 0.0, 1.0},  /* origin */ 
	{5.0, 0.0, 0.0, 1.0},  /* maxx */ 
	{0.0, 5.0, 0.0, 1.0}, /* maxy */ 
	{0.0, 0.0, 5.0, 1.0}  /* maxz */ 

};
GLfloat colors [][3] = {
  {0.0, 0.0, 0.0},  /* black   */
  {1.0, 0.0, 0.0},  /* red     */
  {1.0, 1.0, 0.0},  /* yellow  */
  {1.0, 0.0, 1.0},  /* magenta */
  {0.0, 1.0, 0.0},  /* green   */
  {0.0, 1.0, 1.0},  /* cyan    */
  {0.0, 0.0, 1.0},  /* blue    */
  {0.5, 0.5, 0.5},  /* 50%grey */
  {1.0, 1.0, 1.0}   /* white   */
};

void Magnify();

using namespace std;

InteractionMode CurrIM;

int buttonDownX, buttonDownY;
int LeftButtonDown, RightButtonDown ;
GLfloat anglex, angley, anglez ;

SelectedObjectCatagory ObjectCatagory ;
int iObjectName ;

int iSeedpointMoveDir = -1 ;			// 0 for x, 1 for y and 2 for z

bool magnifying = false;
GLfloat positionX = 0.0;
GLfloat positionY = 0.0;

Streamline				streamlines ;				//streamlines object.
VolumeRenderer		volRenderer ;					//volume renderer
VectorFieldRenderer		vectorField ;				//vector field renderer
Vector 					***velocityField;			//vectorfield data
Vector 					***eigenValues; 
Vector 					***eigenVectors; 

GLfloat 				*eigenfloats; 				// eigen values in float form

int slices = 0;
int rows = 0;
int cols = 0;

void draw_axes( void ) {
  glLineWidth( 5.0 );
  glBegin(GL_LINES); 
  {
    glColor3fv(colors[1]);
    glVertex4fv(vertices_axes[0]);
    glVertex4fv(vertices_axes[1]);
		
    glColor3fv(colors[4]);
    glVertex4fv(vertices_axes[0]);
    glVertex4fv(vertices_axes[2]);
		
    glColor3fv(colors[6]);
    glVertex4fv(vertices_axes[0]);
    glVertex4fv(vertices_axes[3]);
  }
  glEnd();
  glLineWidth( 1.0 );
	
}

// =======================
// = Keyboard Operations =
// =======================

// The keyboard callback
void keyboard(unsigned char key, int x, int y)
{
   switch (key) {
		case 'e':				// explore mode
		case 'E':				// move, rotate, zoom
			printf("Entering EXPLORE_MODE.\n") ;
			CurrIM = ExploreMode ;
			break;
		case 's': // select mode
		case 'S':
			printf("Entering SELECT_MODE.\n") ;
			CurrIM = SelectMode ;
			break ;
		case '1':														// enable/disable volume rendering
			if (volRenderer.bVisible) 
				volRenderer.bVisible = false ;
			else
				volRenderer.bVisible = true ;
			glutPostRedisplay() ;			
			break ;		  
		case '2':													// enable/disable vector field view			
			if( vectorField.bVisible )
				vectorField.bVisible = false ;
			else
				vectorField.bVisible = true ;
			glutPostRedisplay() ;			
			break ;								  
		case '3':													// enable/disable streamlines
			if(streamlines.bVisible)
				streamlines.bVisible  = false ;		
			else
				streamlines.bVisible  = true ;
			glutPostRedisplay() ;			
			break ;
		case 'x':  // toggle grid x
		case 'X':
			if( CurrIM == SelectMode && ObjectCatagory == SEEDPOINT ){
				iSeedpointMoveDir = 0 ;	
			}else{
				vectorField.ToggleSliceVisibilityX();
			}
			glutPostRedisplay();			
			break ;
		case 'y':	// toggle grid y
		case 'Y':
			if( CurrIM == SelectMode && ObjectCatagory == SEEDPOINT ){
				iSeedpointMoveDir = 1 ;	
			}else{
				vectorField.ToggleSliceVisibilityY();
			}
			glutPostRedisplay();
			break ;
		case 'z':		// toggle grid z
		case 'Z':
			if( CurrIM == SelectMode && ObjectCatagory == SEEDPOINT ){
				iSeedpointMoveDir = 2 ;	
			}
			
			else{
				vectorField.ToggleSliceVisibilityZ();
			}
			
			glutPostRedisplay();
			break ;
		case 'm': 		// enable mag glass
		case 'M':
			magnifying = !magnifying;
			glutPostRedisplay();
			break;
		case 'l':		// increase radius
		case 'L':
			radius += 0.01;
			glutPostRedisplay();
			break;
		case 'k':			// decrease radius
		case 'K':
			radius -= 0.01;
			glutPostRedisplay();
			break;
		case '+':  			// zoom in
			if(magnifying){ 
				zoom += 0.5;
			}	
			glutPostRedisplay() ;		
			break ;
		case '-':	// zoom out
			
			if(magnifying){
				zoom -= 0.5;
				if(zoom < 1.0)
					zoom = 1.0;
			}	
			glutPostRedisplay() ;		
			break ;
		
		case '0':
			if(shift){
				vectorField.IncreaseSliceNumber();
				glutPostRedisplay();
			}
		break;
		case '9':
			if(shift){
				vectorField.DecreaseSliceNumber();
			}
		glutPostRedisplay();
		break;
		
		case 27:
			exit(0);
			break ;
			
		case 't':
		case 'T':
			shift = !shift;
			glutPostRedisplay();
			break;
   }

}

// ====================
// = Mouse Operations =
// ====================

void processHits (GLint hits, GLuint buffer[])
{
		
   unsigned int i, j;
   GLuint names, *ptr;

   float fMinMin = 10.0f, fMaxMin= 10.0f ;		// since the z value should be in the range 0 to 1
   int iMinZWinner=-1, iMaxZWinner=-1 ;
   SelectedObjectCatagory ObjCat = NO_OBJECT ;

   //printf ("hits = %d\n", hits);
   ptr = (GLuint *) buffer;

   for (i = 0; i < (unsigned)hits; i++) {	/*  for each hit  */
      names = *ptr;
      //printf(" number of names for hit = %d\n", names); 
	  ptr++;
	  if( names != 2 ){
		printf("We have problem here.\n") ;
		exit(1) ;
	  }
	  float z1 = (float) *ptr/0x7fffffff ; ptr++;
	  float z2 = (float) *ptr/0x7fffffff ; ptr++;
	  int iIndex ;

      //printf(" z1 is %g\n", z1 ); 
      //printf(" z2 is %g\n", z2 ) ; 
		  	
	  int oc = *ptr ; ptr++ ;
	  iIndex = *ptr ; ptr++ ;

	  if( oc == SEEDPOINT ){		// SEEDPOINT HAS HIGHER PRIORITY
		  ObjCat = (SelectedObjectCatagory)oc ;
		  iMinZWinner = iIndex ;
		  printf("Seedpoint selected.\n") ;
	  }								//if( oc == SEEDPOINT )
	  else{
		  if( z1 < fMinMin ){
			fMinMin = z1 ;
			iMinZWinner = iIndex ;
			ObjCat = (SelectedObjectCatagory)oc ;
		  }
		  if( z2 < fMaxMin){
			fMaxMin = z2;
			iMaxZWinner = iIndex ;
		  }
	  }//else
   } //for (i = 0; i < (unsigned)hits; i++) 

	ObjectCatagory =  ObjCat;
	iObjectName  = iMinZWinner ;
	printf("Catagory = %d and ObjecSelected = %d\n",ObjCat, iMinZWinner) ;
}

// selects an object at location X,Y in screen coordinate
// single object selection. Catagory of the object will be set in oc. Right now we have just one catagory. STREAMLINE
// and name of the object will be set in on
#define BUFSIZE					512			// buffer to process mouse click hits
void SelectObjectAtXY(int x, int y){
	GLint hits;
	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);
	GLuint selectBuf[BUFSIZE];

	glSelectBuffer (BUFSIZE, selectBuf);
	(void) glRenderMode (GL_SELECT);

	glInitNames();									// initializing the name stack

	glMatrixMode (GL_PROJECTION);
	glPushMatrix ();
	glLoadIdentity ();
			/*  create 5x5 pixel picking region near cursor location	*/
	gluPickMatrix ((GLdouble) x, (GLdouble) (viewport[3] - y), 
				  5.0, 5.0, viewport);
			
	glMatrixMode (GL_MODELVIEW);
	glPushName(STREAMLINE);							// signature of a streamline		
	streamlines.SelectModeStreamlinesDraw() ; 
	glPopName() ;									// pop streamline signature

	glPushName(SEEDPOINT);							// signature of a seedpoint		
	streamlines.SelectModeSeedPointDraw() ; 
	glPopName() ;									// pop seedpoint signature

	glMatrixMode (GL_PROJECTION);
	glPopMatrix ();
	//glFlush ();

	hits = glRenderMode (GL_RENDER);
	if( hits > 0 )		
		processHits (hits, selectBuf);	
	else{		// no hit
		ObjectCatagory = NO_OBJECT ;
	}
}

void MouseMoveToRotate(int diffX, int diffY){
	// map mousemotion into degrees
	//normalize mouse motion
	float yNormalizer = 360*3/window_width  ;  // rotate 360/3 = 120 degree if your mouse drag covers half of the screen
	float xNormalizer = 360*3/window_height ; // rotate 120 degree if your mouse drag covers half of the screen

	angley += diffX/yNormalizer ;
	anglex += diffY/xNormalizer ;
	
	volRenderer.RotateObject(anglex, angley, anglez ) ;
	volRenderer.SetProjectionResolution(true, true) ;			// Fastpass and intensity 
	
	glutPostRedisplay() ;
}

void MouseMoveToSeedpointMove(float diffX, float diffY){
	fPoint fsp;
	streamlines.GetSeedPoint(fsp) ;
	//printf("%f\n",diffY) ;
	switch(iSeedpointMoveDir){
		case 0:			//along x axis
			fsp.x -= (diffY*2) ;
			break ;
		case 1:
			fsp.y -= (diffY*2) ;	
			break ;
		case 2:
			fsp.z -= (diffY*2) ;
			break ;
	}
	streamlines.MoveSeedPoint(fsp) ; 
} 

void MouseMoveToZoom(int diffX, int diffY){

	volRenderer.ChangeProjectorLocation(diffY) ; 
	volRenderer.SetProjectionResolution(true, true) ;			// Fastpass and intensity 
		
	glutPostRedisplay() ;
}

// mouse motion callback

void motion(int x, int y){	
		
	if( CurrIM == ExploreMode){
		if(LeftButtonDown == 1 ){				// rotate the volume
			MouseMoveToRotate(x-buttonDownX, y-buttonDownY) ;
			buttonDownX = x ;
			buttonDownY = y ;
		}
		else if( RightButtonDown == 1){			// move the camera, zoom in or out
			MouseMoveToZoom(x-buttonDownX, y-buttonDownY) ;
			buttonDownX = x ;
			buttonDownY = y ;	
		}
	} //if( InteractionMode == ExploreMode)
	else if( CurrIM == SelectMode){
		if(LeftButtonDown == 1 ){				// temporary streamline selection or seedpoint move
			if( streamlines.bSeedpointSelected && iSeedpointMoveDir != -1 ){// seedpoint was already selected
				MouseMoveToSeedpointMove(x-buttonDownX, y-buttonDownY) ;
				buttonDownX = x ;
				buttonDownY = y ;
				glutPostRedisplay() ;
			}
			else{
				SelectObjectAtXY(x, y);
				switch(ObjectCatagory){
					case STREAMLINE:
						//printf("ObjectName %d\n",iObjectName) ;
						streamlines.bSeedpointSelected = false ;
						streamlines.TemporarySelectedStreamline(iObjectName) ;
						glutPostRedisplay() ;
						break ;
					case SEEDPOINT:
						if(!streamlines.bSeedpointSelected){		// user just moved in to seed point
							buttonDownX = x ;
							buttonDownY = y ;
							streamlines.bSeedpointSelected = true ;
							glutPostRedisplay() ;
						}
						break ;
					case NO_OBJECT:
						streamlines.TempSelectedStreamline = -1 ;
						streamlines.bSeedpointSelected = false ;
						glutPostRedisplay() ;
						break ;
					default:
						printf("Unidentified object.\n") ;
						break ;
				}  //switch(iObjectCatagory){	
			} // else
		} //if(LeftButtonDown == 1 ){			
	} //else if( CurrIM == SelectMode){
}

// the mouse callback
void mouse(int button, int state, int x ,int y){	

	if( CurrIM == ExploreMode ){				// xplore mode mouse interaction
		if (state == GLUT_DOWN ) { 
			 switch (button) {
			 case GLUT_LEFT_BUTTON:
				  LeftButtonDown = 1 ;
				  buttonDownX = x ;
				  buttonDownY = y ;
				  break;
			 case GLUT_MIDDLE_BUTTON:
				  break;
			 case GLUT_RIGHT_BUTTON:
				  RightButtonDown = 1 ;
				  buttonDownX = x ;
				  buttonDownY = y ;
				  break;
			 }
		} //if (state == GLUT_DOWN ) 
		else if( state == GLUT_UP ){		
			LeftButtonDown = RightButtonDown = 0 ;
			volRenderer.SetProjectionResolution(false, false) ; 
			glutPostRedisplay() ;
		}
    } //if( CurrIM == ExploreMode )
	else if(CurrIM == SelectMode ){		
		if (state == GLUT_DOWN ) { 
			switch (button) {
				case GLUT_LEFT_BUTTON:
					LeftButtonDown = 1 ;
					SelectObjectAtXY(x, y);
					switch(ObjectCatagory){
						case STREAMLINE:
							streamlines.TemporarySelectedStreamline(iObjectName);
							glutPostRedisplay() ;
							break ;
						case SEEDPOINT:
							buttonDownX = x ;
							buttonDownY = y ;
							streamlines.bSeedpointSelected = true ;
							glutPostRedisplay() ;
							break ;
						case NO_OBJECT:
							break;
						default:
							printf("Unidentified object.\n") ;							
							break ;
					}//switch(iObjectCatagory){  			
				break;
			}//switch (button) {
		} //if (state == GLUT_DOWN ) 
		else if( state == GLUT_UP ){		
			LeftButtonDown = RightButtonDown = 0 ;						
			//clear if any streamline was temporarily selected
			streamlines.TempSelectedStreamline = -1 ;
			switch(ObjectCatagory){				
				case STREAMLINE:					
					streamlines.UpdateStreamlineStatus(iObjectName) ;					
					glutPostRedisplay() ;
					break ;
				case SEEDPOINT:
					if( streamlines.bSeedpointMoved )
						streamlines.UpdateSeedPoint() ;
					streamlines.bSeedpointMoved = false ;
					streamlines.bSeedpointSelected = false ;
					ObjectCatagory = NO_OBJECT ;
					iSeedpointMoveDir = -1 ;
					glutPostRedisplay() ;
					break ;
				case NO_OBJECT:
					break ;
			}//switch(ObjectCatagory){  			
		}//else if( state == GLUT_UP )
	} //else if(CurrIM == SelectMode )
}

void Magnify(){
				
	// go back to 2D first	
	glDisable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0); 
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();	
	
//	glPushMatrix();
//	glTranslatef(radius, 0.0, 0.0);
	
	glBegin(GL_POINTS); // draw circle around area
	
	glColor3f(1.0, 1.0, 1.0); // white
	
    for(int i = 0; i <= 360; i++) { // increment thru the degrees
		
        glVertex2f( (float)(positionX/window_width) + radius * cos(i * M_PI/180.0) ,
                   ( (float)(positionY/window_height) ) + radius * sin(i * M_PI/180.0) );
    }
	
    glEnd();
	
//	glPopMatrix();
	
	 glClearStencil(0); // clear stencil
	 	  
	 glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE); // disable drawing
	 glEnable(GL_STENCIL_TEST); // endable stencil test
	 	 
	 glStencilFunc(GL_ALWAYS, 1, 0xFFFFFFFF); // always past test, replace with 1
	 glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE); // always replace
	 		
	// stencil magnification area
			
	glBegin(GL_TRIANGLE_FAN);

    glVertex2f( (float)(positionX/window_width), (float)(positionY/window_height) ); // origin
	
    for(int i = 0; i <= 360; i++) { // increment thru the degrees
		
        glVertex2f( (float)(positionX/window_width) + radius * cos(i * M_PI/180.0) ,
                   ( (float)(positionY/window_height) ) + radius * sin(i * M_PI/180.0) );
    }
	
    glEnd();
	
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE); // enable drawing
	//  
	glStencilFunc(GL_EQUAL, 1, 0xFFFFFFFF); // replace everything within stencil
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP); // keep attributes of drawing
	
	// put back to normal
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glEnable(GL_DEPTH_TEST);		

	if(vectorField.bVisible){
		
		if(shift){ // grid
			vectorField.setZoom(zoom);
		
		}else{	// along streamlines
		
			vectorField.setZoom(zoom, mag, false);
			vectorField.drawZoomedQuadric();
		} // end else
	} // end if visible
		
	glDisable(GL_STENCIL_TEST); // disable stencil test
	
}

// passive motion callback
void mouseMotion(int x, int y){
	
	if(!magnifying)
		return;
		
	/* mouse coordinate */	
	positionX = x;
	positionY = window_height - y;
	
	glutPostRedisplay();	
}

// ====================
// = Display Callback =
// ====================

/// The display callback.
void display()
{	
	glClear(GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW) ;
	glLoadIdentity ();             /* clear the matrix */
      
	/* viewing transformation  */
	gluLookAt(0.0, 0.0, volRenderer.projector.vSourceLocation.VectorZ, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	
	//volRenderer.RenderVolImage() ;			// render the volume image 
	
	glClear(GL_DEPTH_BUFFER_BIT) ;			// so that this volume image does not block other objects
											// Since we have enabled depth test	
		
	glRotatef(anglex,1.0,0.0,0.0) ;
	glRotatef(angley,0.0,1.0,0.0) ;
	glRotatef(anglez,0.0,0.0,1.0) ;
	
	double boxLenX = volRenderer.params.CtFileVoxelSpacing.VectorX ;
	double boxLenY = volRenderer.params.CtFileVoxelSpacing.VectorY ;
	double boxLenZ = volRenderer.params.CtFileVoxelSpacing.VectorZ ;
	
	DrawAxisAndBox(boxLenX, boxLenY, boxLenZ) ;						// draw axis and bounding box
		
	// get the volume center in the coordinate center
	float tz = (slices)*boxLenZ/2 ;
	float ty = (rows)*boxLenY/2 ;
	float tx = (cols)*boxLenX/2 ;		
	
	glTranslatef(-tx, -ty, -tz) ;

	streamlines.RenderStreamlines(eigenfloats) ;		// display streamlines
	
	if(vectorField.bVisible){
		if(shift){
			if(!vectorField.rendered){
				vectorField.RenderSlice();
			}
			vectorField.drawQuadric(shift) ;
		}else
			vectorField.drawQuadric(shift) ;	// draw vector field
	}
	
	if(magnifying){
		Magnify();
		glutSetCursor(GLUT_CURSOR_NONE);
	}
	else
		glutSetCursor(GLUT_CURSOR_TOP_LEFT_CORNER);

	glPopMatrix();

	glutSwapBuffers() ;
}
