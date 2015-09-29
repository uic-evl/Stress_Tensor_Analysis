/* written by Abedul Haque, December 2009		*/
/* for CS3610 class project						*/

//#include "VBO.h"
#include "VectorFieldRenderer.hh"
#include <cmath>
#include <iostream>
#include <float.h>

extern GLfloat shine;
extern GLfloat emi[4];
extern GLfloat amb[4];
extern GLfloat diff[4];
extern GLfloat spec[4];

/* vector field with eigen values */

void VectorFieldRenderer::Configure( Vector*** vf, Vector ***ev, Streamline *s, int slice, int row, int col ){
	
	this->VectorField = vf ;
	this->EigenValues = ev;

	this->streamline = s;
	this->scale = false;
	this->zoomScale = 1.0;

	this->iSlices = slice;
	this->iRows = row ;
	this->iCols = col ;
	
	this->rendered = false;	
		
	this->skip = 1;
	
}

void VectorFieldRenderer::SetSpacing(double dx, double dy, double dz){
	
	this->boxLenX = dx;
	this->boxLenY = dy;
	this->boxLenZ = dz;
	
}

void VectorFieldRenderer::cleanup(){
	
	glDeleteLists(index, streamline->getSize());
}
void VectorFieldRenderer::cleanupzoom(){
	
	glDeleteLists(scaleIndex, streamline->getSize());	
}

void VectorFieldRenderer::drawQuadric(bool shift){
	
	glListBase(index);	
	if(!shift)
		glCallLists(streamline->getSize()/skip-1, GL_UNSIGNED_INT, lists);
	else
		glCallLists(this->iRows*this->iCols*this->iSlices/125, GL_UNSIGNED_INT, lists);
}
void VectorFieldRenderer::drawZoomedQuadric(){
	
	glListBase(scaleIndex);	
	glCallLists(streamline->getSize()/skip-1, GL_UNSIGNED_INT, scaleList);

}

void VectorFieldRenderer::storeQuadric(GLfloat verts[20][20][4], int rs, int vs){
	
	float U[2][3];
	float V[2][3];
	
	float N[2][3];
	GLfloat finalN[3];
	
	glEnable(GL_LIGHTING);
		
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, &shine);
	
	glBegin(GL_TRIANGLE_FAN);
	
	for(int i = 0; i < rs; i++)
		for(int j = 0; j < vs; j++){
			
			U[0][0] = verts[i][j+1][0]-verts[i][j][0];
			U[0][1] = verts[i][j+1][1]-verts[i][j][1];
			U[0][2] = verts[i][j+1][2]-verts[i][j][2];

			V[0][0] = verts[i+1][j+1][0]-verts[i][j][0];
			V[0][1] = verts[i+1][j+1][1]-verts[i][j][1];
			V[0][2] = verts[i+1][j+1][2]-verts[i][j][2];

			N[0][0] = (U[0][1] * V[0][2]) - (U[0][2]*V[0][1]);
			N[0][1] = (U[0][2] * V[0][0]) - (U[0][0]*V[0][2]);
			N[0][2] = (U[0][0] * V[0][1]) - (U[0][1]*V[0][0]);

			U[1][0] = verts[i][j+1][0]-verts[i+1][j][0];
			U[1][1] = verts[i][j+1][1]-verts[i+1][j][1];
			U[1][2] = verts[i][j+1][2]-verts[i+1][j][2];

			V[1][0] = verts[i+1][j+1][0]-verts[i+1][j][0];
			V[1][1] = verts[i+1][j+1][1]-verts[i+1][j][1];
			V[1][2] = verts[i+1][j+1][2]-verts[i+1][j][2];

			N[1][0] = (U[1][1] * V[1][2]) - (U[1][2]*V[1][1]);
			N[1][1] = (U[1][2] * V[1][0]) - (U[1][0]*V[1][2]);
			N[1][2] = (U[1][0] * V[1][1]) - (U[1][1]*V[1][0]);

			finalN[0] = (N[0][0]+ N[1][0])/2.0;
			finalN[1] = (N[0][1]+ N[1][1])/2.0;
			finalN[2] = (N[0][2]+ N[1][2])/2.0;

		//	normalize(&finalN);
			glNormal3fv(finalN);
			
			glVertex4fv(verts[i][j]);
			glVertex4fv(verts[i][j+1]);
			glVertex4fv(verts[i+1][j+1]);
			glVertex4fv(verts[i+1][j]);
		}

	glEnd();
	
	glDisable(GL_LIGHTING);
	
}

void VectorFieldRenderer::DrawSuperQuadric(Vector e){

	/* draw the superquadrics along the stream lines */

	float denom = e.X() + e.Y() + e.Z(); // common denmoinator
	float cL = ( e.X() - e.Y() )/denom; // linear
	float cP = ( 2.0 * (e.Y() - e.Z() ) )/denom; // planar
	float cS = ( 3.0 * e.Z() )/denom; // spherical

	float alpha = 0.0f;
	float beta = 0.0f;
	float gamma = 2.0f; //anisotropy

	float theta = 0.0f;
	float phi = 0.0f;

	int radialS = 6;
	int verticalS = 6;

	GLfloat verts[20][20][4]; 

	if(cL >= cP){ // if more linear in shape

		alpha = pow( static_cast<float>(1.0-cP), gamma);
		beta = pow( static_cast<float>(1.0-cL), gamma);

		for(int i = 0; i <= radialS; i++){
			theta = (2.0 * MY_PI * i)/radialS;			
			for(int j = 0; j <= verticalS; j++){

				phi = (MY_PI * j)/verticalS;

			verts[i][j][0] = copysign( 1.0, cos(phi) ) * pow( fabs( cos(phi) ), beta);
			verts[i][j][1] = -copysign( 1.0, sin(theta) ) * pow( fabs( sin(theta) ), alpha ) *
						copysign( 1.0, sin(phi) ) * pow( fabs( sin(phi) ), beta );
			verts[i][j][2] = copysign( 1.0, cos(theta) ) * pow( fabs( cos(theta) ), alpha ) *
						copysign( 1.0, sin(phi) ) * pow( fabs( sin(phi) ), beta );
			
			verts[i][j][3] = 1.0;
			
			} // end j

		} // end i

	storeQuadric(verts, radialS, verticalS);

	}else if(cL < cP){
		
		alpha = pow( static_cast<float>(1.0-cL), gamma);
		beta = pow( static_cast<float>(1.0-cP), gamma);

		for(int i = 0; i <= radialS; i++){
			theta = (2.0 * MY_PI * i)/radialS;			
			for(int j = 0; j <= verticalS; j++){

				phi = (MY_PI * j)/verticalS;

				verts[i][j][0] = copysign( 1.0, cos(theta) ) * pow( fabs( cos(theta) ), alpha ) *
						copysign( 1.0, sin(phi) ) * pow( fabs( sin(phi) ), beta );
				verts[i][j][1] = copysign( 1.0, sin(theta) ) * pow( fabs( sin(theta) ), alpha ) *
						copysign( 1.0, sin(phi) ) * pow( fabs( sin(phi) ), beta );
				verts[i][j][2] = copysign( 1.0, cos(phi) ) * pow( fabs( cos(phi) ), beta);
				
				verts[i][j][3] = 1.0;

			} // end j

		} // end i

		storeQuadric(verts, radialS, verticalS);

	} // end else
}

void VectorFieldRenderer::setZoom(GLfloat zoom, float mag, bool redo){
	
	if(this->zoomScale == zoom && !redo) // no need to reconfigure
		return;
			
	this->zoomScale = zoom; // assign new zoom
	scale = true;
		
	this->RenderStreamlines(mag);	
	scale = false;		
}
	 
void VectorFieldRenderer::RenderStreamlines(float mag){

	GLfloat posX, posY, posZ;
	Vector dirV;
	Vector eigV;
	Vector strmpt;

	GLuint locInd;
	GLuint *locList;
	
	int x, y, z;

	GLfloat *verts = (GLfloat*)malloc( (streamline->getSize()) * 4.0 * sizeof(float) + 1);
	
	locInd = glGenLists((int)streamline->getSize()/skip);
	locList = (GLuint *)malloc( sizeof(GLuint)*streamline->getSize()/skip );
		
	int ind = 0;	
		
	for(int i = 0; i < streamline->getSize(); i++){

		if(i%skip != 0)
			continue;
		
		if(ind == (int)(streamline->getSize()/skip - 1) )
			break;

		strmpt = new Vector(streamline->getStreamlinePoint(i)); // get the streamline point
	
		/* get the x, y, and z positions of the streamline point */
		x = verts[i*4+0] = (int)ceil(strmpt.X());
		y = verts[i*4+1] = (int)ceil(strmpt.Y());
		z = verts[i*4+2] = (int)ceil(strmpt.Z());
		verts[i*4+3] = 1.0;

		dirV = this->VectorField[z][y][x]  ;
		eigV = this->EigenValues[z][y][x] ;

		dirV.Normalize() ; 
		posX = x*boxLenX + boxLenX/2 ;
		posY = y*boxLenY + boxLenY/2 ;
		posZ = z*boxLenZ + boxLenZ/2 ;

		glNewList(locInd+ind, GL_COMPILE);
		locList[ind] = ind;
		ind++;

		glPushMatrix() ;
		glTranslatef(posX,posY,posZ) ;
		
		if(scale)
			glScalef(zoomScale, zoomScale, zoomScale);
		glColor3f(1.0, 1.0, 1.0);
//		glColor3f(fabs(posX * eigV.X() * mag), fabs(posY * eigV.Y() * mag), fabs(posZ * eigV.Z() * mag) );
		this->DrawSuperQuadric(eigV);				
		glPopMatrix() ; 
	
		glEndList();
	
	} // end for i
	
	if(scale){ // if scaled
		this->scaleList = locList;
		this->scaleIndex = locInd;
		printf("configure scaled\n");
	}else{
		this->lists = locList;
		this->index = locInd;	
	}
	
	free(verts);	
}

